#!/usr/bin/env Rscript
# =============================================================================
# 04_detection_matched_null.R
#
# WHY THIS SCRIPT EXISTS
# The detection-matched null for SLC16A3 currently has four different values in
# this repository and in the manuscript:
#
#   source                                detection   n matched   percentile
#   data/38_all_genes.csv (recomputed)      0.066        2,090        0.33
#   data/38_all_genes.csv (stored column)   0.066          -          0.70
#   output/tables/SuppTable_detection_null  0.0669       2,117        0.43
#   output/tables/detection_null.csv        0.0612       2,042        0.39
#
# The Figure 6 caption cites the third; Figure 6C renders the first, because the
# figure script recomputes the band at render time from data/38_all_genes.csv.
# The Limitations section cites 6.1%, which is the whole-astrocyte detection
# rate used in the microglia comparison, not the bin-restricted rate used for
# the null. Both are astrocytic SLC16A3 detection rates over different cell sets.
#
# This script computes the quantity once, from raw data, with every filter
# stated explicitly, so that a single set of numbers can be propagated to the
# body, the Figure 6 caption, Figure 6C and the supplementary table.
#
# Usage:
#   Rscript R/audit/04_detection_matched_null.R selftest
#   Rscript R/audit/04_detection_matched_null.R
# =============================================================================

source(file.path("R", "sensitivity", "config_sensitivity.R"))
source(file.path("R", "sensitivity", "00_sensitivity_utils.R"))

## ---- Filters, declared explicitly ------------------------------------
# Two detection rates are produced, because the manuscript uses both:
#   det_all  - over every astrocyte in the 84-donor cohort
#              (this is the number quoted in the microglia comparison)
#   det_bin  - over astrocytes in Bins 0.2-0.9 only, i.e. the same cells that
#              enter the early-to-late change; this is the correct denominator
#              for a null matched to that change
BIN_EXCLUDE   <- 0.1          # Bin 0.1 is excluded from bin-level analyses
EARLY_BINS    <- c(0.2, 0.3, 0.4)
LATE_BINS     <- c(0.6, 0.7, 0.8)
BAND_TOL      <- 0.25         # +/- 25% of the target detection rate
MIN_EXPRESSED <- 0            # genes with zero detection are dropped

selftest <- function() {
  cat("=== SELF-TEST ===\n"); set.seed(SEED); f <- character(0)

  # Detection counting must equal a dense ground truth, per gene.
  nC <- 41L; nG <- 19L
  D <- matrix(0, nC, nG); D[sample(nC * nG, 150)] <- sample(1:9, 150, TRUE)
  h <- mock_h5(D); ip <- h$all("X/indptr")
  det <- numeric(nG)
  for (a in seq(1L, nC, by = 5L)) {
    b <- min(a + 4L, nC)
    r <- read_block(h, ip, a, b, rep(NA_integer_, nG), 0L)
    if (!length(r$gx)) next
    tb <- tabulate(r$gi, nG); det <- det + tb
  }
  t1 <- identical(as.numeric(det), as.numeric(colSums(D != 0)))
  cat(sprintf("  T1 detection count == colSums(D != 0)   : %s\n", ifelse(t1, "PASS", "FAIL")))
  if (!t1) f <- c(f, "T1")

  # Band selection and percentile must reproduce an analytic answer.
  dt <- c(0.10, 0.08, 0.12, 0.05, 0.20); pc_ <- c(-50, -10, 0, -90, +5)
  idx <- which(dt >= 0.10 * (1 - BAND_TOL) & dt <= 0.10 * (1 + BAND_TOL))
  t2 <- setequal(idx, c(1, 2, 3)) && abs(mean(pc_[setdiff(idx, 1)] < pc_[1]) - 0) < 1e-9
  cat(sprintf("  T2 band selection + percentile analytic : %s\n", ifelse(t2, "PASS", "FAIL")))
  if (!t2) f <- c(f, "T2")

  mk <- list(has = function(n) n %in% c("obs/F", "obs/__categories/F"),
             all = function(n) if (n == "obs/F") c(0L, -1L, 1L) else c("N", "Y"))
  t3 <- identical(obs_cat(mk, "F", 3L), c("N", NA, "Y"))
  cat(sprintf("  T3 obs_cat preserves NA codes           : %s\n", ifelse(t3, "PASS", "FAIL")))
  if (!t3) f <- c(f, "T3")

  ok <- length(f) == 0
  cat(sprintf("=== SELF-TEST %s%s ===\n\n", ifelse(ok, "PASS", "FAIL"),
              ifelse(ok, "", paste0(" (", paste(f, collapse = ", "), ")"))))
  ok
}

main <- function() {
  if (!file.exists(H5AD)) { cat(sprintf("Not found: %s\n", H5AD)); quit(status = 1) }
  dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
  h5 <- make_h5(H5AD); on.exit(try(h5$close(), silent = TRUE))

  genes <- as.character(h5$all("var/_index")); nG <- length(genes)
  cps   <- as.numeric(h5$all("obs/Continuous Pseudo-progression Score")); nC <- length(cps)
  sub   <- obs_cat(h5, "Subclass", nC)
  don   <- obs_cat(h5, "Donor ID", nC)
  donors <- sort(unique(don[!is.na(cps)]))
  isA <- sub == "Astrocyte" & !is.na(cps) & don %in% donors
  bin <- round(cps, 1)
  cat(sprintf("Astrocytes in cohort: %d | of which in bins %.1f-0.9: %d\n",
              sum(isA), min(EARLY_BINS), sum(isA & bin != BIN_EXCLUDE)))

  indptr <- as.numeric(h5$all("X/indptr"))
  det_all <- numeric(nG); det_bin <- numeric(nG)
  n_all <- 0L; n_bin <- 0L
  sum_by_bin <- matrix(0, nG, 9); n_by_bin <- numeric(9)   # bins 0.1 .. 0.9
  starts <- seq(1L, nC, by = BLOCK_SIZE)
  step <- max(1L, length(starts) %/% 20L)
  for (k in seq_along(starts)) {
    a <- starts[k]; b <- min(a + BLOCK_SIZE - 1L, nC)
    sel <- which(isA[a:b])
    if (!length(sel)) next
    r <- read_block(h5, indptr, a, b, rep(NA_integer_, nG), 0L)
    if (!length(r$gx)) next
    m <- r$cid %in% sel
    if (!any(m)) next
    gi <- r$gi[m]; gx <- r$gx[m]; ci <- a + r$cid[m] - 1L
    det_all <- det_all + tabulate(gi, nG); n_all <- n_all + length(sel)
    inb <- bin[ci] != BIN_EXCLUDE
    if (any(inb)) {
      det_bin <- det_bin + tabulate(gi[inb], nG)
      n_bin <- n_bin + sum(bin[a + sel - 1L] != BIN_EXCLUDE)
    }
    bi <- round(bin[ci] * 10)
    ok <- bi >= 1 & bi <= 9
    if (any(ok)) {
      lin <- gi[ok] + (bi[ok] - 1L) * nG
      ag <- rowsum(gx[ok], lin, reorder = FALSE)
      i <- as.integer(rownames(ag)); sum_by_bin[i] <- sum_by_bin[i] + as.numeric(ag)
    }
    bs <- round(bin[a + sel - 1L] * 10)
    n_by_bin <- n_by_bin + tabulate(bs[bs >= 1 & bs <= 9], 9)
    if (k %% step == 0 || k == length(starts))
      cat(sprintf("  block %d/%d\n", k, length(starts)))
  }
  det_all <- det_all / n_all
  det_bin <- det_bin / n_bin
  mean_by_bin <- sweep(sum_by_bin, 2, pmax(n_by_bin, 1), "/")
  early <- rowMeans(mean_by_bin[, round(EARLY_BINS * 10), drop = FALSE])
  late  <- rowMeans(mean_by_bin[, round(LATE_BINS  * 10), drop = FALSE])
  pct <- 100 * (late / early - 1)

  keep <- which(det_bin > MIN_EXPRESSED & is.finite(pct))
  i4 <- match("SLC16A3", genes)
  d0 <- det_bin[i4]
  band <- keep[det_bin[keep] >= d0 * (1 - BAND_TOL) & det_bin[keep] <= d0 * (1 + BAND_TOL)]
  band <- setdiff(band, i4)
  pctl <- 100 * mean(pct[band] < pct[i4], na.rm = TRUE)

  cat("\n=== Authoritative detection-matched null ===\n")
  cat(sprintf("  SLC16A3 detection, all astrocytes      : %.4f  (%.1f%%)\n", det_all[i4], 100 * det_all[i4]))
  cat(sprintf("  SLC16A3 detection, bins %.1f-0.9        : %.4f  (%.1f%%)\n",
              min(EARLY_BINS), det_bin[i4], 100 * det_bin[i4]))
  cat(sprintf("  early-to-late change                   : %.1f%%\n", pct[i4]))
  cat(sprintf("  matched band (+/-%.0f%%)                  : %.4f - %.4f\n",
              BAND_TOL * 100, d0 * (1 - BAND_TOL), d0 * (1 + BAND_TOL)))
  cat(sprintf("  genes in band                          : %d\n", length(band)))
  cat(sprintf("  SLC16A3 percentile within the band     : %.2f\n", pctl))
  cat("\n  Use the bin-restricted rate for the null and the all-astrocyte rate\n")
  cat("  for the microglia comparison, and say which is which in both places.\n")

  out <- data.frame(gene = genes, detect_all = det_all, detect_bin = det_bin,
                    pct_change = pct, stringsAsFactors = FALSE)
  write.csv(out[keep, ], file.path(OUT_DIR, "detection_matched_null_all_genes.csv"), row.names = FALSE)
  write.csv(data.frame(gene = "SLC16A3", detect_all = det_all[i4], detect_bin = d0,
                       pct_change = pct[i4], n_band = length(band), percentile = pctl),
            file.path(OUT_DIR, "detection_matched_null_summary.csv"), row.names = FALSE)
  cat(sprintf("\nWritten to %s/\n", OUT_DIR))
}

args <- commandArgs(trailingOnly = TRUE)
if (!selftest()) { cat("SELF-TEST FAILED - stopping before any file is opened.\n"); quit(status = 1) }
if (length(args) && args[1] == "selftest") { cat("selftest mode: done\n"); quit(status = 0) }
main()
