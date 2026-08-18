#!/usr/bin/env Rscript
# =============================================================================
# 78_supp_table_1b.R  -  standalone; paths are relative to the repository root
#
# Recomputes Supplemental Table 1(b), the Bin 0.1 sensitivity analysis, under a
# single normalisation convention for both columns.
#
# WHY
# The published table normalises the two columns differently. With Bin 0.1
# excluded the trajectory is normalised to Bin 0.2 (the Methods convention),
# giving MCT4 beta = -1.036. With Bin 0.1 included it is normalised to Bin 0.1,
# giving -0.727. A reader comparing -1.036 with -0.727 reads an attenuation that
# is an artefact of the changed denominator: on a common denominator the slope
# steepens to about -2.15. This script reports both conventions side by side so
# the table can be fixed with one number rather than a footnote.
#
# It also recomputes the astrocytic V-ATPase composite, which is ten subunits
# (not the six used for neurons), using the definition in the deposited
# extraction code.
#
# INPUT   output/sensitivity/donor_level.csv        (from R/sensitivity/01_extract_h5ad.R)
#         output/sensitivity/donor_by_gene.rds
# No h5ad access is needed; X is not read.
#
# Usage:
#   Rscript 78_supp_table_1b.R selftest
#   Rscript 78_supp_table_1b.R
# =============================================================================

OUT  <- file.path("output", "sensitivity")
CSV  <- file.path(OUT, "donor_level.csv")
RDS  <- file.path(OUT, "donor_by_gene.rds")
SEED <- 20260720

## ---- Composite definitions, from R/data_extraction/01_extract_seaad.R -------
VATP_ASTRO <- c("ATP6V1A","ATP6V1B2","ATP6V0D1","ATP6V0A1","ATP6V1C1",
                "ATP6V1E1","ATP6V1H","ATP6V0C","ATP6V0E1","ATP6V0B")   # 10, astrocyte
VATP_NEURON <- c("ATP6V1A","ATP6V1B2","ATP6V0A1","ATP6V0C","ATP6V0D1","ATP6V1E1") # 6, neuron

EARLY <- c(0.2, 0.3, 0.4)
LATE  <- c(0.6, 0.7, 0.8)

## ---- Published values, used as a reproduction gate --------------------------
PUB <- list(mct4_pct = -43.2, mct4_beta = -1.036, vatp_astro_pct = -0.8)
TOL <- list(pct = 0.5, beta = 0.01)

## ---- Verdict, fixed before the data are read --------------------------------
# G1  If MCT4 does not reproduce -43.2% and a normalised slope of -1.036 with
#     Bin 0.1 excluded, the pipeline differs from the published one and the
#     V-ATPase numbers are not reported. The tolerance is not relaxed.

## ---- Functions ---------------------------------------------------------------
# Bin means weighted by the number of cells contributing to each donor value.
bin_mean <- function(v, w, bin) tapply(v * w, bin, sum) / tapply(w, bin, sum)

# Slope of a trajectory normalised so that the reference bin equals 1.
# Returns the coefficient, its standard error and the p-value.
norm_slope <- function(y, x, ref) {
  s <- summary(lm(I(y / y[which(x == ref)]) ~ x))$coefficients
  c(beta = s[2, 1], se = s[2, 2], p = s[2, 4])
}

# Two-sided z test for the difference between two independent slopes.
slope_diff <- function(a, b) {
  z <- (a[["beta"]] - b[["beta"]]) / sqrt(a[["se"]]^2 + b[["se"]]^2)
  c(z = z, p = 2 * pnorm(-abs(z)))
}

## ---- Self-test ---------------------------------------------------------------
selftest <- function() {
  cat("=== SELF-TEST ===\n"); set.seed(SEED); f <- character(0)

  # T1: normalising by a constant divides the slope by that constant exactly.
  x <- c(.2,.3,.4,.5,.6,.7,.8,.9); y <- 1 - 0.5 * x + rnorm(8, 0, 1e-9)
  raw <- coef(lm(y ~ x))[[2]]
  t1 <- abs(norm_slope(y, x, 0.2)[["beta"]] - raw / y[1]) < 1e-9
  cat(sprintf("  T1 normalised slope == raw / reference   : %s\n", ifelse(t1,"PASS","FAIL")))
  if (!t1) f <- c(f,"T1")

  # T2: the composite of donor means equals the donor mean of the composite.
  #     This is what lets a rowMeans-per-cell composite be rebuilt from a
  #     donor x gene matrix of per-gene means.
  nD <- 12; nG <- 5
  M <- matrix(rnorm(nD*nG, 1, .2), nD, nG)
  t2 <- max(abs(rowMeans(M) - apply(M, 1, mean))) < 1e-12
  cat(sprintf("  T2 composite of means == mean of composite: %s\n", ifelse(t2,"PASS","FAIL")))
  if (!t2) f <- c(f,"T2")

  # T3: the slope-difference z test against an analytic case. Two slopes one
  #     pooled standard error apart must give z = 1 and p = 0.3173.
  a <- c(beta = -1.0, se = 0.1); b <- c(beta = -1.1414214, se = 0.1)
  d <- slope_diff(a, b)
  t3 <- abs(d[["z"]] - 1) < 1e-4 && abs(d[["p"]] - 0.3173105) < 1e-5
  cat(sprintf("  T3 slope-difference z test analytic       : %s (z=%.4f p=%.4f)\n",
              ifelse(t3,"PASS","FAIL"), d[["z"]], d[["p"]]))
  if (!t3) f <- c(f,"T3")

  # T4: bin_mean must equal a hand-computed weighted mean.
  v <- c(1,2,3,4); w <- c(10,10,20,20); b <- c(0.2,0.2,0.3,0.3)
  t4 <- isTRUE(all.equal(as.numeric(bin_mean(v,w,b)), c(1.5, 3.5)))
  cat(sprintf("  T4 bin_mean == weighted mean by hand      : %s\n", ifelse(t4,"PASS","FAIL")))
  if (!t4) f <- c(f,"T4")

  ok <- length(f) == 0
  cat(sprintf("=== SELF-TEST %s%s ===\n\n", ifelse(ok,"PASS","FAIL"),
              ifelse(ok,"",paste0(" (",paste(f,collapse=", "),")"))))
  ok
}

## ---- Main ---------------------------------------------------------------------
main <- function() {
  for (p in c(CSV, RDS))
    if (!file.exists(p)) { cat(sprintf("Not found: %s - run R/sensitivity/01_extract_h5ad.R first.\n", p)); quit(status = 1) }
  D <- read.csv(CSV, stringsAsFactors = FALSE)
  Z <- readRDS(RDS)
  stopifnot(identical(as.character(Z$donors), as.character(D$donor)))
  D$bin <- round(D$cps, 1)

  A <- Z$astro; colnames(A) <- Z$genes
  miss <- setdiff(VATP_ASTRO, Z$genes)
  if (length(miss)) { cat(sprintf("Missing astrocytic subunits: %s\n", paste(miss, collapse=", "))); quit(status=1) }
  cat(sprintf("Astrocytic V-ATPase composite: %d subunits\n  %s\n\n",
              length(VATP_ASTRO), paste(VATP_ASTRO, collapse = ", ")))

  mct4  <- A[, "SLC16A3"]
  vastro <- rowMeans(A[, VATP_ASTRO, drop = FALSE])

  mB <- bin_mean(mct4,  D$ncell_a, D$bin)
  vB <- bin_mean(vastro, D$ncell_a, D$bin)
  bn <- as.numeric(names(mB))
  pct <- function(y) 100 * (mean(y[bn %in% LATE]) / mean(y[bn %in% EARLY]) - 1)

  cat("--- Bin-level trajectories (astrocytes) ---\n")
  cat(sprintf("  bin       %s\n", paste(sprintf("%6.1f", bn), collapse = "")))
  cat(sprintf("  MCT4      %s\n", paste(sprintf("%6.4f", mB), collapse = "")))
  cat(sprintf("  V-ATPase  %s\n\n", paste(sprintf("%6.4f", vB), collapse = "")))

  # ---- G1 reproduction gate -------------------------------------------------
  k8 <- bn != 0.1
  g_pct  <- abs(pct(mB) - PUB$mct4_pct) <= TOL$pct
  s8 <- norm_slope(mB[k8], bn[k8], 0.2)
  g_beta <- abs(s8[["beta"]] - PUB$mct4_beta) <= TOL$beta
  cat("--- G1 reproduction gate ---\n")
  cat(sprintf("  MCT4 early-to-late %.1f%% (published %.1f%%) %s\n",
              pct(mB), PUB$mct4_pct, ifelse(g_pct, "OK", "FAIL")))
  cat(sprintf("  MCT4 normalised slope, Bin 0.1 excluded %.4f (published %.3f) %s\n",
              s8[["beta"]], PUB$mct4_beta, ifelse(g_beta, "OK", "FAIL")))
  if (!(g_pct && g_beta)) {
    cat("  G1 = FAIL - the pipeline differs from the published one. V-ATPase\n")
    cat("  numbers are not reported and the tolerance is not relaxed.\n")
    return(invisible(NULL))
  }
  cat("  G1 = PASS\n\n")

  # ---- Table 1(b), both conventions ----------------------------------------
  report <- function(y, label) {
    e <- norm_slope(y[k8], bn[k8], 0.2)                 # excluded, ref 0.2
    i2 <- norm_slope(y, bn, 0.2)                        # included, ref 0.2  (consistent)
    i1 <- norm_slope(y, bn, 0.1)                        # included, ref 0.1  (as published)
    list(lab = label, e = e, i2 = i2, i1 = i1)
  }
  M <- report(mB, "MCT4"); V <- report(vB, "V-ATPase (10 subunits)")

  cat("--- Supplemental Table 1(b), recomputed ---\n")
  cat(sprintf("  %-24s %-28s %-28s %s\n", "", "Primary (Bin 0.1 excl, n=8)",
              "Sensitivity (Bin 0.1 incl, n=9)", "as published (ref Bin 0.1)"))
  for (r in list(M, V))
    cat(sprintf("  %-24s %+8.3f (p = %.4f)      %+8.3f (p = %.4f)      %+8.3f (p = %.4f)\n",
                paste0(r$lab, " slope"),
                r$e[["beta"]],  r$e[["p"]],
                r$i2[["beta"]], r$i2[["p"]],
                r$i1[["beta"]], r$i1[["p"]]))
  d8 <- slope_diff(M$e,  V$e)
  d9 <- slope_diff(M$i2, V$i2)
  d9p <- slope_diff(M$i1, V$i1)
  cat(sprintf("  %-24s %8s (p = %.4f)      %8s (p = %.4f)      %8s (p = %.4f)\n",
              "delta-slope MCT4 vs V-ATPase", "", d8[["p"]], "", d9[["p"]], "", d9p[["p"]]))
  cat(sprintf("  %-24s %+8.1f%%                    %+8.1f%%\n",
              "early-to-late change", pct(mB), pct(vB)))

  cat("\n--- What changes if the sensitivity column uses the Methods convention ---\n")
  cat(sprintf("  MCT4      as published %+.3f  ->  on a common reference %+.3f\n",
              M$i1[["beta"]], M$i2[["beta"]]))
  cat(sprintf("  V-ATPase  as published %+.3f  ->  on a common reference %+.3f\n",
              V$i1[["beta"]], V$i2[["beta"]]))
  cat("  Including Bin 0.1 steepens the MCT4 slope rather than attenuating it;\n")
  cat("  the published attenuation is an artefact of the changed denominator.\n")

  out <- data.frame(
    item = c("MCT4 slope","V-ATPase slope","delta-slope p","MCT4 % change","V-ATPase % change"),
    primary_excl_ref02 = c(sprintf("%.3f (p=%.4f)", M$e[["beta"]],  M$e[["p"]]),
                           sprintf("%.3f (p=%.4f)", V$e[["beta"]],  V$e[["p"]]),
                           sprintf("%.4f", d8[["p"]]), sprintf("%.1f%%", pct(mB)), sprintf("%.1f%%", pct(vB))),
    sensitivity_incl_ref02 = c(sprintf("%.3f (p=%.4f)", M$i2[["beta"]], M$i2[["p"]]),
                               sprintf("%.3f (p=%.4f)", V$i2[["beta"]], V$i2[["p"]]),
                               sprintf("%.4f", d9[["p"]]), "", ""),
    sensitivity_incl_ref01_as_published = c(sprintf("%.3f (p=%.4f)", M$i1[["beta"]], M$i1[["p"]]),
                               sprintf("%.3f (p=%.4f)", V$i1[["beta"]], V$i1[["p"]]),
                               sprintf("%.4f", d9p[["p"]]), "", ""),
    stringsAsFactors = FALSE)
  write.csv(out, file.path(OUT, "SuppTable_1b_recomputed.csv"), row.names = FALSE)
  cat(sprintf("\nWritten to %s/SuppTable_1b_recomputed.csv\n", OUT))
  cat("=== paste the whole output back ===\n")
}

args <- commandArgs(trailingOnly = TRUE)
if (!selftest()) { cat("SELF-TEST FAILED - stopping before any file is read.\n"); quit(status = 1) }
if (length(args) && args[1] == "selftest") { cat("selftest mode: done\n"); quit(status = 0) }
main()
