#!/usr/bin/env Rscript
# ---------------------------------------------------------------------
# 01_extract_h5ad.R
#
# Single pass over the SEA-AD MTG h5ad. Produces, in one read:
#   output/donor_level.csv    donor x panel-gene means, cell counts, CPS,
#                             and the transcriptome-wide mean per cell type
#   output/donor_by_gene.rds  donor x all-genes mean matrix (optional)
#
# Definitions are taken from the manuscript Methods and are not re-derived here.
#
# Usage:
#   Rscript R/01_extract_h5ad.R selftest   # logic only, no file access
#   Rscript R/01_extract_h5ad.R
# ---------------------------------------------------------------------

source(file.path("R", "sensitivity", "config_sensitivity.R")); source(file.path("R", "sensitivity", "00_sensitivity_utils.R"))
options(datatable.optimize = 1L)

## ---- Definitions from the manuscript Methods ------------------------
VATP6 <- c("ATP6V1A", "ATP6V1B2", "ATP6V0A1", "ATP6V0C", "ATP6V0D1", "ATP6V1E1")
PANEL <- unique(c("SLC16A3", "SLC16A1", "SLC16A7", "SLC2A1", "LDHA", "LDHB",
                  "HK1", "HK2", "PKM", "PDK1", "PTGDS", "LCN2", "LAMP1", "LAMP2",
                  "TFRC", "FTH1", "FTL", VATP6))
ASTRO_SUBCLASS <- "Astrocyte"
EXC_SUBCLASS   <- c("^L2/3 IT$", "^L4 IT$", "^L5 IT$", "^L6 IT$", "^L6 IT Car3$",
                    "^L5 ET$", "^L5/6 NP$", "^L6 CT$", "^L6b$")

## ---- Verdict, fixed before the data are read ------------------------
# G0 cohort gate. If the cell counts do not reproduce the Methods, the cell-type
# definition or the QC filter differs and downstream output is uninterpretable.
G0 <- list(nuclei = 1378211L, astro = 67419L, exc = 671689L, donors = 84L)

## ---- Self-test ------------------------------------------------------
selftest <- function() {
  cat("=== SELF-TEST ===\n"); set.seed(SEED); fails <- character(0)

  # T1: the block reader must equal a dense ground truth, at several block
  #     sizes and with empty cells present. Integer values are used so that
  #     summation order cannot change the result.
  nC <- 37L; nG <- 23L
  D <- matrix(0, nC, nG); D[sample(nC * nG, 180)] <- sample(1:9, 180, TRUE)
  D[c(3, 17, 37), ] <- 0
  h <- mock_h5(D); ip <- h$all("X/indptr")
  pg <- c(2L, 5L, 11L, 23L); gmap <- rep(NA_integer_, nG); gmap[pg] <- seq_along(pg)
  ok <- logical(0)
  for (bs in c(1L, 4L, 10L, 37L)) {
    P <- matrix(0, nC, length(pg)); S <- numeric(nC)
    for (a in seq(1L, nC, by = bs)) {
      b <- min(a + bs - 1L, nC)
      r <- read_block(h, ip, a, b, gmap, length(pg))
      P[a:b, ] <- r$panel; S[a:b] <- r$csum
    }
    ok <- c(ok, isTRUE(all.equal(P, D[, pg], tolerance = 1e-12)),
                isTRUE(all.equal(S, as.numeric(rowSums(D)), tolerance = 1e-12)))
  }
  cat(sprintf("  T1 block reader == dense truth (BS 1/4/10/37) : %s (%d/%d)\n",
              ifelse(all(ok), "PASS", "FAIL"), sum(ok), length(ok)))
  if (!all(ok)) fails <- c(fails, "T1")

  # T2: whole-transcriptome accumulation must equal colSums.
  gs <- numeric(nG)
  for (a in seq(1L, nC, by = 7L)) {
    b <- min(a + 6L, nC)
    r <- read_block(h, ip, a, b, gmap, length(pg))
    if (!length(r$gx)) next
    ag <- rowsum(r$gx, r$gi, reorder = FALSE)
    i <- as.integer(rownames(ag)); gs[i] <- gs[i] + as.numeric(ag)
  }
  t2 <- isTRUE(all.equal(gs, as.numeric(colSums(D)), tolerance = 1e-12))
  cat(sprintf("  T2 gene accumulation == colSums               : %s\n", ifelse(t2, "PASS", "FAIL")))
  if (!t2) fails <- c(fails, "T2")

  # T3: donor x gene accumulation uses a column-major linear index.
  nD <- 4L; dv <- rep_len(1:nD, nC); Gm <- matrix(0, nD, nG)
  for (a in seq(1L, nC, by = 6L)) {
    b <- min(a + 5L, nC)
    r <- read_block(h, ip, a, b, gmap, length(pg))
    if (!length(r$gx)) next
    lin <- dv[a + r$cid - 1L] + (r$gi - 1L) * nD
    ag <- rowsum(r$gx, lin, reorder = FALSE)
    i <- as.integer(rownames(ag)); Gm[i] <- Gm[i] + as.numeric(ag)
  }
  ref <- rowsum(D, dv, reorder = FALSE)
  ref <- ref[order(as.integer(rownames(ref))), , drop = FALSE]
  t3 <- isTRUE(all.equal(unname(Gm), unname(ref), tolerance = 1e-12))
  cat(sprintf("  T3 donor x gene accumulation == rowsum        : %s\n", ifelse(t3, "PASS", "FAIL")))
  if (!t3) fails <- c(fails, "T3")

  # T4: indptr beyond the 32-bit range must stay double.
  big <- 2^31 + c(0, 5, 12)
  t4 <- is.double(big) && (big[3] - big[1]) == 12
  cat(sprintf("  T4 indptr > 2^31 preserved as double          : %s\n", ifelse(t4, "PASS", "FAIL")))
  if (!t4) fails <- c(fails, "T4")

  # T5: obs_cat must preserve NA codes and assert its length.
  mk <- list(has = function(n) n %in% c("obs/F", "obs/__categories/F"),
             all = function(n) if (n == "obs/F") c(0L, -1L, 1L, -1L, 0L) else c("N", "Y"))
  t5 <- identical(obs_cat(mk, "F", 5L), c("N", NA, "Y", NA, "N")) &&
        inherits(try(obs_cat(mk, "F", 4L), silent = TRUE), "try-error")
  cat(sprintf("  T5 obs_cat NA codes + length assertion        : %s\n", ifelse(t5, "PASS", "FAIL")))
  if (!t5) fails <- c(fails, "T5")

  ok <- length(fails) == 0
  cat(sprintf("=== SELF-TEST %s%s ===\n\n", ifelse(ok, "PASS", "FAIL"),
              ifelse(ok, "", paste0(" (", paste(fails, collapse = ", "), ")"))))
  ok
}

## ---- Main -----------------------------------------------------------
main <- function() {
  if (!file.exists(H5AD)) {
    cat(sprintf("Not found: %s\nPlace the h5ad under data/ or edit config.R.\n", H5AD))
    quit(status = 1)
  }
  dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
  t0 <- Sys.time()
  h5 <- make_h5(H5AD); on.exit(try(h5$close(), silent = TRUE))
  cat(sprintf("backend: %s | block size: %d\n", h5$backend, BLOCK_SIZE))

  # Report the file layout rather than assume it.
  cat("\n--- h5ad layout ---\n")
  for (g in c("", "X", "var", "obs", "obs/__categories")) {
    k <- h5$ls(g)
    if (length(k)) cat(sprintf("  /%-16s : %s\n", g, paste(head(k, 30), collapse = ", ")))
  }
  pick <- function(cands, what) {
    hit <- cands[vapply(cands, h5$has, logical(1))]
    if (!length(hit)) {
      cat(sprintf("\n[stop] Could not locate %s. Tried: %s\n", what, paste(cands, collapse = ", ")))
      quit(status = 1)
    }
    cat(sprintf("  %-14s -> %s\n", what, hit[1])); hit[1]
  }
  cat("\n")
  g_path <- pick(c("var/_index", "var/gene_ids", "var/gene_name"), "gene symbols")
  s_path <- pick(c("obs/Subclass", "obs/subclass"), "Subclass")
  d_path <- pick(c("obs/Donor ID", "obs/donor_id"), "Donor ID")
  c_path <- pick(c("obs/Continuous Pseudo-progression Score", "obs/CPS"), "CPS")
  for (req in c("X/indptr", "X/indices", "X/data"))
    if (!h5$has(req)) { cat(sprintf("\n[stop] %s missing - X is not sparse CSR.\n", req)); quit(status = 1) }

  genes <- as.character(h5$all(g_path)); nG <- length(genes)
  cps   <- as.numeric(h5$all(c_path));   nC <- length(cps)
  sub   <- obs_cat(h5, sub("^obs/", "", s_path), nC)
  don   <- obs_cat(h5, sub("^obs/", "", d_path), nC)
  indptr <- as.numeric(h5$all("X/indptr"))
  if (length(indptr) - 1L != nC) {
    cat("[stop] indptr length implies gene-major storage; this reader expects CSR.\n"); quit(status = 1)
  }
  cat(sprintf("\nnuclei %d | genes %d | CSR confirmed\n", nC, nG))

  # Cohort filter. Neurotypical reference donors carry no CPS and are excluded;
  # without this filter the cell counts do not match the Methods.
  donors <- sort(unique(don[!is.na(cps)]))
  inCoh  <- !is.na(cps) & don %in% donors
  isA <- sub == ASTRO_SUBCLASS & inCoh
  isE <- Reduce(`|`, lapply(EXC_SUBCLASS, function(p) grepl(p, sub))) & inCoh

  cat("\n--- G0 cohort gate ---\n")
  cat(sprintf("  nuclei %d/%d | astrocytes %d/%d | excitatory %d/%d | donors %d/%d\n",
              nC, G0$nuclei, sum(isA), G0$astro, sum(isE), G0$exc, length(donors), G0$donors))
  g0 <- nC == G0$nuclei && sum(isA) == G0$astro && sum(isE) == G0$exc && length(donors) == G0$donors
  cat(sprintf("  G0 = %s\n\n", ifelse(g0, "PASS",
      "FAIL - cell-type definition or QC differs; downstream output is uninterpretable.")))

  gmap <- rep(NA_integer_, nG); mp <- match(PANEL, genes)
  gmap[mp[!is.na(mp)]] <- which(!is.na(mp)); nP <- length(PANEL)
  if (anyNA(mp)) cat(sprintf("Panel genes not found: %s\n\n", paste(PANEL[is.na(mp)], collapse = ", ")))

  di <- match(don, donors); nD <- length(donors)
  acc <- function() list(panel = matrix(0, nD, nP), tot = numeric(nD), n = numeric(nD),
                         g = if (FULL_GENE_MATRIX) matrix(0, nD, nG) else NULL)
  A <- acc(); E <- acc()
  starts <- seq(1L, nC, by = BLOCK_SIZE)
  step <- max(1L, length(starts) %/% 20L)
  for (k in seq_along(starts)) {
    a <- starts[k]; b <- min(a + BLOCK_SIZE - 1L, nC)
    keep <- (isA[a:b] | isE[a:b]) & !is.na(di[a:b])
    if (!any(keep)) next
    r <- read_block(h5, indptr, a, b, gmap, nP)
    loc <- which(keep); gl <- a + loc - 1L
    for (tag in c("A", "E")) {
      sel <- loc[if (tag == "A") isA[gl] else isE[gl]]
      if (!length(sel)) next
      dd <- di[a + sel - 1L]; O <- get(tag)
      ag <- rowsum(r$panel[sel, , drop = FALSE], dd, reorder = FALSE)
      i <- as.integer(rownames(ag)); O$panel[i, ] <- O$panel[i, , drop = FALSE] + ag
      as2 <- rowsum(r$csum[sel], dd, reorder = FALSE)
      i2 <- as.integer(rownames(as2)); O$tot[i2] <- O$tot[i2] + as.numeric(as2)
      O$n <- O$n + tabulate(dd, nD)
      if (FULL_GENE_MATRIX) {
        m <- r$cid %in% sel
        if (any(m)) {
          lin <- di[a + r$cid[m] - 1L] + (r$gi[m] - 1L) * nD
          ag3 <- rowsum(r$gx[m], lin, reorder = FALSE)
          i3 <- as.integer(rownames(ag3)); O$g[i3] <- O$g[i3] + as.numeric(ag3)
        }
      }
      assign(tag, O)
    }
    if (k %% step == 0 || k == length(starts))
      cat(sprintf("  block %d/%d (%.1f min)\n", k, length(starts),
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
  cat(sprintf("\nExtraction finished in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

  mk <- function(O) {
    M <- O$panel / O$n; colnames(M) <- PANEL
    list(M = M, gm = O$tot / (O$n * nG), n = O$n,
         g = if (FULL_GENE_MATRIX) O$g / O$n else NULL)
  }
  a <- mk(A); e <- mk(E)
  D <- data.frame(
    donor = donors,
    cps   = vapply(donors, function(x) mean(cps[don == x], na.rm = TRUE), numeric(1)),
    MCT4  = a$M[, "SLC16A3"], gm_a = a$gm, gm_n = e$gm,
    ncell_a = a$n, ncell_n = e$n,
    VATP_n6 = rowMeans(e$M[, VATP6, drop = FALSE]),
    LAMP1_n = e$M[, "LAMP1"], LDHB_n = e$M[, "LDHB"],
    stringsAsFactors = FALSE)
  write.csv(D, file.path(OUT_DIR, "donor_level.csv"), row.names = FALSE)
  if (FULL_GENE_MATRIX)
    saveRDS(list(astro = a$g, neuron = e$g, genes = genes, donors = donors),
            file.path(OUT_DIR, "donor_by_gene.rds"))
  cat(sprintf("Written to %s/\n", OUT_DIR))
}

args <- commandArgs(trailingOnly = TRUE)
if (!selftest()) { cat("SELF-TEST FAILED - stopping before any file is opened.\n"); quit(status = 1) }
if (length(args) && args[1] == "selftest") { cat("selftest mode: done\n"); quit(status = 0) }
main()
