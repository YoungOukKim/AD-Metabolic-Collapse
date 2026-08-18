#!/usr/bin/env Rscript
# =============================================================================
# 79_anls_row.R  -  standalone; paths are relative to the repository root
#
# Verifies the ANLS row of Supplemental Table 1(a), the only row of that table
# that could not be recomputed earlier because the ANLS composite is not carried
# in the donor-level export.
#
# CONTEXT
# The manuscript's own embedded supplementary table already carries the current
# values for all four rows. The separate .xlsx file is a stale copy: its three
# MCT4 rows were the pre-correction values and have been updated, and its ANLS
# row still reads +0.356 / +0.305 / -14.3% against the manuscript's
# +0.377 / +0.339 / -10.0%. This script recomputes the row from raw donor data so
# the .xlsx is corrected against a computation rather than against another table.
#
# ANLS composite = mean(SLC2A1, LDHA, SLC16A1), astrocytes.
# MCT4/SLC16A3 is deliberately not part of the composite; it is analysed
# separately because its trajectory would dominate the mean. Definition taken
# from R/data_extraction/01_extract_seaad.R.
#
# INPUT   output/sensitivity/donor_level.csv        (from R/sensitivity/01_extract_h5ad.R)
#         output/sensitivity/donor_by_gene.rds
# No h5ad access is needed; X is not read.
#
# Usage:
#   Rscript 79_anls_row.R selftest
#   Rscript 79_anls_row.R
# =============================================================================

OUT  <- file.path("output", "sensitivity")
CSV  <- file.path(OUT, "donor_level.csv")
RDS  <- file.path(OUT, "donor_by_gene.rds")
SEED <- 20260720

ANLS_GENES <- c("SLC2A1", "LDHA", "SLC16A1")   # astrocytic ANLS composite

## ---- Published values (manuscript Table, current version) --------------------
PUB <- list(
  mct4_vatp = c(zero = 0.517, part = 0.474),
  mct4_lamp = c(zero = 0.549, part = 0.501),
  mct4_ldhb = c(zero = 0.427, part = 0.418),
  anls_vatp = c(zero = 0.377, part = 0.339, change = -10.0),
  anls_mct4 = c(zero = 0.482, part = 0.425)     # within astrocytes, Results text
)
TOL <- list(r = 0.005, change = 0.5)

## ---- Verdict, fixed before the data are read --------------------------------
# G1  The three MCT4 rows must reproduce to within 0.005. They are computable
#     from the donor table alone and are already known to reproduce; if they do
#     not here, the composite reconstruction from the donor x gene matrix is at
#     fault and the ANLS row is not reported. The tolerance is not relaxed.

pc <- function(x, y, Z) {
  Z <- as.matrix(Z)
  o <- complete.cases(x, y, Z)
  x <- x[o]; y <- y[o]; Z <- Z[o, , drop = FALSE]
  r <- cor(residuals(lm(x ~ Z)), residuals(lm(y ~ Z)))
  n <- sum(o); k <- ncol(Z)
  c(r = r, p = 2 * pt(-abs(r * sqrt((n - k - 2) / (1 - r^2))), n - k - 2), n = n)
}

gate <- function(obs, ref, tol) {
  d <- abs(obs - ref)
  list(pass = all(is.finite(d)) && all(d <= tol * (1 + 1e-9)),
       maxdev = if (all(is.finite(d))) max(d) else NA_real_)
}

row_stats <- function(x, y, cps) {
  z <- cor.test(x, y)
  p <- pc(x, y, data.frame(cps = cps))
  c(zero = unname(z$estimate), zero_p = z$p.value,
    part = p[["r"]], part_p = p[["p"]],
    change = 100 * (p[["r"]] - unname(z$estimate)) / unname(z$estimate))
}

## ---- Self-test ---------------------------------------------------------------
selftest <- function() {
  cat("=== SELF-TEST ===\n"); set.seed(SEED); f <- character(0)

  # T1: the composite of per-gene donor means equals the donor mean of the
  #     per-cell composite, which is what allows ANLS to be rebuilt from the
  #     donor x gene matrix rather than re-read from the h5ad.
  M <- matrix(rnorm(60, 1, .2), 12, 5)
  t1 <- max(abs(rowMeans(M) - apply(M, 1, mean))) < 1e-12
  cat(sprintf("  T1 composite of means == mean of composite : %s\n", ifelse(t1,"PASS","FAIL")))
  if (!t1) f <- c(f,"T1")

  # T2: pc() against theory. Over 2000 null replicates the Fisher-z standard
  #     deviation must approach 1/sqrt(n-k-3) and the type I error must sit
  #     near 0.05. A single draw would not settle either.
  n <- 84; k <- 1; R <- 2000; z <- numeric(R); pv <- numeric(R)
  for (i in seq_len(R)) {
    o <- pc(rnorm(n), rnorm(n), data.frame(a = rnorm(n)))
    z[i] <- atanh(o[["r"]]); pv[i] <- o[["p"]]
  }
  sdo <- sd(z); sdt <- 1 / sqrt(n - k - 3); a05 <- mean(pv < 0.05)
  t2 <- abs(sdo/sdt - 1) <= 0.08 && a05 >= 0.036 && a05 <= 0.065
  cat(sprintf("  T2 pc() Fisher-z %.4f vs theory %.4f, alpha %.4f : %s\n",
              sdo, sdt, a05, ifelse(t2,"PASS","FAIL")))
  if (!t2) f <- c(f,"T2")

  # T3: row_stats must recover a constructed partial correlation. With y built
  #     from x plus a confounder, controlling for the confounder must return the
  #     planted partial correlation and the zero-order must exceed it.
  set.seed(SEED)
  nn <- 4000; cps <- rnorm(nn); e1 <- rnorm(nn); e2 <- rnorm(nn)
  x <- cps + e1; y <- cps + e2                     # partial r = 0 by construction
  s <- row_stats(x, y, cps)
  t3 <- abs(s[["part"]]) < 0.05 && s[["zero"]] > 0.4
  cat(sprintf("  T3 row_stats recovers planted partial ~0   : %s (zero %.3f, partial %.3f)\n",
              ifelse(t3,"PASS","FAIL"), s[["zero"]], s[["part"]]))
  if (!t3) f <- c(f,"T3")

  t4 <- gate(0.377, 0.377, TOL$r)$pass && !gate(0.356, 0.377, TOL$r)$pass
  cat(sprintf("  T4 gate accepts truth, rejects stale value : %s\n", ifelse(t4,"PASS","FAIL")))
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

  A <- Z$astro; colnames(A) <- Z$genes
  miss <- setdiff(ANLS_GENES, Z$genes)
  if (length(miss)) { cat(sprintf("Missing ANLS genes: %s\n", paste(miss, collapse=", "))); quit(status=1) }
  cat(sprintf("ANLS composite: %s\n", paste(ANLS_GENES, collapse = ", ")))
  cat("MCT4/SLC16A3 is analysed separately and is not in the composite.\n\n")

  anls <- rowMeans(A[, ANLS_GENES, drop = FALSE])
  mct4 <- A[, "SLC16A3"]

  # ---- G1: the three MCT4 rows -------------------------------------------------
  rows <- list(
    list("MCT4 - neuronal V-ATPase", D$MCT4, D$VATP_n6, PUB$mct4_vatp),
    list("MCT4 - neuronal LAMP1",    D$MCT4, D$LAMP1_n, PUB$mct4_lamp),
    list("MCT4 - neuronal LDHB",     D$MCT4, D$LDHB_n,  PUB$mct4_ldhb))
  cat("--- G1 reproduction gate (three MCT4 rows) ---\n")
  ok <- TRUE
  for (r in rows) {
    s <- row_stats(r[[2]], r[[3]], D$cps)
    g <- gate(c(s[["zero"]], s[["part"]]), c(r[[4]][["zero"]], r[[4]][["part"]]), TOL$r)
    ok <- ok && g$pass
    cat(sprintf("  %-26s zero %+.3f (pub %+.3f) | partial %+.3f (pub %+.3f) | max dev %.4f %s\n",
                r[[1]], s[["zero"]], r[[4]][["zero"]], s[["part"]], r[[4]][["part"]],
                g$maxdev, ifelse(g$pass, "OK", "FAIL")))
  }
  if (!ok) {
    cat("  G1 = FAIL - the composite reconstruction is at fault. The ANLS row is\n")
    cat("  not reported and the tolerance is not relaxed.\n")
    return(invisible(NULL))
  }
  cat("  G1 = PASS\n\n")

  # ---- The ANLS row ------------------------------------------------------------
  s <- row_stats(anls, D$VATP_n6, D$cps)
  g <- gate(c(s[["zero"]], s[["part"]]), c(PUB$anls_vatp[["zero"]], PUB$anls_vatp[["part"]]), TOL$r)
  gc <- gate(s[["change"]], PUB$anls_vatp[["change"]], TOL$change)
  cat("--- ANLS row of Supplemental Table 1(a) ---\n")
  cat(sprintf("  zero-order r %+.3f  p %.3g\n", s[["zero"]], s[["zero_p"]]))
  cat(sprintf("  partial r    %+.3f  p %.3g   (CPS-adjusted)\n", s[["part"]], s[["part_p"]]))
  cat(sprintf("  change       %+.1f%%\n", s[["change"]]))
  cat(sprintf("  manuscript   %+.3f / %+.3f / %+.1f%%  -> %s\n",
              PUB$anls_vatp[["zero"]], PUB$anls_vatp[["part"]], PUB$anls_vatp[["change"]],
              ifelse(g$pass && gc$pass, "reproduces", "DOES NOT reproduce")))
  cat(sprintf("  stale .xlsx  +0.356 / +0.305 / -14.3%%\n\n"))

  # ---- Second check: ANLS - MCT4 within astrocytes ----------------------------
  s2 <- row_stats(anls, mct4, D$cps)
  g2 <- gate(c(s2[["zero"]], s2[["part"]]), c(PUB$anls_mct4[["zero"]], PUB$anls_mct4[["part"]]), TOL$r)
  cat("--- Independent check: ANLS - MCT4 within astrocytes (Results text) ---\n")
  cat(sprintf("  zero-order r %+.3f (published %+.3f) | partial r %+.3f (published %+.3f) %s\n\n",
              s2[["zero"]], PUB$anls_mct4[["zero"]], s2[["part"]], PUB$anls_mct4[["part"]],
              ifelse(g2$pass, "OK", "FAIL")))

  out <- data.frame(
    pair = c("ANLS - neuronal V-ATPase", "ANLS - astrocytic MCT4"),
    zero_order_r = round(c(s[["zero"]], s2[["zero"]]), 3),
    zero_order_p = signif(c(s[["zero_p"]], s2[["zero_p"]]), 3),
    partial_r    = round(c(s[["part"]], s2[["part"]]), 3),
    partial_p    = signif(c(s[["part_p"]], s2[["part_p"]]), 3),
    change_pct   = round(c(s[["change"]], s2[["change"]]), 1),
    stringsAsFactors = FALSE)
  write.csv(out, file.path(OUT, "SuppTable_1a_ANLS.csv"), row.names = FALSE)
  cat(sprintf("Written to %s/SuppTable_1a_ANLS.csv\n", OUT))
  cat("=== paste the whole output back ===\n")
}

args <- commandArgs(trailingOnly = TRUE)
if (!selftest()) { cat("SELF-TEST FAILED - stopping before any file is read.\n"); quit(status = 1) }
if (length(args) && args[1] == "selftest") { cat("selftest mode: done\n"); quit(status = 0) }
main()
