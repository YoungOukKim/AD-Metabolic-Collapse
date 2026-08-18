#!/usr/bin/env Rscript
# ---------------------------------------------------------------------
# 03_specificity_screen.R
#
# Screening tool. Given a candidate gene, ask whether its CPS slope and its
# coupling to MCT4 stand out against a null of genes matched on mean expression.
#
# Scope limit, stated up front: the cached matrix holds donor-level means, so
# per-cell detection rate cannot be recovered. This null is matched on
# expression level, not on detection rate. It is a screen. A gene that passes
# here still has to clear a cell-level detection-matched null before it can
# enter a manuscript.
#
# Usage:
#   Rscript R/03_specificity_screen.R selftest
#   Rscript R/03_specificity_screen.R [GENE ...]      # default: BSG
# ---------------------------------------------------------------------

source(file.path("R", "sensitivity", "config_sensitivity.R")); source(file.path("R", "sensitivity", "00_sensitivity_utils.R"))
CSV <- file.path(OUT_DIR, "donor_level.csv")
RDS <- file.path(OUT_DIR, "donor_by_gene.rds")

## ---- Verdict, fixed before the data are read ------------------------
# B1  CPS slope at or below the 5th percentile of the expression-matched null.
# B2  MCT4 coupling (partial r given CPS) above the 95th percentile of the same null.
# Both must pass. SLC16A3 is screened first as a positive control: if it does
# not land in the extreme tail, the screen is not calibrated and the candidate
# result is not interpreted.
B_LOW <- 0.05; B_HIGH <- 0.95; MATCH_TOL <- 0.20; MIN_NULL <- 200

selftest <- function() {
  cat("=== SELF-TEST ===\n"); set.seed(SEED); f <- character(0)

  # slopes_all must agree with lm() gene by gene.
  n <- 40; g <- 25
  cps <- runif(n); gm <- rnorm(n); G <- matrix(rnorm(n * g), n, g)
  s1 <- slopes_all(G, cps, gm)
  s2 <- vapply(seq_len(g), function(j) coef(lm(G[, j] ~ cps + gm))[["cps"]], numeric(1))
  t1 <- max(abs(s1 - s2)) < 1e-9
  cat(sprintf("  T1 slopes_all == lm, gene by gene      : %s (max dev %.2e)\n",
              ifelse(t1, "PASS", "FAIL"), max(abs(s1 - s2))))
  if (!t1) f <- c(f, "T1")

  # matched_pctl must reproduce an analytic answer: with statistics 1..100 and
  # the target at rank 50, the percentile is 49/99.
  t2 <- abs(matched_pctl(as.numeric(1:100), 50L, rep(1, 100), min_n = 50)$pctl - 49 / 99) < 1e-9
  cat(sprintf("  T2 matched_pctl == analytic answer     : %s\n", ifelse(t2, "PASS", "FAIL")))
  if (!t2) f <- c(f, "T2")

  # The screen must separate a planted signal from the null. A gene built to
  # track the exposure should sit in the lower tail; a random gene should not.
  set.seed(SEED)
  nD <- 84; nG <- 1000
  cpsv <- sort(runif(nD)); gmv <- rnorm(nD)
  Gm <- matrix(rnorm(nD * nG, 0.5, 0.1), nD, nG)
  Gm[, 1] <- 0.5 - 0.3 * cpsv + rnorm(nD, 0, 0.01)   # planted decline
  sl <- slopes_all(Gm, cpsv, gmv); ex <- colMeans(Gm)
  t3 <- matched_pctl(sl, 1L, ex)$pctl <= B_LOW && matched_pctl(sl, 500L, ex)$pctl > B_LOW
  cat(sprintf("  T3 screen finds planted signal, not noise : %s\n", ifelse(t3, "PASS", "FAIL")))
  if (!t3) f <- c(f, "T3")

  n <- 84; k <- 2; R <- 2000; z <- numeric(R); pv <- numeric(R)
  for (i in seq_len(R)) {
    o <- pc(rnorm(n), rnorm(n), cbind(rnorm(n), rnorm(n)))
    z[i] <- atanh(o[["r"]]); pv[i] <- o[["p"]]
  }
  t4 <- abs(sd(z) / (1 / sqrt(n - k - 3)) - 1) <= 0.08 &&
        mean(pv < 0.05) >= 0.036 && mean(pv < 0.05) <= 0.065
  cat(sprintf("  T4 pc() against theory over 2000 reps  : %s\n", ifelse(t4, "PASS", "FAIL")))
  if (!t4) f <- c(f, "T4")

  ok <- length(f) == 0
  cat(sprintf("=== SELF-TEST %s%s ===\n\n", ifelse(ok, "PASS", "FAIL"),
              ifelse(ok, "", paste0(" (", paste(f, collapse = ", "), ")"))))
  ok
}

main <- function(targets) {
  for (p in c(CSV, RDS))
    if (!file.exists(p)) { cat(sprintf("Not found: %s - run 01_extract_h5ad.R with FULL_GENE_MATRIX = TRUE.\n", p)); quit(status = 1) }
  D <- read.csv(CSV, stringsAsFactors = FALSE)
  Z <- readRDS(RDS)
  stopifnot(identical(as.character(Z$donors), as.character(D$donor)))

  G <- Z$astro
  expr <- colMeans(G, na.rm = TRUE)
  keep <- which(is.finite(expr) & expr > 0)
  G <- G[, keep, drop = FALSE]; expr <- expr[keep]; genes <- Z$genes[keep]
  sl <- slopes_all(G, D$cps, D$gm_a)
  iM <- match("SLC16A3", genes)
  if (is.na(iM)) { cat("SLC16A3 absent - cannot calibrate the screen.\n"); quit(status = 1) }

  cat("--- Expression-matched null screen (not detection-matched) ---\n")
  report <- function(g, label) {
    i <- match(g, genes)
    if (is.na(i)) { cat(sprintf("  %-24s not present in the matrix\n", label)); return(NULL) }
    mp <- matched_pctl(sl, i, expr, MATCH_TOL, MIN_NULL)
    r  <- pc(G[, iM], G[, i], D[, "cps", drop = FALSE])[["r"]]
    e0 <- expr[i]
    idx <- which(expr >= e0 * (1 - MATCH_TOL) & expr <= e0 * (1 + MATCH_TOL))
    if (length(idx) < MIN_NULL) idx <- head(order(abs(expr - e0)), MIN_NULL)
    idx <- setdiff(idx, c(i, iM))
    rn <- vapply(idx, function(j) pc(G[, iM], G[, j], D[, "cps", drop = FALSE])[["r"]], numeric(1))
    pr <- mean(rn < r, na.rm = TRUE)
    cat(sprintf("  %-24s mean expr %.4f | CPS slope %+.5f -> %.2f pctl of %d | MCT4 partial r %+.3f -> %.2f pctl\n",
                label, e0, sl[i], mp$pctl, mp$n, r, pr))
    c(slope_pctl = mp$pctl, r_pctl = pr)
  }
  ctrl <- report("SLC16A3", "SLC16A3 (positive control)")
  if (!is.null(ctrl) && ctrl[["slope_pctl"]] > B_LOW) {
    cat("\n  The positive control did not reach the lower tail. The screen is not\n")
    cat("  calibrated on this matrix and candidate results are not interpreted.\n")
    return(invisible(NULL))
  }
  cat("\n")
  for (g in targets) {
    v <- report(g, g)
    if (is.null(v)) next
    b1 <- v[["slope_pctl"]] <= B_LOW; b2 <- v[["r_pctl"]] > B_HIGH
    cat(sprintf("    B1 slope in lower %.0f%%: %s | B2 coupling above %.0f%%: %s -> %s\n\n",
                B_LOW * 100, ifelse(b1, "PASS", "FAIL"), B_HIGH * 100, ifelse(b2, "PASS", "FAIL"),
                ifelse(b1 && b2,
                       "PASS - eligible for a cell-level detection-matched null.",
                       "FAIL - closed. The threshold is not moved.")))
  }
}

args <- commandArgs(trailingOnly = TRUE)
if (!selftest()) { cat("SELF-TEST FAILED - stopping before any file is read.\n"); quit(status = 1) }
if (length(args) && args[1] == "selftest") { cat("selftest mode: done\n"); quit(status = 0) }
main(if (length(args)) args else "BSG")
