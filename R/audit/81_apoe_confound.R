#!/usr/bin/env Rscript
# =============================================================================
# 81_apoe_confound.R  -  standalone; paths are relative to the repository root
#
# WHY
# Ayton et al. (Nat Commun 6:6760, 2015) show that APOE genotype sets CSF
# ferritin: epsilon-4 carriers run 22% higher (p = 1.1e-8), and ferritin tracks
# CSF ApoE with partial R2 = 0.341. If APOE governs brain iron handling, then
# any iron-axis result in a cohort whose APOE composition shifts with disease
# stage can be an APOE composition effect rather than a disease effect.
#
# P2 reports an iron pathway (TFRC, FTH1, FTL, HMOX1) alongside the MCT4 result.
# This script asks whether either survives adjustment for APOE, and whether the
# cohort's APOE composition itself moves with the progression axis.
#
# The reviewer question this answers, before it is asked: "is your iron axis an
# APOE composition effect?"
#
# INPUT   data/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad   (obs only)
#         output/sensitivity/donor_level.csv
#         output/sensitivity/donor_by_gene.rds
# X is not read.
#
# Usage:
#   Rscript 81_apoe_confound.R selftest
#   Rscript 81_apoe_confound.R
# =============================================================================

H5AD <- file.path("data", "SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad")
OUT  <- file.path("output", "sensitivity")
CSV  <- file.path(OUT, "donor_level.csv")
RDS  <- file.path(OUT, "donor_by_gene.rds")
SEED <- 20260720

IRON <- c("TFRC", "FTH1", "FTL", "HMOX1", "SLC40A1", "SLC11A2")

## ---- Published values, used as a reproduction gate --------------------------
PUB <- list(slope = -0.0737, partial_r = c(0.474, 0.501, 0.418))
TOL <- list(slope = 0.0003, partial_r = 0.005)

## ---- Verdict, fixed before the data are read --------------------------------
# G1  The donor table must reproduce the published slope and partial
#     correlations. If not, nothing below is interpreted.
#
# A1  APOE composition along the axis.
#     If epsilon-4 carrier fraction is associated with CPS at p < 0.05, the
#     composition path is open and every downstream result must be reported
#     with and without APOE adjustment.
#
# A2  Material confounding, judged on the estimate rather than on significance.
#     A marker is called APOE-confounded if adjustment for epsilon-4 status
#     moves its CPS slope by more than 20% of the unadjusted value, or moves a
#     partial correlation by more than 0.05.
#     These thresholds are fixed here and are not adjusted after the fact.
CONFOUND_SLOPE <- 0.20
CONFOUND_R     <- 0.05

## ---- Functions ---------------------------------------------------------------
pc <- function(x, y, Z) {
  Z <- as.matrix(Z); o <- complete.cases(x, y, Z)
  x <- x[o]; y <- y[o]; Z <- Z[o, , drop = FALSE]
  r <- cor(residuals(lm(x ~ Z)), residuals(lm(y ~ Z)))
  n <- sum(o); k <- ncol(Z)
  c(r = r, p = 2 * pt(-abs(r * sqrt((n - k - 2) / (1 - r^2))), n - k - 2))
}

gate <- function(obs, ref, tol) {
  d <- abs(obs - ref)
  list(pass = all(is.finite(d)) && all(d <= tol * (1 + 1e-9)),
       maxdev = if (all(is.finite(d))) max(d) else NA_real_)
}

obs_cat <- function(h5, col, n_expect = NULL) {
  p1 <- paste0("obs/", col)
  dec <- function(codes, cats) {
    out <- rep(NA_character_, length(codes))
    ok <- !is.na(codes) & codes >= 0
    out[ok] <- as.character(cats)[codes[ok] + 1L]; out
  }
  v <- if (h5$has(paste0(p1, "/codes")))
         dec(as.integer(h5$all(paste0(p1, "/codes"))), h5$all(paste0(p1, "/categories")))
       else if (h5$has(paste0("obs/__categories/", col)))
         dec(as.integer(h5$all(p1)), h5$all(paste0("obs/__categories/", col)))
       else as.character(h5$all(p1))
  if (!is.null(n_expect) && length(v) != n_expect)
    stop(sprintf("obs/%s has length %d but the file has %d cells.", col, length(v), n_expect))
  v
}

make_h5 <- function(path) {
  if (!requireNamespace("hdf5r", quietly = TRUE)) stop("hdf5r required")
  f <- hdf5r::H5File$new(path, "r")
  nm <- function(g) tryCatch(if (g == "") names(f) else names(f[[g]]), error = function(e) character(0))
  list(ls = nm,
       has = function(n) { p <- strsplit(n, "/")[[1]]
         for (i in seq_along(p)) { par <- if (i == 1) "" else paste(p[1:(i-1)], collapse = "/")
           if (!(p[i] %in% nm(par))) return(FALSE) }; TRUE },
       all = function(n) f[[n]][], close = function() f$close_all())
}

# Slope of y on cps, optionally adjusting for extra covariates.
slope_of <- function(y, d, extra = NULL) {
  rhs <- if (is.null(extra)) "cps" else paste("cps +", paste(extra, collapse = " + "))
  s <- summary(lm(as.formula(paste("y ~", rhs)), data = cbind(d, y = y)))$coefficients
  c(beta = s["cps", 1], p = s["cps", 4])
}

## ---- Self-test ---------------------------------------------------------------
selftest <- function() {
  cat("=== SELF-TEST ===\n"); set.seed(SEED); f <- character(0)

  mk <- list(has = function(n) n %in% c("obs/G", "obs/__categories/G"),
             all = function(n) if (n == "obs/G") c(0L, -1L, 1L) else c("3/3", "3/4"))
  t1 <- identical(obs_cat(mk, "G", 3L), c("3/3", NA, "3/4"))
  cat(sprintf("  T1 obs_cat preserves NA codes                : %s\n", ifelse(t1, "PASS", "FAIL")))
  if (!t1) f <- c(f, "T1")

  # T2: the confound rule must fire on a planted confound and stay silent on a
  #     clean one. Without this the whole script could return "no confound"
  #     because the detector does not work.
  n <- 84
  cps <- runif(n)
  e4  <- rbinom(n, 1, plogis(-2 + 4 * cps))          # carrier fraction rises with CPS
  y_conf  <- 1 - 0.5 * e4 + rnorm(n, 0, 0.05)        # driven by APOE only
  y_clean <- 1 - 0.5 * cps + rnorm(n, 0, 0.05)       # driven by CPS only
  d <- data.frame(cps = cps, e4 = e4)
  mv <- function(y) {
    a <- slope_of(y, d)[["beta"]]; b <- slope_of(y, d, "e4")[["beta"]]
    abs(a - b) / abs(a)
  }
  t2 <- mv(y_conf) > CONFOUND_SLOPE && mv(y_clean) < CONFOUND_SLOPE
  cat(sprintf("  T2 confound rule fires on planted, not clean : %s (%.2f vs %.2f)\n",
              ifelse(t2, "PASS", "FAIL"), mv(y_conf), mv(y_clean)))
  if (!t2) f <- c(f, "T2")

  # T3: pc() against theory over 2000 replicates rather than a single draw.
  R <- 2000; z <- numeric(R); pv <- numeric(R); nn <- 84; k <- 2
  for (i in seq_len(R)) {
    o <- pc(rnorm(nn), rnorm(nn), cbind(rnorm(nn), rnorm(nn)))
    z[i] <- atanh(o[["r"]]); pv[i] <- o[["p"]]
  }
  sdo <- sd(z); sdt <- 1 / sqrt(nn - k - 3); a05 <- mean(pv < 0.05)
  t3 <- abs(sdo / sdt - 1) <= 0.08 && a05 >= 0.036 && a05 <= 0.065
  cat(sprintf("  T3 pc() Fisher-z %.4f vs theory %.4f, a=%.4f : %s\n",
              sdo, sdt, a05, ifelse(t3, "PASS", "FAIL")))
  if (!t3) f <- c(f, "T3")

  t4 <- gate(PUB$partial_r, PUB$partial_r, TOL$partial_r)$pass &&
       !gate(PUB$partial_r + c(0, 0, 0.02), PUB$partial_r, TOL$partial_r)$pass
  cat(sprintf("  T4 gate accepts truth, rejects error         : %s\n", ifelse(t4, "PASS", "FAIL")))
  if (!t4) f <- c(f, "T4")

  ok <- length(f) == 0
  cat(sprintf("=== SELF-TEST %s%s ===\n\n", ifelse(ok, "PASS", "FAIL"),
              ifelse(ok, "", paste0(" (", paste(f, collapse = ", "), ")"))))
  ok
}

## ---- Main ---------------------------------------------------------------------
main <- function() {
  for (p in c(CSV, RDS, H5AD))
    if (!file.exists(p)) { cat(sprintf("Not found: %s\n", p)); quit(status = 1) }
  D <- read.csv(CSV, stringsAsFactors = FALSE)
  Z <- readRDS(RDS); A <- Z$astro; colnames(A) <- Z$genes
  stopifnot(identical(as.character(Z$donors), as.character(D$donor)))

  h5 <- make_h5(H5AD); on.exit(try(h5$close(), silent = TRUE))
  cps <- as.numeric(h5$all("obs/Continuous Pseudo-progression Score")); nC <- length(cps)
  don <- obs_cat(h5, "Donor ID", nC)

  # The obs column name is reported rather than assumed.
  cand <- c("APOE Genotype", "APOE4 Status", "APOE genotype", "apoe_genotype", "APOE")
  hit <- cand[vapply(cand, function(x) h5$has(paste0("obs/", x)), logical(1))]
  if (!length(hit)) {
    cat("APOE column not found. obs contains:\n")
    print(head(h5$ls("obs"), 60)); quit(status = 1)
  }
  cat(sprintf("APOE column: obs/%s\n", hit[1]))
  gt <- obs_cat(h5, hit[1], nC)
  cat("Genotype distribution (cells):\n"); print(table(gt, useNA = "ifany"))

  # Donor-level genotype, with a check that it is constant within donor.
  inCoh <- !is.na(cps)
  mixed <- names(which(tapply(gt[inCoh], don[inCoh],
                              function(z) length(unique(z[!is.na(z)]))) > 1))
  if (length(mixed)) cat(sprintf("\n[warn] genotype varies within donor: %s\n", paste(mixed, collapse = ", ")))
  gd <- tapply(gt[inCoh], don[inCoh], function(z) { u <- unique(z[!is.na(z)]); if (length(u) == 1) u else NA })
  D$apoe <- unname(gd[as.character(D$donor)])
  D$e4 <- ifelse(is.na(D$apoe), NA, as.integer(grepl("4", D$apoe)))
  cat(sprintf("\nDonors: %d | genotype known %d | epsilon-4 carriers %d (%.1f%%)\n",
              nrow(D), sum(!is.na(D$e4)), sum(D$e4 == 1, na.rm = TRUE),
              100 * mean(D$e4 == 1, na.rm = TRUE)))

  # ---- G1 -------------------------------------------------------------------
  s <- summary(lm(MCT4 ~ cps, data = D))$coefficients["cps", c(1, 4)]
  A3 <- vapply(c("VATP_n6", "LAMP1_n", "LDHB_n"),
               function(v) pc(D$MCT4, D[[v]], D[, "cps", drop = FALSE])[["r"]], numeric(1))
  gS <- gate(s[1], PUB$slope, TOL$slope); gA <- gate(A3, PUB$partial_r, TOL$partial_r)
  cat(sprintf("\n--- G1 ---\n  slope %+.4f (dev %.4f) %s | partial %+.3f/%+.3f/%+.3f (dev %.3f) %s\n",
              s[1], gS$maxdev, ifelse(gS$pass, "OK", "FAIL"),
              A3[1], A3[2], A3[3], gA$maxdev, ifelse(gA$pass, "OK", "FAIL")))
  if (!(gS$pass && gA$pass)) { cat("  G1 = FAIL - nothing below is interpreted.\n"); return(invisible(NULL)) }
  cat("  G1 = PASS\n")

  d <- D[!is.na(D$e4), ]

  # ---- A1: is the composition path open? ------------------------------------
  a1 <- summary(glm(e4 ~ cps, data = d, family = binomial))$coefficients["cps", c(1, 4)]
  cat(sprintf("\n--- A1 APOE composition along the axis ---\n"))
  cat(sprintf("  logistic e4 ~ cps : beta %+.3f, p %.3g -> %s\n", a1[1], a1[2],
              ifelse(a1[2] < 0.05, "composition path OPEN",
                     "no detectable composition shift")))
  for (q in list(c(0, .33), c(.33, .66), c(.66, 1))) {
    k <- d$cps >= quantile(d$cps, q[1]) & d$cps <= quantile(d$cps, q[2])
    cat(sprintf("    CPS %.2f-%.2f : %d donors, e4 %.0f%%\n",
                quantile(d$cps, q[1]), quantile(d$cps, q[2]), sum(k), 100 * mean(d$e4[k])))
  }

  # ---- A2: does adjustment move the estimates? ------------------------------
  cat("\n--- A2 CPS slope, unadjusted vs epsilon-4 adjusted ---\n")
  cat(sprintf("  %-10s %11s %11s %9s  %s\n", "marker", "unadj beta", "adj beta", "shift", "verdict"))
  res <- list()
  marks <- c(list(MCT4 = d$MCT4), setNames(lapply(IRON, function(g)
    if (g %in% colnames(A)) A[match(d$donor, D$donor), g] else NULL), IRON))
  for (nm in names(marks)) {
    y <- marks[[nm]]
    if (is.null(y)) { cat(sprintf("  %-10s  not in matrix\n", nm)); next }
    u <- slope_of(y, d); a <- slope_of(y, d, "e4")
    sh <- abs(u[["beta"]] - a[["beta"]]) / abs(u[["beta"]])
    v <- if (sh > CONFOUND_SLOPE) "APOE-CONFOUNDED" else "robust"
    cat(sprintf("  %-10s %+11.5f %+11.5f %8.1f%%  %s\n", nm, u[["beta"]], a[["beta"]], 100 * sh, v))
    res[[length(res) + 1]] <- data.frame(marker = nm, beta_unadj = u[["beta"]], p_unadj = u[["p"]],
      beta_adj = a[["beta"]], p_adj = a[["p"]], shift_pct = 100 * sh, verdict = v,
      stringsAsFactors = FALSE)
  }

  cat("\n--- A2 cross-cellular partial correlations ---\n")
  cat(sprintf("  %-22s %9s %9s %8s  %s\n", "pair", "CPS-adj", "+e4", "delta", "verdict"))
  for (v in c("VATP_n6", "LAMP1_n", "LDHB_n")) {
    r1 <- pc(d$MCT4, d[[v]], d[, "cps", drop = FALSE])[["r"]]
    r2 <- pc(d$MCT4, d[[v]], d[, c("cps", "e4")])[["r"]]
    vd <- if (abs(r1 - r2) > CONFOUND_R) "APOE-CONFOUNDED" else "robust"
    cat(sprintf("  MCT4-%-17s %+9.3f %+9.3f %8.3f  %s\n", v, r1, r2, r2 - r1, vd))
    res[[length(res) + 1]] <- data.frame(marker = paste0("partial MCT4-", v),
      beta_unadj = r1, p_unadj = NA, beta_adj = r2, p_adj = NA,
      shift_pct = 100 * abs(r1 - r2), verdict = vd, stringsAsFactors = FALSE)
  }

  out <- do.call(rbind, res)
  write.csv(out, file.path(OUT, "apoe_confound.csv"), row.names = FALSE)
  n_conf <- sum(out$verdict == "APOE-CONFOUNDED")
  cat(sprintf("\n  markers called APOE-confounded: %d / %d\n", n_conf, nrow(out)))
  cat("  Thresholds were fixed before the data were read and are not adjusted.\n")
  cat(sprintf("\nWritten to %s/apoe_confound.csv\n", OUT))
  cat("=== paste the whole output back ===\n")
}

args <- commandArgs(trailingOnly = TRUE)
if (!selftest()) { cat("SELF-TEST FAILED - stopping before any file is read.\n"); quit(status = 1) }
if (length(args) && args[1] == "selftest") { cat("selftest mode: done\n"); quit(status = 0) }
main()
