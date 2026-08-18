#!/usr/bin/env Rscript
# ---------------------------------------------------------------------
# 02_sa_sensitivity.R
#
# Severely-affected donor sensitivity analysis, run from the cached output of
# 01_extract_h5ad.R plus the obs columns of the h5ad. X is not read again.
#
# The eleven donors designated 'severely affected' by the atlas study are
# flagged in the released metadata as obs/Severely Affected Donor. They are
# retained in the primary analysis; this script removes them and reports what
# changes.
#
# Usage:
#   Rscript R/02_sa_sensitivity.R selftest
#   Rscript R/02_sa_sensitivity.R
# ---------------------------------------------------------------------

source(file.path("R", "sensitivity", "config_sensitivity.R")); source(file.path("R", "sensitivity", "00_sensitivity_utils.R"))
CSV <- file.path(OUT_DIR, "donor_level.csv")

## ---- Published values, used as a reproduction gate ------------------
PUB <- list(slope = -0.0737, partial_r = c(0.474, 0.501, 0.418), bin_pct = -43.2)
TOL <- list(slope = 0.0003, partial_r = 0.005)

## ---- Verdict, fixed before the data are read ------------------------
# G1 reproduction gate. If the cached table does not reproduce the published
# donor-level values, the sensitivity result is not interpretable and is not
# printed. Thresholds are not relaxed to obtain a pass.

selftest <- function() {
  cat("=== SELF-TEST ===\n"); set.seed(SEED); f <- character(0)

  mk <- list(has = function(n) n %in% c("obs/F", "obs/__categories/F"),
             all = function(n) if (n == "obs/F") c(0L, -1L, 1L, -1L, 0L) else c("N", "Y"))
  t1 <- identical(obs_cat(mk, "F", 5L), c("N", NA, "Y", NA, "N")) &&
        inherits(try(obs_cat(mk, "F", 4L), silent = TRUE), "try-error")
  cat(sprintf("  T1 obs_cat NA codes + length assertion : %s\n", ifelse(t1, "PASS", "FAIL")))
  if (!t1) f <- c(f, "T1")

  # The gate must accept a true match and reject a false one. Test values are
  # kept away from the threshold so that the outcome is decided by logic.
  t2 <- gate(PUB$partial_r, PUB$partial_r, TOL$partial_r)$pass &&
       !gate(PUB$partial_r + c(0, 0, 0.02), PUB$partial_r, TOL$partial_r)$pass &&
        gate(-0.0738, PUB$slope, TOL$slope)$pass &&
       !gate(-0.0747, PUB$slope, TOL$slope)$pass
  cat(sprintf("  T2 gate accepts truth, rejects error   : %s\n", ifelse(t2, "PASS", "FAIL")))
  if (!t2) f <- c(f, "T2")

  # pc() is checked against theory rather than a single draw: over 2000
  # replicates the Fisher-z standard deviation must approach 1/sqrt(n-k-3).
  n <- 84; k <- 2; R <- 2000; z <- numeric(R); pv <- numeric(R)
  for (i in seq_len(R)) {
    o <- pc(rnorm(n), rnorm(n), cbind(rnorm(n), rnorm(n)))
    z[i] <- atanh(o[["r"]]); pv[i] <- o[["p"]]
  }
  sdo <- sd(z); sdt <- 1 / sqrt(n - k - 3); a05 <- mean(pv < 0.05)
  t3 <- abs(sdo / sdt - 1) <= 0.08 && a05 >= 0.036 && a05 <= 0.065
  cat(sprintf("  T3 Fisher-z sd %.4f vs theory %.4f, alpha %.4f : %s\n",
              sdo, sdt, a05, ifelse(t3, "PASS", "FAIL")))
  if (!t3) f <- c(f, "T3")

  # mediate() must recover known proportions on constructed data, in both the
  # full and the partial case. The mediator is deliberately not collinear with
  # the exposure: a near-deterministic mediator would make the adjusted
  # coefficient unstable and the test uninformative.
  # A single draw would not settle this: with n = 1000 the adjusted coefficient
  # carries a standard error of a few percent, so one replicate can land several
  # points away from the truth by chance. The judgement is made on the median of
  # 200 replicates instead.
  set.seed(SEED)
  R4 <- 200; nn <- 1000
  pf <- numeric(R4); pp <- numeric(R4)
  for (i in seq_len(R4)) {
    cps <- rnorm(nn)
    M   <- cps + rnorm(nn)                 # mediator, half signal half noise
    dd  <- data.frame(cps = cps, MCT4 = M,
                      Yf = M + rnorm(nn),          # full: no direct path
                      Yp = M + cps + rnorm(nn))    # partial: equal direct path
    pf[i] <- mediate(dd, "Yf")[["prop"]]
    pp[i] <- mediate(dd, "Yp")[["prop"]]
  }
  mf <- median(pf); mp <- median(pp)
  t4 <- abs(mf - 100) < 2 && abs(mp - 50) < 2
  cat(sprintf("  T4 mediate() median over %d reps: %.1f%% (100) and %.1f%% (50) : %s\n",
              R4, mf, mp, ifelse(t4, "PASS", "FAIL")))
  if (!t4) f <- c(f, "T4")

  ok <- length(f) == 0
  cat(sprintf("=== SELF-TEST %s%s ===\n\n", ifelse(ok, "PASS", "FAIL"),
              ifelse(ok, "", paste0(" (", paste(f, collapse = ", "), ")"))))
  ok
}

summarise <- function(s) {
  q <- summary(lm(MCT4 ~ cps, data = s))$coefficients["cps", c(1, 4)]
  A <- vapply(c("VATP_n6", "LAMP1_n", "LDHB_n"),
              function(v) pc(s$MCT4, s[[v]], s[, "cps", drop = FALSE])[["r"]], numeric(1))
  C <- vapply(c("VATP_n6", "LAMP1_n", "LDHB_n"),
              function(v) pc(s$MCT4, s[[v]], s[, c("cps", "gm_a", "gm_n")])[["r"]], numeric(1))
  b  <- tapply(s$MCT4 * s$ncell_a, round(s$cps, 1), sum) /
        tapply(s$ncell_a, round(s$cps, 1), sum)
  bn <- as.numeric(names(b))
  list(n = nrow(s), slope = q[1], p = q[2], A = A, C = C,
       pct = 100 * (mean(b[bn %in% c(0.6, 0.7, 0.8)]) / mean(b[bn %in% c(0.2, 0.3, 0.4)]) - 1))
}

main <- function() {
  if (!file.exists(CSV)) { cat(sprintf("Not found: %s - run 01_extract_h5ad.R first.\n", CSV)); quit(status = 1) }
  if (!file.exists(H5AD)) { cat(sprintf("Not found: %s\n", H5AD)); quit(status = 1) }
  D <- read.csv(CSV, stringsAsFactors = FALSE)

  h5 <- make_h5(H5AD); on.exit(try(h5$close(), silent = TRUE))
  cps <- as.numeric(h5$all("obs/Continuous Pseudo-progression Score")); nC <- length(cps)
  don <- obs_cat(h5, "Donor ID", nC)
  flag <- obs_cat(h5, "Severely Affected Donor", nC)
  cat("Severely-affected flag distribution (NA = neurotypical reference):\n")
  print(table(flag, useNA = "ifany"))

  inCoh <- !is.na(cps); fv <- flag[inCoh]; dv <- don[inCoh]
  mixed <- names(which(tapply(fv, dv, function(z) length(unique(z[!is.na(z)]))) > 1))
  cat(sprintf("Donors with mixed flag values: %s\n",
              if (length(mixed)) paste(mixed, collapse = ", ") else "none (donor-level attribute)"))
  SA <- sort(unique(dv[fv %in% c("Y", "Yes", "True", "TRUE", "true", "1")]))
  cat(sprintf("\nSeverely-affected donors (n = %d):\n  %s\n", length(SA), paste(SA, collapse = ", ")))
  cat(sprintf("CPS range of those donors: %.3f - %.3f (cohort %.3f - %.3f)\n",
              min(D$cps[D$donor %in% SA]), max(D$cps[D$donor %in% SA]), min(D$cps), max(D$cps)))

  a <- summarise(D)
  gS <- gate(a$slope, PUB$slope, TOL$slope); gA <- gate(a$A, PUB$partial_r, TOL$partial_r)
  cat("\n--- G1 reproduction gate ---\n")
  cat(sprintf("  slope %+.4f (published %+.4f, max dev %.4f) %s\n",
              a$slope, PUB$slope, gS$maxdev, ifelse(gS$pass, "OK", "FAIL")))
  cat(sprintf("  partial r %+.3f/%+.3f/%+.3f (published %+.3f/%+.3f/%+.3f, max dev %.3f) %s\n",
              a$A[1], a$A[2], a$A[3], PUB$partial_r[1], PUB$partial_r[2], PUB$partial_r[3],
              gA$maxdev, ifelse(gA$pass, "OK", "FAIL")))
  cat(sprintf("  bin-level MCT4 %.1f%% (published %.1f%%)\n", a$pct, PUB$bin_pct))
  if (!(gS$pass && gA$pass)) {
    cat("  G1 = FAIL - sensitivity results are not printed. The threshold is not moved.\n")
    return(invisible(NULL))
  }
  cat("  G1 = PASS\n")

  cat("\n--- Sensitivity (columns: neuronal V-ATPase / LAMP1 / LDHB) ---\n")
  for (lab in c("all donors", "SA excluded")) {
    s <- if (lab == "all donors") D else D[!(D$donor %in% SA), ]
    r <- summarise(s)
    cat(sprintf("  %-12s n=%-3d slope %+.4f (p=%.3g) | CPS-adj %+.3f/%+.3f/%+.3f | genome-wide %+.3f/%+.3f/%+.3f | bin %.1f%%\n",
                lab, r$n, r$slope, r$p, r$A[1], r$A[2], r$A[3], r$C[1], r$C[2], r$C[3], r$pct))
  }

  cat("\n--- Mediation, donor level (proportion mediated) ---\n")
  cat(sprintf("  %-26s %10s %10s\n", "", "all (84)", "SA excl (73)"))
  E <- D[!(D$donor %in% SA), ]
  for (Y in c("VATP_n6", "LAMP1_n"))
    cat(sprintf("  forward  MCT4 -> %-9s %9.1f%% %9.1f%%\n", Y,
                mediate(D, Y)[["prop"]], mediate(E, Y)[["prop"]]))
  rev_med <- function(d, M) {
    tot <- coef(lm(MCT4 ~ cps, data = d))[["cps"]]
    adj <- coef(lm(as.formula(paste("MCT4 ~ cps +", M)), data = d))[["cps"]]
    100 * (tot - adj) / tot
  }
  for (M in c("VATP_n6", "LAMP1_n"))
    cat(sprintf("  reverse  %-9s -> MCT4 %9.1f%% %9.1f%%\n", M, rev_med(D, M), rev_med(E, M)))

  write.csv(data.frame(donor = D$donor, severely_affected = D$donor %in% SA),
            file.path(OUT_DIR, "severely_affected_donors.csv"), row.names = FALSE)
  cat(sprintf("\nWritten to %s/severely_affected_donors.csv\n", OUT_DIR))
}

args <- commandArgs(trailingOnly = TRUE)
if (!selftest()) { cat("SELF-TEST FAILED - stopping before any file is read.\n"); quit(status = 1) }
if (length(args) && args[1] == "selftest") { cat("selftest mode: done\n"); quit(status = 0) }
main()
