#!/usr/bin/env Rscript
# =============================================================================
# 80_partial_df_audit.R  -  standalone; paths are relative to the repository root
#
# WHY
# Every partial correlation p-value in the manuscript was computed on n - 2
# degrees of freedom. For a partial correlation controlling k variables the
# correct value is n - k - 2. With one control variable that is an off-by-one;
# the published p-values are therefore very slightly anti-conservative. No
# conclusion changes, but a reviewer recomputing them will find the difference,
# and the specifications differ from row to row (different covariate sets,
# different sample sizes after CPS thresholding), so the correction cannot be
# applied by hand to the text.
#
# This script recomputes every specification from the donor-level data, prints
# the published value beside the corrected one, and flags any that change
# materially. A reproduction gate on the correlation coefficients runs first: if
# the r values do not reproduce, the specification has been mis-transcribed and
# no p-value is reported for it.
#
# INPUT   output/sensitivity/donor_level.csv        (from R/sensitivity/01_extract_h5ad.R)
#         output/sensitivity/donor_by_gene.rds      (ANLS composite only)
#
# Usage:
#   Rscript 80_partial_df_audit.R selftest
#   Rscript 80_partial_df_audit.R
# =============================================================================

OUT  <- file.path("output", "sensitivity")
CSV  <- file.path(OUT, "donor_level.csv")
RDS  <- file.path(OUT, "donor_by_gene.rds")
SEED <- 20260720
ANLS_GENES <- c("SLC2A1", "LDHA", "SLC16A1")

TOL_R <- 0.006          # reproduction tolerance on r
MATERIAL <- 0.05        # a change that crosses this threshold is flagged

## ---- Correlation with an explicit degrees-of-freedom argument ---------------
# df_mode "published" reproduces the manuscript (n - 2 regardless of k).
# df_mode "correct"   uses n - k - 2, which is the definition for a partial
# correlation controlling k variables. For a zero-order correlation k = 0 and
# the two agree, which is why rows 1-3 of the table need no change.
corr <- function(x, y, Z = NULL, df_mode = "correct") {
  if (is.null(Z)) {
    o <- complete.cases(x, y); r <- cor(x[o], y[o]); n <- sum(o); k <- 0L
  } else {
    Z <- as.matrix(Z); o <- complete.cases(x, y, Z)
    r <- cor(residuals(lm(x[o] ~ Z[o, , drop = FALSE])),
             residuals(lm(y[o] ~ Z[o, , drop = FALSE])))
    n <- sum(o); k <- ncol(Z)
  }
  df <- if (df_mode == "published") n - 2 else n - k - 2
  t  <- r * sqrt(df / (1 - r^2))
  c(r = r, p = 2 * pt(-abs(t), df), n = n, k = k, df = df)
}

# Some specifications residualise each variable on its own cell type's global
# mean rather than entering both means as shared covariates. The two are
# different procedures and give different coefficients. Degrees of freedom are
# not standard for this case; we use n minus the number of distinct covariates
# minus 2, and flag the row so the choice is visible rather than buried.
corr_own <- function(x, y, Zx, Zy, df_mode = "correct") {
  o <- complete.cases(x, y, Zx, Zy)
  rx <- residuals(lm(x[o] ~ as.matrix(Zx)[o, , drop = FALSE]))
  ry <- residuals(lm(y[o] ~ as.matrix(Zy)[o, , drop = FALSE]))
  r <- cor(rx, ry); n <- sum(o)
  k <- length(unique(c(colnames(Zx), colnames(Zy))))
  df <- if (df_mode == "published") n - 2 else n - k - 2
  t <- r * sqrt(df / (1 - r^2))
  c(r = r, p = 2 * pt(-abs(t), df), n = n, k = k, df = df)
}

## ---- Self-test ---------------------------------------------------------------
selftest <- function() {
  cat("=== SELF-TEST ===\n"); set.seed(SEED); f <- character(0)

  # T1: with no covariates the two modes must agree exactly.
  x <- rnorm(84); y <- 0.4 * x + rnorm(84)
  a <- corr(x, y); b <- corr(x, y, df_mode = "published")
  t1 <- abs(a[["p"]] - b[["p"]]) < 1e-15 && a[["df"]] == 82
  cat(sprintf("  T1 zero-order: both modes identical, df=82 : %s\n", ifelse(t1,"PASS","FAIL")))
  if (!t1) f <- c(f,"T1")

  # T2: with k covariates the correct df must be n-k-2 and the published df
  #     must stay at n-2, giving a strictly smaller p.
  Z <- data.frame(a = rnorm(84), b = rnorm(84), c = rnorm(84))
  a <- corr(x, y, Z); b <- corr(x, y, Z, "published")
  t2 <- a[["df"]] == 84 - 3 - 2 && b[["df"]] == 82 && b[["p"]] < a[["p"]]
  cat(sprintf("  T2 k=3: df %d vs %d, published p smaller     : %s\n",
              a[["df"]], b[["df"]], ifelse(t2,"PASS","FAIL")))
  if (!t2) f <- c(f,"T2")

  # T3: the null distribution of the corrected test must be calibrated. Over
  #     2000 replicates the type I error at 0.05 must land near 0.05, and the
  #     published-df test must be measurably anti-conservative. A single draw
  #     would not distinguish the two.
  R <- 2000; pc_ <- numeric(R); pp_ <- numeric(R)
  for (i in seq_len(R)) {
    xx <- rnorm(84); yy <- rnorm(84); ZZ <- data.frame(a = rnorm(84), b = rnorm(84), c = rnorm(84))
    pc_[i] <- corr(xx, yy, ZZ)[["p"]]; pp_[i] <- corr(xx, yy, ZZ, "published")[["p"]]
  }
  ac <- mean(pc_ < 0.05); ap <- mean(pp_ < 0.05)
  t3 <- ac >= 0.036 && ac <= 0.065 && ap >= ac
  cat(sprintf("  T3 alpha correct %.4f vs published %.4f     : %s\n", ac, ap, ifelse(t3,"PASS","FAIL")))
  if (!t3) f <- c(f,"T3")

  # T4: reproduce a published pair analytically. r = 0.474, n = 84, k = 1 must
  #     give p = 5.2e-6 on df = 82 and 6.0e-6 on df = 81.
  r <- 0.474
  p82 <- 2 * pt(-abs(r * sqrt(82/(1-r^2))), 82)
  p81 <- 2 * pt(-abs(r * sqrt(81/(1-r^2))), 81)
  t4 <- abs(p82 - 5.25e-6) < 2e-7 && abs(p81 - 6.00e-6) < 2e-7
  cat(sprintf("  T4 analytic: df82 %.2e, df81 %.2e            : %s\n", p82, p81, ifelse(t4,"PASS","FAIL")))
  if (!t4) f <- c(f,"T4")

  # T5: corr_own must equal the residual correlation computed by hand. This
  #     tests the implementation against its definition. An earlier version of
  #     this test asserted that the own-covariate and joint procedures give
  #     different answers; that is a claim about the data, not about the code,
  #     and it is false for constructions in which each covariate is
  #     informative only for its own variable.
  set.seed(SEED)
  n5 <- 84; ga <- rnorm(n5); gn <- rnorm(n5); cp <- rnorm(n5)
  xx <- cp + 2*ga + rnorm(n5); yy <- cp + 2*gn + rnorm(n5)
  Zx <- data.frame(cps = cp, gm_a = ga); Zy <- data.frame(cps = cp, gm_n = gn)
  w <- corr_own(xx, yy, Zx, Zy)
  manual <- cor(residuals(lm(xx ~ cp + ga)), residuals(lm(yy ~ cp + gn)))
  t5 <- abs(w[["r"]] - manual) < 1e-12 && w[["df"]] == n5 - 3 - 2
  cat(sprintf("  T5 corr_own == hand-computed residual r, df %d : %s\n",
              w[["df"]], ifelse(t5,"PASS","FAIL")))
  if (!t5) f <- c(f,"T5")

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
  Z <- readRDS(RDS); A <- Z$astro; colnames(A) <- Z$genes
  D$ANLS <- rowMeans(A[, ANLS_GENES, drop = FALSE])

  # Every specification that appears in the manuscript, stated explicitly.
  # thr = CPS threshold applied before the correlation.
  S <- list(
    list(id="zero MCT4-VATP",      x="MCT4", y="VATP_n6", Z=NULL,                     thr=0,    pub_r=0.517, pub_p=4.86e-7),
    list(id="zero MCT4-LAMP1",     x="MCT4", y="LAMP1_n", Z=NULL,                     thr=0,    pub_r=0.549, pub_p=6.30e-8),
    list(id="zero MCT4-LDHB",      x="MCT4", y="LDHB_n",  Z=NULL,                     thr=0,    pub_r=0.427, pub_p=5.02e-5),
    list(id="CPS  MCT4-VATP",      x="MCT4", y="VATP_n6", Z="cps",                    thr=0,    pub_r=0.474, pub_p=5.13e-6),
    list(id="CPS  MCT4-LAMP1",     x="MCT4", y="LAMP1_n", Z="cps",                    thr=0,    pub_r=0.501, pub_p=1.21e-6),
    list(id="CPS  MCT4-LDHB",      x="MCT4", y="LDHB_n",  Z="cps",                    thr=0,    pub_r=0.418, pub_p=7.62e-5),
    list(id="CPS  ANLS-VATP",      x="ANLS", y="VATP_n6", Z="cps",                    thr=0,    pub_r=0.339, pub_p=1.58e-3),
    list(id="CPS  ANLS-MCT4",      x="ANLS", y="MCT4",    Z="cps",                    thr=0,    pub_r=0.435, pub_p=NA),
    list(id="GW   MCT4-VATP",      x="MCT4", y="VATP_n6", Z=c("cps","gm_a","gm_n"),   thr=0,    pub_r=0.23,  pub_p=0.04),
    list(id="GW   MCT4-LAMP1",     x="MCT4", y="LAMP1_n", Z=c("cps","gm_a","gm_n"),   thr=0,    pub_r=0.31,  pub_p=0.004),
    list(id="GW   MCT4-LDHB",      x="MCT4", y="LDHB_n",  Z=c("cps","gm_a","gm_n"),   thr=0,    pub_r=0.16,  pub_p=0.14),
    list(id="CPS>=0.2 MCT4-VATP",  x="MCT4", y="VATP_n6", Z="cps",                    thr=0.2,  pub_r=0.461, pub_p=1.7e-5),
    list(id="CPS>=0.3 MCT4-VATP",  x="MCT4", y="VATP_n6", Z="cps",                    thr=0.3,  pub_r=0.479, pub_p=2.3e-5),
    list(id="own-gm  MCT4-VATP",   x="MCT4", y="VATP_n6", own=TRUE, thr=0,   pub_r=0.281, pub_p=9.5e-3),
    list(id="own-gm>=0.15",        x="MCT4", y="VATP_n6", own=TRUE, thr=0.15, pub_r=NA,    pub_p=NA),
    list(id="own-gm>=0.2",         x="MCT4", y="VATP_n6", own=TRUE, thr=0.2,  pub_r=0.261, pub_p=0.019),
    list(id="own-gm>=0.3",         x="MCT4", y="VATP_n6", own=TRUE, thr=0.3,  pub_r=0.271, pub_p=0.022)
  )

  cat(sprintf("%-22s %6s %6s %4s %4s  %11s %11s %11s  %s\n",
              "specification", "pub r", "obs r", "n", "df", "published p", "corrected p", "ratio", "flag"))
  res <- list(); nrep <- 0; nchg <- 0
  for (s in S) {
    d <- if (s$thr > 0) D[D$cps >= s$thr, ] else D
    if (isTRUE(s$own)) {
      a <- corr_own(d[[s$x]], d[[s$y]], d[, c("cps","gm_a")], d[, c("cps","gm_n")], "correct")
      b <- corr_own(d[[s$x]], d[[s$y]], d[, c("cps","gm_a")], d[, c("cps","gm_n")], "published")
    } else {
      Zm <- if (is.null(s$Z)) NULL else d[, s$Z, drop = FALSE]
      a <- corr(d[[s$x]], d[[s$y]], Zm, "correct")
      b <- corr(d[[s$x]], d[[s$y]], Zm, "published")
    }
    rep_ok <- is.na(s$pub_r) || abs(a[["r"]] - s$pub_r) <= TOL_R
    nrep <- nrep + rep_ok
    flag <- if (!rep_ok) "r MISMATCH"
            else if (!is.na(s$pub_p) && ((s$pub_p < MATERIAL) != (a[["p"]] < MATERIAL))) "CROSSES 0.05"
            else ""
    if (flag == "CROSSES 0.05") nchg <- nchg + 1
    cat(sprintf("%-22s %6s %+6.3f %4d %4d  %11.3g %11.3g %11.2f  %s\n",
                s$id, ifelse(is.na(s$pub_r),"  -",sprintf("%+6.3f",s$pub_r)), a[["r"]], a[["n"]], a[["df"]],
                b[["p"]], a[["p"]], a[["p"]]/b[["p"]], flag))
    res[[length(res)+1]] <- data.frame(spec=s$id, published_r=ifelse(is.na(s$pub_r),NA,s$pub_r), observed_r=round(a[["r"]],4),
      n=a[["n"]], k=a[["k"]], df_published=b[["df"]], df_correct=a[["df"]],
      p_published_recomputed=signif(b[["p"]],3), p_corrected=signif(a[["p"]],3),
      p_in_manuscript=s$pub_p, flag=flag, stringsAsFactors=FALSE)
  }
  cat(sprintf("\n  specifications reproduced on r: %d / %d\n", nrep, length(S)))
  cat(sprintf("  p-values that cross 0.05 after correction: %d\n", nchg))
  cat("  Any 'r MISMATCH' means the specification in this script does not match\n")
  cat("  the one used in the manuscript; that row's p-value is not to be used.\n")

  out <- do.call(rbind, res)
  write.csv(out, file.path(OUT, "partial_df_audit.csv"), row.names = FALSE)
  cat(sprintf("\nWritten to %s/partial_df_audit.csv\n", OUT))
  cat("=== paste the whole output back ===\n")
}

args <- commandArgs(trailingOnly = TRUE)
if (!selftest()) { cat("SELF-TEST FAILED - stopping before any file is read.\n"); quit(status = 1) }
if (length(args) && args[1] == "selftest") { cat("selftest mode: done\n"); quit(status = 0) }
main()
