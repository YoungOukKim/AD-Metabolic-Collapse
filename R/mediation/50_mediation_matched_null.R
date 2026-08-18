# =============================================================================
#  P2 revision  |  Script 05
#  Mediation specificity: detection- and effect-size-matched null
#
#  WHY THIS EXISTS
#  ---------------
#  Table 3 reports the fraction of each CPS -> outcome effect that is removed by
#  adjusting for astrocytic MCT4 (SLC16A3). A mediated fraction is NOT
#  interpretable on its own: adjusting for ANY variable that is strongly
#  associated with CPS will remove part of a weaker CPS -> outcome effect, purely
#  because of the difference in signal-to-noise ratio.
#
#  This script therefore benchmarks each mediated fraction against a null of
#  astrocytic genes matched to SLC16A3 on BOTH
#     (i)  detection rate            (SLC16A3 = 0.0669)
#     (ii) CPS effect size |t|       (SLC16A3 = 5.78)
#  Matching on |t| alone is not sufficient: it fills the null with low-detection
#  lncRNAs, whose noise attenuates them as mediators and biases the null
#  DOWNWARD, i.e. in favour of MCT4.
#
#  INPUT   : objects built by the upstream P2 scripts (E, DT, Cn, u, dcps, genes)
#  OUTPUT  : tables/Table3_mediation.csv                (manuscript Table 3)
#            tables/SuppTable_mediation_specificity.csv  (full test, incl. dropped row)
#            data/S_matched_null_controls.csv
#            data/S_matched_null_distribution.csv
#            data/Figure6_mediation_panel_data.csv
#
#  All paths are relative to the repository root. Run with the repo as wd.
# =============================================================================

## ---- 0. session guard -------------------------------------------------------
## A previous session had `sd` bound to a character string. Passing a masked
## symbol as FUN to apply() fails silently or fatally. Check before anything.
# Paths are RELATIVE to the repository root. Set the working directory there,
# or override with the environment variables SEAAD_H5AD / ROSMAP_ASTRO /
# ROSMAP_CLIN / P2_OUT_DIR. Raw data are not redistributable; see README.md.

.needed <- c("mean", "sd", "median", "quantile", "lm", "coef", "colMeans")
.bad <- Filter(function(f) {
  v <- get0(f, inherits = TRUE)
  is.null(v) || !is.function(v)
}, .needed)
if (length(.bad)) {
  stop("Masked base functions in this session: ", paste(.bad, collapse = ", "),
       "\nRestart R or rm() the offending objects before sourcing.")
}

## Project seed is 42 everywhere else. This block is the single documented
## exception: the bootstrap confidence intervals reported in Table 3
## (boot_lo / boot_hi) were generated with the seed below, so changing it
## would change those two columns and the deposited table would no longer be
## reproduced by this script. The mediated fractions, the matched-null
## membership and the empirical p values are deterministic and do not depend
## on it.
set.seed(20260713)
OUT_T <- "tables"; OUT_D <- "data"
dir.create(OUT_T, showWarnings = FALSE); dir.create(OUT_D, showWarnings = FALSE)

.req <- c("E", "DT", "u", "dcps", "genes")
.missing <- .req[!vapply(.req, exists, logical(1), inherits = TRUE)]
if (length(.missing)) {
  stop("Missing upstream objects: ", paste(.missing, collapse = ", "),
       "\nRun the P2 extraction script first (it builds E, DT, u, dcps, genes).")
}

## ---- 1. canonical definitions (hardcoded; do NOT read from workspace) -------
## The workspace carries several similarly-named variants (N_ACT vs ACT_NEUR,
## A_MITO vs MITO_ASTRO) that give different answers. These are the sets that
## reproduce the published Table 3 values.
MCT4_GENE <- "SLC16A3"

SETS <- list(
  "Neuronal activity genes"           = list(cmp = "exc",
    g = c("FOS", "NPAS4", "ARC", "EGR1", "JUNB"),                 published = 96.4),
  "Neuronal V-ATPase (6 subunits)"    = list(cmp = "exc",
    g = c("ATP6V1A", "ATP6V1B2", "ATP6V0A1", "ATP6V0C",
          "ATP6V0D1", "ATP6V1E1"),                                 published = 123.3),
  "Neuronal LAMP1"                    = list(cmp = "exc",
    g = c("LAMP1"),                                                published = 117.1),
  "Astrocytic mitochondrial ETC"      = list(cmp = "astro",
    g = c("NDUFS1", "NDUFV1", "NDUFA1", "UQCRC1",
          "COX4I1", "ATP5F1A", "ATP5F1B"),                         published = 45.7)
)

MATCH_TOL   <- c(0.10, 0.15, 0.20, 0.30, 0.50)  # widen until >= MIN_CTRL
MIN_CTRL    <- 200
N_BOOT      <- 2000

## ---- 2. helpers -------------------------------------------------------------
mod <- function(compartment, gs, idx) {
  g <- intersect(gs, colnames(E[[compartment]]))
  if (!length(g)) stop("no genes found for module")
  rowMeans(E[[compartment]][idx, g, drop = FALSE])
}

## mediated fraction: proportion of the CPS->y effect removed by adjusting for m
pct_med <- function(y, m, x) {
  o <- is.finite(y) & is.finite(m) & is.finite(x)
  if (sum(o) < 20) return(NA_real_)
  bt <- coef(lm(y[o] ~ x[o]))[2]
  if (!is.finite(bt) || abs(bt) < 1e-12) return(NA_real_)
  bd <- coef(lm(y[o] ~ x[o] + m[o]))[2]
  as.numeric((bt - bd) / bt * 100)
}

tstat <- function(v, x) {
  if (!is.finite(stats::sd(v, na.rm = TRUE)) ||
      stats::sd(v, na.rm = TRUE) == 0) return(NA_real_)
  abs(summary(lm(v ~ x))$coefficients[2, 3])
}

ax <- dcps[u]
M4 <- E[["astro"]][u, MCT4_GENE]

## ---- 3. G1 sanity: reproduce the six published percentages from raw ---------
## The script must not run if the pipeline is not the P2 pipeline.
pctv <- function(v, I_E, I_L) {
  e <- base::mean(v[I_E]); l <- base::mean(v[I_L])
  if (!is.finite(e) || abs(e) < 1e-9) NA_real_ else (l - e) / e * 100
}
if (all(vapply(c("BM", "I_E", "I_L"), exists, logical(1), inherits = TRUE))) {
  G1 <- data.frame(
    gene      = c("SLC16A3", "HK2", "LDHA", "PDK1", "SLC16A1", "SLC2A1"),
    published = c(-43.2, -35.2, -20.8, -16.5, -11.1, -7.4)
  )
  G1$recomputed <- vapply(G1$gene, function(g) round(pctv(BM[, g], I_E, I_L), 1), 0)
  G1$ok <- abs(G1$recomputed - G1$published) <= 1
  cat("\n[G1] P2 percentages reproduced from raw\n"); print(G1)
  if (!all(G1$ok)) stop("G1 FAILED - pipeline does not reproduce P2. STOP.")
  cat("  G1 >>> PASS\n")
} else {
  warning("G1 objects (BM, I_E, I_L) not in session; sanity check skipped.")
}

## ---- 4. build the matched null ---------------------------------------------
d4 <- base::mean(DT[["astro"]][u, MCT4_GENE])
t4 <- tstat(M4, ax)
cat(sprintf("\n[MATCH] %s  detection = %.4f   |t| = %.2f\n", MCT4_GENE, d4, t4))

Da <- colMeans(DT[["astro"]][u, , drop = FALSE])
Ta <- apply(E[["astro"]][u, , drop = FALSE], 2, function(z) tstat(z, ax))

ok <- integer(0)
for (w in MATCH_TOL) {
  ok <- which(is.finite(Da) & is.finite(Ta) &
              names(Da) != MCT4_GENE &
              abs(Da - d4) / d4 < w &
              abs(Ta - t4) / t4 < w)
  if (length(ok) >= MIN_CTRL) break
}
ctrl <- names(Da)[ok]
cat(sprintf("[MATCH] window = +/-%.0f%%   matched controls = %d\n",
            w * 100, length(ctrl)))
if (length(ctrl) < 50) stop("Too few matched controls; null is not usable.")

write.csv(
  data.frame(gene = ctrl,
             detection = round(Da[ok], 4),
             abs_t_CPS = round(Ta[ok], 3),
             row.names = NULL),
  file.path(OUT_D, "S_matched_null_controls.csv"), row.names = FALSE)

## ---- 5. observed vs null, per outcome --------------------------------------
res  <- list()
dist <- list()

for (nm in names(SETS)) {
  s <- SETS[[nm]]
  y <- mod(s$cmp, s$g, u)

  bt <- summary(lm(y ~ ax))$coefficients[2, c(1, 4)]
  bd <- summary(lm(y ~ ax + M4))$coefficients[2, c(1, 4)]
  obs <- pct_med(y, M4, ax)

  ## reproduce check against the published Table 3 value
  repro <- abs(obs - s$published) <= 0.5

  ## matched null
  nul <- vapply(ctrl, function(g) pct_med(y, E[["astro"]][u, g], ax), 0)
  nul <- nul[is.finite(nul)]

  pctile <- base::mean(nul < obs) * 100
  pemp   <- (sum(nul >= obs) + 1) / (length(nul) + 1)

  ## bootstrap CI on the observed mediated fraction (donor-level resampling)
  bs <- replicate(N_BOOT, {
    i <- sample(seq_along(y), replace = TRUE)
    pct_med(y[i], M4[i], ax[i])
  })
  bs <- bs[is.finite(bs)]
  ci <- stats::quantile(bs, c(0.025, 0.975), na.rm = TRUE)

  verdict <- if (pemp < 0.01) "RETAIN" else if (pemp < 0.05) "RETAIN (weak)" else "DROP"

  res[[nm]] <- data.frame(
    outcome        = nm,
    n              = length(u),
    beta_CPS       = round(bt[1], 4),
    p_CPS          = signif(bt[2], 3),
    beta_CPS_adj   = round(bd[1], 4),
    p_CPS_adj      = signif(bd[2], 3),
    pct_mediated   = round(obs, 1),
    published      = s$published,
    reproduced     = repro,
    boot_lo        = round(ci[1], 1),
    boot_hi        = round(ci[2], 1),
    null_n         = length(nul),
    null_median    = round(stats::median(nul), 1),
    null_p95       = round(stats::quantile(nul, 0.95), 1),
    null_max       = round(max(nul), 1),
    percentile     = round(pctile, 1),
    p_empirical    = signif(pemp, 3),
    verdict        = verdict,
    row.names      = NULL
  )
  dist[[nm]] <- data.frame(outcome = nm, gene = names(nul),
                           pct_mediated = round(nul, 2), row.names = NULL)

  cat(sprintf("  %-32s obs=%6.1f%% (pub %5.1f, repro=%s)  null_med=%6.1f%%  p95=%6.1f%%  pctile=%5.1f  p=%.4f  -> %s\n",
              nm, obs, s$published, ifelse(repro, "OK", "MISMATCH"),
              stats::median(nul), stats::quantile(nul, 0.95), pctile, pemp, verdict))
}

TAB <- do.call(rbind, res)
TAB <- TAB[order(TAB$p_empirical), ]          # strongest evidence first
DIST <- do.call(rbind, dist)

if (!all(TAB$reproduced)) {
  warning("At least one row does not reproduce the published Table 3 value. ",
          "Do not update the manuscript until this is resolved.")
}

## Table 3 as it appears in the manuscript: forward rows that clear the matched
## null (ordered by empirical p), then the reverse direction.
KEEP <- TAB[TAB$verdict != "DROP", ]
NAME <- c("Neuronal activity genes"        = "Neuronal activity genes",
          "Neuronal V-ATPase (6 subunits)" = "Neuronal V-ATPase (6 subunits, P2)",
          "Neuronal LAMP1"                 = "Neuronal LAMP1")
FWD <- data.frame(
  Direction = "MCT4 as mediator",
  Outcome   = unname(NAME[KEEP$outcome]),
  Mediator  = "Astrocytic MCT4",
  n = KEEP$n, beta_CPS = KEEP$beta_CPS, p = KEEP$p_CPS,
  beta_CPS_adj = KEEP$beta_CPS_adj, p_adj = KEEP$p_CPS_adj,
  pct_mediated = KEEP$pct_mediated, null_median = KEEP$null_median,
  null_p95 = KEEP$null_p95, percentile = KEEP$percentile,
  p_empirical = KEEP$p_empirical)
REV <- read.csv("tables/Table3_reverse.csv")   # unchanged reverse rows, from script 49
write.csv(rbind(FWD, REV[, names(FWD)]),
          file.path(OUT_T, "Table3_mediation.csv"), row.names = FALSE)
write.csv(TAB, file.path(OUT_T, "SuppTable_mediation_specificity.csv"), row.names = FALSE)
write.csv(DIST, file.path(OUT_D, "S_matched_null_distribution.csv"),  row.names = FALSE)

## ---- 6. Figure 6 mediation panel data  (was Figure 6 before renumbering) --------------------------------------
## The published panel plots pct_mediated as bars. Bars alone imply that a larger
## percentage is stronger evidence, which is false: the null median differs by
## outcome. The panel must overlay the matched-null median and 95th percentile.
FIG <- TAB[, c("outcome", "pct_mediated", "boot_lo", "boot_hi",
               "null_median", "null_p95", "percentile", "p_empirical", "verdict")]
write.csv(FIG, file.path(OUT_D, "Figure6_mediation_panel_data.csv"), row.names = FALSE)

cat("\n[DONE] tables/ and data/ updated\n")
print(TAB[, c("outcome", "pct_mediated", "null_median", "null_p95",
              "percentile", "p_empirical", "verdict")])

cat("\nDECISION RULE (fixed before running):\n")
cat("  p_empirical < 0.01  -> retain in Table 3\n")
cat("  0.01 - 0.05         -> retain, state specificity limit in Limitations\n")
cat("  > 0.05              -> remove the row from Table 3\n")
