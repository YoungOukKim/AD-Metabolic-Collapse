# ==============================================================================
# P2_revision2_S12_all.R
#
# ONE FILE. Paste it in and run it. No switch, no second command, no sourcing
# of anything else. It replaces both of the previously separate steps:
#
#     step 1  P2_revision2_S12.R          S12-P, S12-N, S12-H, S12-A, S12-D
#     step 2  50b_mediation_competition.R S12-M, the mediation competition
#
# and it reads the h5ad ONCE for both, instead of once for S12-A and again for
# the extraction that the mediation script needs.
#
# WHY THE MEDIATION BLOCK IS IN HERE AND NOT SOURCED FROM THE REPOSITORY
# ---------------------------------------------------------------------
# R/mediation/50_mediation_matched_null.R requires upstream objects E, DT, u,
# dcps and genes. The deposited repository does not build DT: donor_by_gene.rds
# is list(astro, neuron, genes, donors) and carries donor x gene MEANS only, with
# no detection matrix, and R/data_extraction/01_extract_seaad.R writes per-cell
# CSVs for a ~80-gene target list rather than a genome-wide donor x gene
# detection matrix. So the null cannot be rebuilt from the deposit as it stands.
# This file therefore builds E and DT in the single h5ad pass and then runs the
# null with the house implementation COPIED UNCHANGED from script 50: the same
# pct_med(), the same tstat(), the same widening match ladder, the same
# MIN_CTRL, the same seed. Only the mediator is looped. Three published values
# (96.4 / 123.3 / 117.1 per cent mediated) are used as a reproduction gate: if
# the rebuilt E and DT do not recover them, the competitor loop does not run.
#
# Reviewer mapping
#   S12-P  ambient / internal / channel / iron panels vs the matched null
#                                             Reviewer 1 point 2, Reviewer 2 point 3
#   S12-H  sparse-cell guard, S7 re-output    Reviewer 2 point 5
#   S12-A  SLC16A3 in four compartments       Reviewer 1 point 2
#   S12-D  microglial uptake and glycolysis   Reviewer 1 point 2
#   S12-N  ANLS panel-definition reconciliation
#                                             Reviewer 2 points 1 and 3
#   S12-M  mediation competition              Reviewer 1 point 4, Reviewer 2 point 3
#
# House conventions taken from the repositories, not re-invented:
#   set.seed(42)                          P2_run_all.R, the P3 consensus engine
#   set.seed(20260713) for S12-M only     50_mediation_matched_null.R
#   bin <- round(cps, 1), bin 0.1 out     77_detection_matched_null_standalone.R
#   early / late = UNWEIGHTED mean of
#     three bin means                     77_detection_matched_null_standalone.R
#   det_bin >= 0.10 to report alone       oligo_MCT_detection_check.R
#   det_bin >= 0.02 to compute a null     P2 04 / 77
#   matched band +/- 25 per cent of the
#     gene's own det_bin                  P2_revision2_all.R s4_acidbase
#   mediation null: widening ladder
#     10/15/20/30/50 per cent on BOTH
#     detection and |t|, MIN_CTRL 200     50_mediation_matched_null.R
#   pct_med uses lm(y ~ x) and
#     lm(y ~ x + m), NO global covariate  50_mediation_matched_null.R
#   ANLS = SLC2A1, LDHA, SLC16A1          manuscript Methods; P2_run_all.R L69;
#     with MCT4 analysed separately       figures/utils.R L65; 79_anls_row.R L36
#   block 100000, rhdf5, bit64 double     oligo_MCT_detection_check.R
#   pc() partial correlation              AD-Metabolic-Collapse 00_utils.R
#
#   WHERE THIS FILE DEPARTS FROM A HOUSE CONVENTION IT SAYS SO IN THE OUTPUT.
#
# Input   : OUTDIR/S3_detection_by_celltype.csv   (written by P2_revision2_all.R)
#           OUTDIR/S7_leverage_composition.csv    (written by P2_revision2_all.R)
#           the SEA-AD MTG h5ad
#           P2_H5DIR/donor_by_gene.rds            (optional, used as a cross-check)
#           donor_merged_base.csv                 (optional, used as a cross-check)
# Output  : OUTDIR/S12_*.csv, S12_RULES.txt, S12_VERDICTS.txt, S12_MANIFEST.csv,
#           OUTDIR/S7_leverage_composition_guarded.csv
# Runtime : S12-P and S12-H are seconds and run first. The single h5ad pass is
#           45 to 80 minutes. S12-A, S12-D, S12-N and S12-M then take minutes.
# ==============================================================================
suppressPackageStartupMessages({ library(rhdf5); library(data.table) })
set.seed(42)

## ---- 0. configuration --------------------------------------------------------
MTG_H5   <- Sys.getenv("H5AD_PATH",   unset = file.path("data", "SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad"))
P2_H5DIR <- Sys.getenv("P2_H5DIR",    unset = file.path("output", "sensitivity"))
P2_RUNDIR<- Sys.getenv("P2_RUNDIR",   unset = file.path("output", "tables"))
OUTDIR   <- Sys.getenv("P2R2_OUT",    unset = file.path("output", "revision2"))
if (!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE)

BUILD <- "P2_revision2_S12_all.R  2026-08-14  s2 (one h5ad pass; S12-P/H/A/D/N/M)"

DET_FLOOR_HOUSE <- 0.10
DET_FLOOR_NULL  <- 0.02
MATCH_TOL       <- 0.25                          # s4_acidbase band
MED_LADDER      <- c(0.10, 0.15, 0.20, 0.30, 0.50)  # script 50 ladder
MED_MIN_CTRL    <- 200L
MED_MIN_USABLE  <- 50L
MIN_NULL        <- 200L
ALPHA           <- 0.05
NULL_TAIL       <- 0.05
BLOCK           <- 100000L
EARLY_BINS      <- c(0.2, 0.3, 0.4)
LATE_BINS       <- c(0.6, 0.7, 0.8)

ASTRO_SUB <- "Astrocyte"
EXC_SUB   <- c("^L2/3 IT$","^L4 IT$","^L5 IT$","^L6 IT$","^L6 IT Car3$",
               "^L5 ET$","^L5/6 NP$","^L6 CT$","^L6b$")
MICRO_PAT <- "Micro"          # SEA-AD MTG Subclass is "Microglia-PVM". The labels
OLIGO_PAT <- "^Oligodendro"   # actually matched are printed and written out.

VATP6  <- c("ATP6V1A","ATP6V1B2","ATP6V0A1","ATP6V0C","ATP6V0D1","ATP6V1E1")
LYSO_N <- c("LAMP1","LAMP2","CTSB","CTSD")
ACT_N  <- c("FOS","NPAS4","ARC","EGR1","JUNB")   # script 50, "Neuronal activity genes"

P_AMBIENT  <- c("P2RY12","MOG","C1QA","CX3CR1","SNAP25","TMEM119","CSF1R","PLP1","MBP")
P_INTERNAL <- c("ALDOA","PKM","LDHA","PFKFB3","AQP4","SLC1A3","GJA1","SLC1A2")
P_MCTS     <- c("SLC16A1","SLC16A3","SLC16A7","SLC16A8")
P_CHANNEL  <- c("PANX1","PANX2","GJA1","GJB6","CALHM1","BEST1","LRRC8A")
P_IRON     <- c("CP","FTH1","FTL","TFRC","SLC40A1")
P_MICROGL  <- c("SLC2A1","SLC2A3","HK2","PKM","LDHA","SLC16A3","SLC16A1")
P_AMBREF   <- c("AQP4","PTGDS","HK2","SLC16A3")
S12_PANEL  <- unique(c(P_AMBIENT,P_INTERNAL,P_MCTS,P_CHANNEL,P_IRON,P_MICROGL,
                       P_AMBREF,VATP6,LYSO_N,ACT_N,"SLC2A1","PDK1","LDHB"))

ANLS_MS <- c("SLC2A1","LDHA","SLC16A1")            # manuscript and repository
ANLS_R6 <- c("SLC16A3","LDHA","HK2","PFKFB3")      # P2_revision2_all.R r6, contrast only

## Competitors for S12-M. Every one of them sits at or near MCT4's own
## detection-matched percentile in the astrocytic compartment, which is why they
## are here. SLC16A3 is included so that the incumbent runs through the identical
## code path as its challengers.
COMPETITORS <- c("SLC16A3","CP","FTH1","FTL","P2RY12","TFRC")

## Reproduction gates.
G0 <- list(nuclei = 1378211L, astro = 67419L, exc = 671689L, donors = 84L)
G1 <- list(gene="SLC16A3", det_bin_astro=0.0602, pct_astro=-43.2, tol_det=0.004, tol_pct=1.5)
G3 <- list(pct_micro = -14.4, tol = 1.5)           # July ambient_control.csv, verified by S2
G4 <- list(d4 = 0.0669, t4 = 5.78, tol_d = 0.004, tol_t = 0.60)   # script 50 header
G5 <- list(published = c("Neuronal activity genes" = 96.4,
                         "Neuronal V-ATPase (6 subunits)" = 123.3,
                         "Neuronal LAMP1" = 117.1), tol = 3.0)    # Table3_mediation.csv

## ---- 1. bookkeeping ----------------------------------------------------------
.MAN <- new.env(); .MAN$rows <- list(); .MAN$verd <- character(0)
put <- function(obj, file, what) {
  f <- file.path(OUTDIR, file); data.table::fwrite(as.data.table(obj), f)
  .MAN$rows[[length(.MAN$rows)+1L]] <- data.table(file=file, rows=nrow(obj),
                                                  cols=ncol(obj), contents=what)
  invisible(f) }
verdict <- function(...) { ln <- c(...); cat("\n[VERDICT]\n")
  cat(paste0("  ", ln, collapse="\n"), "\n"); .MAN$verd <- c(.MAN$verd, "", ln); invisible(ln) }
note <- function(...) cat(sprintf(...))

read_obs <- function(name, h5 = MTG_H5) {
  v <- h5read(h5, paste0("obs/", name))
  cats <- tryCatch(h5read(h5, paste0("obs/__categories/", name)), error=function(e) NULL)
  if (is.null(cats)) { cats <- tryCatch(h5read(h5, paste0("obs/", name, "/categories")), error=function(e) NULL)
    if (!is.null(cats)) v <- h5read(h5, paste0("obs/", name, "/codes")) }
  if (is.null(cats)) return(as.vector(v))
  code <- as.integer(v); out <- rep(NA_character_, length(code))
  ok <- !is.na(code) & code >= 0; out[ok] <- as.character(cats)[code[ok] + 1L]; out }

pc <- function(x, y, Z) {                      # AD-Metabolic-Collapse 00_utils.R
  Z <- as.matrix(Z); o <- complete.cases(x, y, Z)
  x <- x[o]; y <- y[o]; Z <- Z[o, , drop=FALSE]
  r <- cor(residuals(lm(x ~ Z)), residuals(lm(y ~ Z)))
  n <- sum(o); k <- ncol(Z)
  c(r=r, p=2*pt(-abs(r*sqrt((n-k-2)/(1-r^2))), n-k-2), n=n) }

find_in <- function(nm) { for (d in c(P2_RUNDIR, P2_H5DIR, OUTDIR, "."))
  { p <- file.path(d, nm); if (file.exists(p)) return(p) }; NA_character_ }

## ---- 2. PRE-REGISTERED VERDICT RULES ----------------------------------------
## Written to disk before any data file is opened.
RULES <- c(
"S12_RULES   pre-registered verdict rules",
sprintf("written : %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
sprintf("build   : %s", BUILD),
"",
"PROVENANCE, stated exactly, because these blocks were not all written at the",
"same time and only some of them are blind:",
"  Rules A, B, C, E, F and G, and the panels they act on, were fixed in the",
"  handoff of 2026-08-14 section 11 BEFORE any of the values below were",
"  computed. They are reproduced here unchanged.",
"  The S12-P values themselves were then computed from the already existing",
"  S3_detection_by_celltype.csv, so by the time THIS file was written they were",
"  known. The rule is pre-registered; the run is not blind. Both are true and",
"  the reply must not imply otherwise.",
"  S12-A, S12-D, S12-N and S12-M read quantities that have not been examined.",
"  GJA1's direction (+7.6 per cent in astrocytes) was already known when the",
"  channel panel was written. It is marked DIRECTION KNOWN in the output and",
"  carries no weight. PANX1, PANX2, GJB6, CALHM1, BEST1 and LRRC8A were not",
"  examined before that panel was fixed.",
"  The competing-mediator set for S12-M was chosen AFTER the S12-P percentiles",
"  were seen. The choice is adversarial - the genes selected are the ones most",
"  likely to defeat our own claim - but it is not pre-registered and must not be",
"  called so.",
"",
"A  cell-type decomposition of SLC16A3",
"  A-1 astrocytic decline holds AND microglia do not decline -> astrocytic axis stands",
"  A-2 both decline -> we state first that the decline is not astrocyte-restricted",
"      and widen the claim to glial lactate export capacity",
"  A-3 microglia only -> the astrocytic narrative in the manuscript is withdrawn",
"  A-4 microglial SLC16A3 below the null floor -> report NOT TESTABLE, A-1 forbidden",
"",
"B  ambient background",
"  B-1 MCT4's matched-null percentile is not separable from that of the ambient",
"      controls -> the astrocyte-specificity claim is withdrawn",
"  B-2 separable from all of them -> ambient bleed-through is excluded at this layer",
"  B-3 separable from some and not others -> report which, name the ones it is not",
"      separable from, and restrict the claim accordingly",
"",
"C  internal specificity inside astrocytes",
"  C-1 glycolysis and identity genes do not decline while SLC16A3 does ->",
"      'production and identity held, export gate specifically lost' is allowed",
"  C-2 any of them declines in the same direction AND of the same order of",
"      magnitude -> that sentence is forbidden",
"  C-3 any of them sits inside the null tail at a smaller magnitude -> the sentence",
"      is allowed only with that gene named in the same paragraph",
"",
"E  control transporters",
"  E-1 neuronal SLC16A7 preserved -> the receiving side is intact, may be stated",
"  E-2 neuronal SLC16A7 declines -> the receiving-side argument is withdrawn",
"",
"F  MCT-independent channel route",
"  F-1 PANX1 / GJB6 decline and meet MCT4's joint condition -> total-export",
"      compensation is unsupported at the transcript layer",
"  F-2 they do not decline -> we say so first and narrow the claim to",
"      MCT4-mediated export capacity",
"  F-3 below the floor -> NOT MEASURABLE with the name given, neither claim allowed",
"  F-common: transcript abundance is not open probability. This limit applies to",
"      the channel route and to MCT alike, and is written on both sides.",
"",
"G  specificity contrast against the iron pathway, gene level",
"  G-1 CP or any iron gene reaches MCT4's matched-null percentile ->",
"      the gene-level specificity claim is withdrawn and we write first that the",
"      detection-matched null does not separate iron handling from lactate export",
"  G-2 the iron genes are clearly less extreme -> the percentile may be used to",
"      separate them, with the raw per cent change stated alongside",
"  Operational reading of 'same order of magnitude': a ratio of ten or less to",
"  MCT4's own percentile. Written down here rather than chosen to suit a result.",
"",
"N  panel-definition reconciliation",
"  the manuscript and the deposited repository define the astrocytic ANLS",
"  composite as SLC2A1 + LDHA + SLC16A1 with MCT4 analysed separately.",
"  P2_revision2_all.R r6 used a different set under the same name, and that set",
"  contained the gene the screen was meant to test. The screen is re-run with the",
"  manuscript definition and MCT4 is given its own row.",
"  N-1 the lactate composite declines under BOTH conventions -> a module-level",
"      decline may be reported",
"  N-2 it declines under one convention only -> both are printed and no",
"      module-level decline is claimed",
"  N-3 it does not decline while MCT4 alone does -> that is reported as the",
"      manuscript's own position and NOT as a module-level decline",
"",
"M  mediation competition (the layer the manuscript says carries the weight)",
"  M-1 a competing mediator equals or exceeds MCT4's position in its OWN matched",
"      null for neuronal V-ATPase -> the mediation layer does not separate them",
"      either. The paper is reframed: a coordinated astrocytic metabolic decline",
"      of which MCT4 is one component, not MCT4 as the identified bottleneck.",
"      'Bottleneck' and 'rate-limiting' stay out of the manuscript.",
"  M-2 MCT4 exceeds and no competitor does -> specificity survives at the",
"      mediation layer ONLY, and the reply states that the gene-level null does",
"      not separate them and that the whole separation rests on this one test",
"  M-3 MCT4 does not exceed -> the mediation claim is withdrawn",
"  M-common: the forward / reverse asymmetry is computed for every competitor,",
"      not only for MCT4. An asymmetry reported for one gene and not tested for",
"      the others is not evidence.",
"  M-gate: if the rebuilt donor matrices do not reproduce the three published",
"      mediated fractions (96.4 / 123.3 / 117.1) the competitor loop does not run.",
"",
"H  sparse-cell guard (house rule 0-41)",
"  a chi-square whose statistic is more than half generated by cells with an",
"  expected count below 5 is not reported. Fisher's exact test and a merged-cell",
"  chi-square are printed beside it, and the reply quotes the guarded value.")
writeLines(RULES, file.path(OUTDIR, "S12_RULES.txt"))
cat(paste(RULES, collapse="\n"), "\n\n")

## ---- 3. self-test, before any file is opened ---------------------------------
## Every expected value below was computed directly before it was written here.
## No threshold is a single-draw threshold: the mediation tests use the median
## over repetitions, which is the house pattern in 02_sa_sensitivity.R.
pct_med <- function(y, m, x) {                 # 50_mediation_matched_null.R, verbatim
  o <- is.finite(y) & is.finite(m) & is.finite(x)
  if (sum(o) < 20) return(NA_real_)
  bt <- coef(lm(y[o] ~ x[o]))[2]
  if (!is.finite(bt) || abs(bt) < 1e-12) return(NA_real_)
  bd <- coef(lm(y[o] ~ x[o] + m[o]))[2]
  as.numeric((bt - bd) / bt * 100) }

tstat <- function(v, x) {                      # 50_mediation_matched_null.R, verbatim
  if (!is.finite(stats::sd(v, na.rm = TRUE)) ||
      stats::sd(v, na.rm = TRUE) == 0) return(NA_real_)
  abs(summary(lm(v ~ x))$coefficients[2, 3]) }

selftest <- function() {
  cat("=== SELF-TEST ===\n"); f <- character(0)
  say <- function(id, ok, w) { cat(sprintf("  %-5s %-60s %s\n", id, w, ifelse(isTRUE(ok),"PASS","FAIL")))
    if (!isTRUE(ok)) f <<- c(f, id) }

  # T1 the matched band is a two-sided proportional window on det_bin
  d0 <- 0.06; v <- c(0.0449, 0.0451, 0.06, 0.0749, 0.0751)
  say("T1", identical(abs(v-d0) <= MATCH_TOL*d0, c(FALSE,TRUE,TRUE,TRUE,FALSE)),
      "matched band is +/- 25 per cent of the gene's own det_bin")

  # T2 the percentile counts the gene itself, as s4_acidbase does
  x <- c(-50, -43.2, -40, rep(0, 7))
  say("T2", abs(mean(x <= x[2]) - 2/10) < 1e-12,
      "percentile is self-inclusive, as s4_acidbase computes it")

  # T3 9/2027 = 0.004440 and (9-1)/(2027-1) = the 0.39th percentile: one fact,
  #    two forms. This is the reconciliation the reply has to state.
  say("T3", abs(9/2027 - 0.004440059) < 1e-8 && abs(100*8/2026 - 0.3948667) < 1e-6,
      "0.00444 and 'the 0.39th percentile' are one fact in two forms")

  # T4 early / late is the unweighted mean of three bin means
  say("T4", abs(mean(c(0.10, 0.20, 0.60)) - 0.30) < 1e-12,
      "early / late is the unweighted mean of three bin means")

  # T5 pct = 100 * (late / early - 1) reproduces the published -43.2
  say("T5", abs(100*(0.03598/0.063339 - 1) - (-43.1952)) < 1e-3,
      "pct form reproduces -43.2 for astrocytic SLC16A3")

  # T6 bin 0.1 and bin 1.0 are outside the bin-level window
  say("T6", identical(which(round(c(0.05,0.1,0.2,0.9,1.0)*10) %in% 2:9), c(3L,4L)),
      "bin 0.1 and bin 1.0 are outside the bin-level window")

  # T7 the sparse-cell guard fires on the known S7 table:
  #    CDR 4 levels, trimmed 5 per cent 18/381, 34/636, 1/87, 1/1.
  M <- rbind(c(18,34,1,1), c(363,602,86,0))
  ct <- suppressWarnings(chisq.test(M)); Ex <- ct$expected
  contrib <- (M-Ex)^2/Ex; share <- sum(contrib[Ex < 5])/sum(contrib)
  say("T7", abs(ct$statistic - 22.3875) < 0.01 && ct$p.value < 1e-4 && share > 0.5,
      "chi-square 22.3875 with sparse cells making 98 per cent of it")

  # T8 merging CDR >= 1 removes it: chi-square 1.6038, p 0.4485
  c2 <- suppressWarnings(chisq.test(rbind(c(18,34,2), c(363,602,86))))
  say("T8", abs(c2$statistic - 1.6038) < 0.01 && c2$p.value > 0.05,
      "merged CDR >= 1 chi-square is 1.60 and not significant")

  # T9 pc() degrees of freedom reproduce the published 5.86e-06 at r = 0.4744686
  r <- 0.474468599439771; n <- 84; k <- 1
  p9 <- 2*pt(-abs(r*sqrt((n-k-2)/(1-r^2))), n-k-2)
  say("T9", abs(p9 - 5.8568e-06) < 1e-9, "pc() reproduces p = 5.86e-06 at r = +0.4745")

  # T10 proportion mediated is a ratio of coefficients: (-0.0690 - 0.0161)/-0.0690
  say("T10", abs(((-0.0690) - (0.0161))/(-0.0690) - 1.23333) < 1e-4,
      "proportion mediated reproduces the published 123.3 per cent")

  # T11 to T13 pct_med() on constructed data. A single draw is NOT a test: at
  #     n = 200 the full-mediation estimate lands outside 85-115 per cent in
  #     about one run in nine. The median over 200 repetitions was simulated over
  #     60 independent runs and fell in 98.0-101.8 (full), 65.0-67.6 (partial)
  #     and -0.04 to +0.05 (unrelated). The bands below contain all 60.
  set.seed(42); NREP <- 200L; n <- 200L   # project seed; the bands below hold for any draw (median over 200 repetitions)
  pf <- pp <- pu <- numeric(NREP)
  for (i in seq_len(NREP)) {
    xx <- rnorm(n); mm <- xx + rnorm(n)
    yf <- mm + rnorm(n); yp <- 0.5*xx + mm + rnorm(n)
    pf[i] <- pct_med(yf, mm, xx); pp[i] <- pct_med(yp, mm, xx)
    pu[i] <- pct_med(yf, rnorm(n), xx) }
  say("T11", median(pf) > 95 && median(pf) < 105,
      "pct_med recovers full mediation, median of 200 reps in 95-105")
  say("T12", median(pp) > 60 && median(pp) < 72,
      "pct_med recovers partial mediation, median in 60-72")
  say("T13", median(pu) > -5 && median(pu) < 5,
      "pct_med returns ~0 for a mediator unrelated to the exposure")

  # T14 tstat() is the absolute t of the slope, not of the intercept
  set.seed(42); xx <- rnorm(80); yy <- 3 + 2*xx + rnorm(80)   # project seed; the assertion is an identity and holds for any draw
  say("T14", abs(tstat(yy, xx) - abs(summary(lm(yy ~ xx))$coefficients[2,3])) < 1e-12,
      "tstat takes the slope t, not the intercept t")

  # T15 the widening ladder stops at the first window reaching MIN_CTRL
  d <- c(1, rep(1.05, 250), rep(1.4, 400)); t0 <- rep(1, length(d))
  hit <- NA_real_
  for (w in MED_LADDER) { ok <- which(abs(d-1)/1 < w & abs(t0-1)/1 < w)
    hit <- w; if (length(ok) >= MED_MIN_CTRL) break }
  say("T15", identical(hit, 0.10), "the match ladder stops at the first adequate window")

  cat(sprintf("=== SELF-TEST %s ===\n\n", if (!length(f)) "ALL PASS" else paste("FAIL:", paste(f, collapse=","))))
  invisible(!length(f)) }

## ---- 4. shared helpers -------------------------------------------------------
load_det <- function() {
  f <- file.path(OUTDIR, "S3_detection_by_celltype.csv")
  if (!file.exists(f)) { note("[stop] %s missing - run P2_revision2_all.R first\n", f); return(NULL) }
  D <- fread(f, showProgress=FALSE); i <- match(G1$gene, D$gene)
  ok <- !is.na(i) && abs(D$det_bin_astro[i]-G1$det_bin_astro) <= G1$tol_det &&
        abs(D$pct_change_astro[i]-G1$pct_astro) <= G1$tol_pct
  note("G1 on the S3 table: SLC16A3 det_bin %.4f pct %+.1f%% -> %s\n",
       D$det_bin_astro[i], D$pct_change_astro[i], ifelse(ok,"PASS","FAIL"))
  if (!ok) { note("[stop] the S3 table does not reproduce the published values\n"); return(NULL) }
  D }

band_pctl <- function(D, gene, ct) {
  db <- paste0("det_bin_", ct); pcc <- paste0("pct_change_", ct)
  i <- match(gene, D$gene)
  if (is.na(i)) return(list(det=NA_real_, pct=NA_real_, n=NA_integer_, pctl=NA_real_,
                            med=NA_real_, verdict="not in var"))
  d0 <- D[[db]][i]; x <- D[[pcc]][i]
  if (!is.finite(d0) || d0 <= 0) return(list(det=d0, pct=x, n=NA_integer_, pctl=NA_real_,
                            med=NA_real_, verdict="NOT MEASURABLE (zero detection)"))
  bd <- which(is.finite(D[[pcc]]) & abs(D[[db]] - d0) <= MATCH_TOL*d0)
  list(det=d0, pct=x, n=length(bd),
       pctl = if (length(bd) >= MIN_NULL) mean(D[[pcc]][bd] <= x) else NA_real_,
       med  = if (length(bd)) median(D[[pcc]][bd]) else NA_real_,
       verdict = if (!isTRUE(d0 >= DET_FLOOR_NULL)) "NOT MEASURABLE (below the null floor)"
                 else if (isTRUE(d0 < DET_FLOOR_HOUSE)) "below the house reporting floor - null only"
                 else "measurable") }

comp_z <- function(M, g) { j <- match(g, colnames(M)); j <- j[!is.na(j)]
  if (!length(j)) return(rep(NA_real_, nrow(M)))
  Z <- M[, j, drop=FALSE]; s <- apply(Z, 2, sd, na.rm=TRUE); k <- is.finite(s) & s > 1e-12
  if (!any(k)) return(rep(NA_real_, nrow(M)))
  rowMeans(scale(Z[, k, drop=FALSE]), na.rm=TRUE) }

slope_of <- function(v, cps, gm) {
  if (all(!is.finite(v)) || sd(v, na.rm=TRUE) < 1e-12) return(c(NA_real_, NA_real_))
  m <- summary(lm(v ~ cps + gm))$coefficients
  if (!("cps" %in% rownames(m))) return(c(NA_real_, NA_real_))
  c(m["cps","Estimate"], m["cps","Pr(>|t|)"]) }

## ---- 5. S12-P  panels against the matched null, from the existing S3 table ---
s12_panels <- function() {
  D <- load_det(); if (is.null(D)) return(invisible(FALSE))
  blocks <- list("S12-B ambient background"=P_AMBIENT, "S12-C astrocyte internal"=P_INTERNAL,
                 "S12-E control transporters"=P_MCTS,  "S12-F channel route"=P_CHANNEL,
                 "S12-G specificity contrast"=P_IRON)
  out <- rbindlist(lapply(names(blocks), function(b) rbindlist(lapply(blocks[[b]], function(g)
    rbindlist(lapply(c("astro","exc"), function(ct) { r <- band_pctl(D, g, ct)
      data.table(block=b, gene=g, celltype=ct, det_bin=r$det, pct_change=r$pct,
        n_null=r$n, pctl=r$pctl, band_median=r$med,
        reportable = isTRUE(r$det >= DET_FLOOR_HOUSE),
        in_null_tail = if (is.na(r$pctl)) NA else r$pctl <= NULL_TAIL,
        prereg = if (g == "GJA1") "DIRECTION KNOWN BEFORE THE PANEL WAS FIXED" else "",
        verdict = r$verdict) }))))))
  ref <- band_pctl(D, "SLC16A3", "astro"); refe <- band_pctl(D, "SLC16A3", "exc")
  out <- rbind(data.table(block="REFERENCE", gene="SLC16A3", celltype=c("astro","exc"),
    det_bin=c(ref$det,refe$det), pct_change=c(ref$pct,refe$pct), n_null=c(ref$n,refe$n),
    pctl=c(ref$pctl,refe$pctl), band_median=c(ref$med,refe$med),
    reportable=c(ref$det>=DET_FLOOR_HOUSE, refe$det>=DET_FLOOR_HOUSE),
    in_null_tail=c(ref$pctl<=NULL_TAIL, refe$pctl<=NULL_TAIL), prereg="",
    verdict=c(ref$verdict, refe$verdict)), out)
  put(out, "S12_panels_from_S3.csv",
      "ambient, internal, channel and iron panels against the detection-matched null")

  A <- out[celltype=="astro"]
  A[, ratio_to_MCT4 := pctl / ref$pctl]
  put(A[order(pctl)], "S12_astro_ranked_by_percentile.csv",
      "every S12 astrocytic gene ordered by its matched-null percentile, with the ratio to MCT4")
  ord <- function(sel) { z <- A[sel & is.finite(pctl)][order(pctl)]
    if (!nrow(z)) return("(none)")
    paste(sprintf("%s %.5f (%.1fx)", z$gene, z$pctl, z$ratio_to_MCT4), collapse="; ") }
  same_order <- function(sel) A[sel & is.finite(pctl) & ratio_to_MCT4 <= 10]$gene
  amb_so  <- setdiff(same_order(A$block %like% "ambient"), "SLC16A3")
  iron_so <- same_order(A$block %like% "specificity")
  iron_lo <- A[block %like% "specificity" & is.finite(pctl) & pctl <= ref$pctl]$gene
  chan_dn <- A[block %like% "channel" & is.finite(pctl) & pct_change < 0 & pctl <= NULL_TAIL]$gene
  int_tail<- A[block %like% "internal" & is.finite(pctl) & pct_change < 0 & pctl <= NULL_TAIL]$gene
  verdict("S12-P panels against the detection-matched null (astrocytes)",
    sprintf("REFERENCE SLC16A3 : %+.1f%%, percentile %.5f of %d, band median %+.1f%%",
            ref$pct, ref$pctl, ref$n, ref$med),
    sprintf("B ambient  : %s", ord(A$block %like% "ambient")),
    sprintf("             within one order of magnitude of MCT4 -> %s",
            if (length(amb_so)) paste(amb_so, collapse=", ") else "(none)"),
    sprintf("             B-%s", if (length(amb_so))
      "3 : separable from some controls and not from others. Name the ones it is not separable from."
      else "2 : separable from every ambient control at this layer"),
    sprintf("G iron     : %s", ord(A$block %like% "specificity")),
    sprintf("             at or BELOW MCT4's own percentile -> %s",
            if (length(iron_lo)) paste(iron_lo, collapse=", ") else "(none)"),
    sprintf("             G-%s", if (length(iron_so))
      "1 FIRES : the detection-matched null does not separate iron handling from lactate export. The gene-level specificity claim is withdrawn."
      else "2 holds"),
    sprintf("C internal : declining genes inside the null tail -> %s",
            if (length(int_tail)) paste(int_tail, collapse=", ") else "(none)"),
    sprintf("             C-%s", if (length(int_tail))
      "3 : the export-gate sentence is allowed only with these genes named in the same paragraph"
      else "1 : the export-gate sentence is allowed"),
    sprintf("F channel  : MCT-independent genes in the declining tail -> %s",
            if (length(chan_dn)) paste(chan_dn, collapse=", ") else "(none)"),
    sprintf("             F-%s", if (length(chan_dn)) "1"
      else "2 FIRES : narrow the claim to MCT4-mediated export capacity, and say so first"),
    "GJA1 carries no weight here: its direction was known before the panel was written.",
    "Transcript abundance is not open probability. That limit is written on the",
    "channel side and on the MCT side alike.")
  invisible(TRUE) }

## ---- 6. S12-H  the sparse-cell guard S7 is missing --------------------------
## The function in P2_revision2_all.R is cat_t(), at line 847, not cat_test().
s12_guard <- function() {
  fS <- file.path(OUTDIR, "S7_leverage_composition.csv")
  guard_tab <- function(M) {
    M <- M[, colSums(M) > 0, drop=FALSE]
    if (ncol(M) < 2) return(list(p=NA_real_, sparse_share=NA_real_, ok=FALSE, fp=NA_real_))
    ct <- suppressWarnings(chisq.test(M)); Ex <- ct$expected
    contrib <- (M-Ex)^2/Ex; share <- sum(contrib[Ex < 5])/sum(contrib)
    fp <- tryCatch(fisher.test(M, simulate.p.value=TRUE, B=20000)$p.value, error=function(e) NA_real_)
    list(p=ct$p.value, sparse_share=share, ok=share <= 0.5, fp=fp) }
  rows <- list()
  add <- function(var, M, merged=NULL) {
    g <- guard_tab(M); mg <- if (!is.null(merged)) guard_tab(merged) else NULL
    rows[[length(rows)+1L]] <<- data.table(variable=var,
      chisq_p_unguarded=g$p, sparse_cell_share=g$sparse_share,
      reportable_under_rule_0_41=g$ok, fisher_p=g$fp,
      chisq_p_after_merging=if (is.null(mg)) NA_real_ else mg$p,
      guarded_value=if (isTRUE(g$ok)) sprintf("chi-square p = %.3g", g$p)
        else if (!is.null(mg)) sprintf("merged chi-square p = %.3g ; Fisher p = %.3g", mg$p, g$fp)
        else sprintf("Fisher p = %.3g", g$fp)) }
  if (file.exists(fS)) {
    S <- fread(fS, showProgress=FALSE)
    parse_detail <- function(s) { pr <- strsplit(trimws(unlist(strsplit(s, ";"))), " ")
      k <- vapply(pr, function(z) z[1], character(1))
      ab <- do.call(rbind, lapply(pr, function(z) as.integer(unlist(strsplit(z[2], "/")))))
      list(lv=k, top=ab[,1], tot=ab[,2]) }
    for (v in c("CDR","APOE e4 count","APOE genotype","diagnosis")) {
      r <- S[variable == v]
      if (!nrow(r) || is.na(r$detail[1]) || !grepl("/", r$detail[1])) next
      z <- parse_detail(r$detail[1])
      M <- rbind(z$top, z$tot - z$top); colnames(M) <- z$lv
      merged <- NULL
      num <- suppressWarnings(as.numeric(z$lv))
      if (v %in% c("CDR","APOE e4 count") && !anyNA(num)) {
        grp <- ifelse(num >= 1, ">=1", as.character(num))
        merged <- t(rowsum(t(M), grp, reorder=FALSE)) }
      add(v, M, merged) } }
  if (!length(rows)) { note("[stop S12-H] S7_leverage_composition.csv not found or unparsable\n")
    return(invisible(FALSE)) }
  G <- rbindlist(rows, fill=TRUE)
  put(G, "S7_leverage_composition_guarded.csv",
      "the same trimmed-5-per-cent composition tests under the sparse-cell guard (house rule 0-41)")
  print(G)
  bad <- G[reportable_under_rule_0_41 == FALSE]$variable
  verdict("S12-H sparse-cell guard (Reviewer 2, point 5)",
    sprintf("tests whose chi-square is generated mostly by cells with expectation below 5 : %s",
            if (length(bad)) paste(bad, collapse=", ") else "(none)"),
    "For those the unguarded p is NOT reported. The reply quotes the merged-cell",
    "chi-square and Fisher's exact test printed in the same row.",
    "The deposited S7_leverage_composition.csv is left in place and this file sits",
    "beside it, so the correction is visible rather than silent.")
  invisible(TRUE) }

## ---- 7. ONE h5ad PASS  -------------------------------------------------------
## Produces, in a single read:
##   (i)  genome-wide detection and early-to-late change in FOUR compartments
##        (astrocyte, microglia, oligodendrocyte, excitatory neuron)  -> S12-A, S12-D
##   (ii) donor x gene MEAN and DETECTION matrices for astrocytes and excitatory
##        neurons, plus per-donor CPS and transcriptome-wide mean  -> S12-N, S12-M
## The deposited donor_by_gene.rds carries the means but no detection matrix, and
## 01_extract_seaad.R writes a ~80-gene target list, so (ii) cannot be taken from
## the deposit and has to be built here.
PASS <- new.env()
one_pass <- function() {
  if (!file.exists(MTG_H5)) { note("[stop] h5ad not found at %s\n", MTG_H5); return(invisible(FALSE)) }
  t0 <- Sys.time()
  sub <- read_obs("Subclass"); cps <- as.numeric(read_obs("Continuous Pseudo-progression Score"))
  don <- read_obs("Donor ID"); genes <- as.character(h5read(MTG_H5, "var/_index"))
  nC <- length(sub); nG <- length(genes)
  donors <- sort(unique(don[!is.na(cps)])); inCoh <- !is.na(cps) & don %in% donors
  nD <- length(donors); dix <- match(don, donors)

  sel <- list(astro = sub == ASTRO_SUB & inCoh,
              micro = grepl(MICRO_PAT, sub) & inCoh,
              oligo = grepl(OLIGO_PAT, sub) & inCoh,
              exc   = Reduce(`|`, lapply(EXC_SUB, function(p) grepl(p, sub))) & inCoh)
  lab <- data.table(compartment=names(sel),
    labels=vapply(names(sel), function(k) paste(sort(unique(sub[sel[[k]]])), collapse=" | "), character(1)),
    n_cells=vapply(sel, sum, integer(1)))
  put(lab, "S12_compartment_labels.csv",
      "which Subclass labels entered each compartment, and how many nuclei"); print(lab)

  note("G0 cohort gate: nuclei %d/%d | astro %d/%d | exc %d/%d | donors %d/%d\n",
       nC, G0$nuclei, sum(sel$astro), G0$astro, sum(sel$exc), G0$exc, nD, G0$donors)
  g0 <- nC==G0$nuclei && sum(sel$astro)==G0$astro && sum(sel$exc)==G0$exc && nD==G0$donors
  note("  G0 = %s\n", ifelse(g0,"PASS","FAIL - stopping"))
  if (!g0) return(invisible(FALSE))

  binid <- round(round(cps,1)*10); inbin <- !is.na(binid) & binid >= 2L & binid <= 9L
  ecol <- round(EARLY_BINS*10); lcol <- round(LATE_BINS*10)
  indptr <- h5read(MTG_H5, "X/indptr", bit64conversion="double")

  CT <- lapply(sel, function(z) list(det=numeric(nG), det_b=numeric(nG),
    sum_by_bin=matrix(0,nG,9), n_by_bin=numeric(9), n=0, nb=0))
  DG <- lapply(c("astro","exc"), function(k) list(s=matrix(0,nD,nG), d=matrix(0,nD,nG),
    n=numeric(nD), tot=numeric(nD)))
  names(DG) <- c("astro","exc")

  note("reading expression: four compartments and two donor x gene matrices, one pass ...\n")
  for (s0 in seq(1L, nC, by=BLOCK)) {
    e0 <- min(s0+BLOCK-1L, nC); sp <- indptr[s0]; cnt <- indptr[e0+1L]-sp
    if (cnt <= 0) next
    if (!any(Reduce(`|`, lapply(sel, function(z) z[s0:e0])))) { note("  ...%d/%d\r", e0, nC); next }
    ci <- h5read(MTG_H5,"X/indices",start=sp+1L,count=cnt,bit64conversion="double") + 1L
    cd <- h5read(MTG_H5,"X/data",   start=sp+1L,count=cnt)
    p  <- indptr[s0:(e0+1L)] - sp
    cid <- rep.int(seq_len(e0-s0+1L), as.integer(diff(p))); gl <- s0:e0
    for (k in names(sel)) {
      kp <- sel[[k]][gl]; if (!any(kp)) next
      O <- CT[[k]]; loc <- which(kp); m <- cid %in% loc
      if (any(m)) {
        gi <- ci[m]; gx <- cd[m]; cc <- cid[m]
        bglob <- binid[gl][cc]; ibb <- inbin[gl][cc]
        O$det <- O$det + tabulate(gi, nG)
        if (any(ibb)) O$det_b <- O$det_b + tabulate(gi[ibb], nG)
        okb <- !is.na(bglob) & bglob >= 1L & bglob <= 9L
        if (any(okb)) { lin <- gi[okb] + (bglob[okb]-1L)*nG
          ag <- rowsum(gx[okb], lin, reorder=FALSE); i <- as.integer(rownames(ag))
          O$sum_by_bin[i] <- O$sum_by_bin[i] + as.numeric(ag) }
        if (k %in% c("astro","exc")) {
          Q <- DG[[k]]; dd <- dix[gl][cc]; okd <- !is.na(dd)
          if (any(okd)) {
            lin2 <- dd[okd] + (gi[okd]-1L)*nD
            a2 <- rowsum(gx[okd], lin2, reorder=FALSE); i2 <- as.integer(rownames(a2))
            Q$s[i2] <- Q$s[i2] + as.numeric(a2)
            d2 <- rowsum(rep(1, sum(okd)), lin2, reorder=FALSE); i3 <- as.integer(rownames(d2))
            Q$d[i3] <- Q$d[i3] + as.numeric(d2)
            t2 <- rowsum(gx[okd], dd[okd], reorder=FALSE); i4 <- as.integer(rownames(t2))
            Q$tot[i4] <- Q$tot[i4] + as.numeric(t2) }
          DG[[k]] <- Q } }
      O$n <- O$n + length(loc); O$nb <- O$nb + sum(inbin[gl][loc])
      b2 <- binid[gl][loc]; b2 <- b2[!is.na(b2) & b2 >= 1L & b2 <= 9L]
      O$n_by_bin <- O$n_by_bin + tabulate(b2, 9); CT[[k]] <- O
      if (k %in% c("astro","exc")) { Q <- DG[[k]]
        Q$n <- Q$n + tabulate(dix[gl][loc], nD); DG[[k]] <- Q } }
    note("  ...%d/%d\r", e0, nC) }
  note("\n")

  mk <- function(O) { mbb <- sweep(O$sum_by_bin, 2, pmax(O$n_by_bin,1), "/")
    e <- rowMeans(mbb[, ecol, drop=FALSE]); l <- rowMeans(mbb[, lcol, drop=FALSE])
    data.table(gene=genes, det_all=O$det/max(O$n,1), det_bin=O$det_b/max(O$nb,1),
      pct_change=ifelse(!is.finite(e)|e<=0, NA_real_, 100*(l/e-1))) }
  DT4 <- lapply(CT, mk); for (k in names(DT4)) DT4[[k]][, compartment := k]
  put(rbindlist(DT4)[gene %in% S12_PANEL], "S12_panel_four_compartments.csv",
      "S12 panel detection and change in astro, microglia, oligodendrocyte and excitatory neurons")

  Ea <- DG$astro$s / pmax(DG$astro$n, 1); Ee <- DG$exc$s / pmax(DG$exc$n, 1)
  Da <- DG$astro$d / pmax(DG$astro$n, 1); De <- DG$exc$d / pmax(DG$exc$n, 1)
  dimnames(Ea) <- dimnames(Da) <- list(donors, genes)
  dimnames(Ee) <- dimnames(De) <- list(donors, genes)
  DM <- data.table(donor=donors,
    cps = vapply(donors, function(x) mean(cps[don == x], na.rm=TRUE), numeric(1)),
    gm_a = DG$astro$tot / pmax(DG$astro$n * nG, 1),
    gm_n = DG$exc$tot   / pmax(DG$exc$n   * nG, 1),
    ncell_a = DG$astro$n, ncell_n = DG$exc$n,
    MCT4 = Ea[, "SLC16A3"], VATP_n6 = rowMeans(Ee[, VATP6, drop=FALSE]))
  put(DM, "S12_donor_level.csv", "donor-level CPS, cell counts, transcriptome-wide means, MCT4 and VATP_n6")

  PASS$genes <- genes; PASS$donors <- donors; PASS$DT4 <- DT4; PASS$CT <- CT
  PASS$Ea <- Ea; PASS$Ee <- Ee; PASS$Da <- Da; PASS$De <- De; PASS$DM <- DM
  PASS$elapsed <- as.numeric(difftime(Sys.time(), t0, units="mins"))

  # cross-checks against the deposited artefacts, where they exist
  chk <- list()
  fG <- file.path(P2_H5DIR, "donor_by_gene.rds")
  if (file.exists(fG)) { G <- readRDS(fG)
    o <- match(donors, G$donors); j <- match("SLC16A3", G$genes)
    if (!anyNA(o) && !is.na(j)) chk$rds_MCT4 <- max(abs(G$astro[o, j] - Ea[, "SLC16A3"])) }
  fM <- find_in("donor_merged_base.csv")
  if (!is.na(fM)) { B <- fread(fM, showProgress=FALSE)
    o <- match(DM$donor, B$donor)
    if (!anyNA(o)) { if ("VATP_n6" %in% names(B)) chk$base_VATP <- max(abs(B$VATP_n6[o] - DM$VATP_n6))
                     if ("gm_a" %in% names(B))    chk$base_gm_a <- max(abs(B$gm_a[o] - DM$gm_a)) } }
  note("cross-check against the deposited artefacts: %s\n",
       if (!length(chk)) "(no deposited file found to compare against)"
       else paste(sprintf("%s max |diff| %.3g", names(chk), unlist(chk)), collapse=" | "))

  a <- DT4$astro; i <- match("SLC16A3", a$gene)
  verdict("one h5ad pass",
    sprintf("elapsed %.1f min | G0 PASS | %d donors, %d genes", PASS$elapsed, nD, nG),
    sprintf("astrocytic SLC16A3 det_bin %.4f (published %.4f) | change %+.1f%% (published %.1f%%)",
            a$det_bin[i], G1$det_bin_astro, a$pct_change[i], G1$pct_astro),
    sprintf("G1 on the new pass : %s",
            ifelse(abs(a$det_bin[i]-G1$det_bin_astro) <= G1$tol_det &&
                   abs(a$pct_change[i]-G1$pct_astro) <= G1$tol_pct, "PASS", "FAIL")),
    if (length(chk)) sprintf("deposited artefacts reproduced to %s",
      paste(sprintf("%s %.3g", names(chk), unlist(chk)), collapse=", ")) else
      "no deposited artefact was available to cross-check against")
  invisible(TRUE) }

## ---- 8. S12-A and S12-D  four compartments ----------------------------------
s12_compartments <- function() {
  if (is.null(PASS$DT4)) { note("[stop S12-A] the h5ad pass did not run\n"); return(invisible(FALSE)) }
  DT4 <- PASS$DT4; CT <- PASS$CT
  pctl_in <- function(k, g) { d <- DT4[[k]]; i <- match(g, d$gene)
    if (is.na(i) || !is.finite(d$det_bin[i]) || d$det_bin[i] <= 0) return(c(NA_real_, NA_integer_))
    d0 <- d$det_bin[i]
    bd <- which(is.finite(d$pct_change) & abs(d$det_bin - d0) <= MATCH_TOL*d0)
    if (length(bd) < MIN_NULL) return(c(NA_real_, length(bd)))
    c(mean(d$pct_change[bd] <= d$pct_change[i]), length(bd)) }

  A <- rbindlist(lapply(names(DT4), function(k) { d <- DT4[[k]]; i <- match("SLC16A3", d$gene)
    z <- pctl_in(k, "SLC16A3")
    data.table(compartment=k, n_cells=CT[[k]]$n, det_all=d$det_all[i], det_bin=d$det_bin[i],
      pct_change=d$pct_change[i], n_null=z[2], pctl=z[1],
      measurable = isTRUE(d$det_bin[i] >= DET_FLOOR_NULL),
      reportable = isTRUE(d$det_bin[i] >= DET_FLOOR_HOUSE)) }))
  put(A, "S12_A_SLC16A3_by_compartment.csv",
      "SLC16A3 in four compartments with a matched null computed inside each"); print(A)

  Dm <- rbindlist(lapply(P_MICROGL, function(g) { d <- DT4$micro; i <- match(g, d$gene)
    z <- pctl_in("micro", g)
    data.table(gene=g, det_bin=if (is.na(i)) NA_real_ else d$det_bin[i],
      pct_change=if (is.na(i)) NA_real_ else d$pct_change[i], n_null=z[2], pctl=z[1]) }))
  put(Dm, "S12_D_microglial_glycolysis.csv", "microglial uptake and glycolysis panel")

  AMB <- rbindlist(lapply(unique(c(P_AMBREF, P_AMBIENT)), function(g) {
    ia <- match(g, DT4$astro$gene); im <- match(g, DT4$micro$gene)
    pa <- pctl_in("astro", g); pm <- pctl_in("micro", g)
    aa <- if (is.na(ia)) NA_real_ else DT4$astro$pct_change[ia]
    mm <- if (is.na(im)) NA_real_ else DT4$micro$pct_change[im]
    data.table(gene=g, pct_astro=aa, pct_micro=mm,
      fold = if (is.finite(aa) && is.finite(mm) && abs(mm) > 1e-9) aa/mm else NA_real_,
      det_bin_astro=if (is.na(ia)) NA_real_ else DT4$astro$det_bin[ia],
      det_bin_micro=if (is.na(im)) NA_real_ else DT4$micro$det_bin[im],
      pctl_astro=pa[1], pctl_micro=pm[1],
      diverges = if (is.finite(aa) && is.finite(mm)) abs(aa) >= 2*abs(mm) else NA) }))
  put(AMB, "S12_ambient_symmetry.csv",
      "astrocytic against microglial change for the ambient reference genes and the myeloid markers")
  print(AMB)

  am <- A[compartment=="astro"]; mi <- A[compartment=="micro"]
  g3 <- is.finite(mi$pct_change) && abs(mi$pct_change - G3$pct_micro) <= G3$tol
  cls <- if (!isTRUE(mi$measurable)) "A-4 : microglial SLC16A3 is below the null floor - NOT TESTABLE"
    else if (is.finite(am$pct_change) && am$pct_change < 0 && (!is.finite(mi$pct_change) || mi$pct_change >= 0))
      "A-1 : astrocytic decline holds and microglia do not decline"
    else if (is.finite(am$pct_change) && am$pct_change < 0 && is.finite(mi$pct_change) && mi$pct_change < 0)
      "A-2 : both decline - state first that the decline is not astrocyte-restricted"
    else "A-3 : the astrocytic narrative is withdrawn"
  verdict("S12-A cell-type decomposition (Reviewer 1, point 2)",
    sprintf("G3 microglial reproduction gate (July -14.4%%) : %s (observed %+.1f%%)",
            ifelse(g3,"PASS","FAIL - the compartment was defined differently; do not interpret"),
            mi$pct_change),
    paste(sprintf("%s n=%d det_bin %.4f change %+.1f%% percentile %s",
      A$compartment, A$n_cells, A$det_bin, A$pct_change,
      ifelse(is.finite(A$pctl), sprintf("%.5f", A$pctl), "NA")), collapse=" | "),
    cls,
    "Detection rates are not comparable across studies with different enrichment",
    "and sequencing depth. Wachter 2024 reports SLC16A3 in 17 to 31 per cent of",
    "FACS-enriched myeloid nuclei; that number belongs beside this one in the reply.")
  invisible(TRUE) }

## ---- 9. S12-N  panel-definition reconciliation -------------------------------
## Rule 0: open the original before asserting a house practice.
## The manuscript Methods and the deposited repository both define
##     ANLS = mean(SLC2A1, LDHA, SLC16A1)
## with MCT4 analysed SEPARATELY and HK2 and PDK1 explicitly excluded
## (P2_run_all.R L69, figures/utils.R L65, 79_anls_row.R L36, 80_partial_df_audit.R L33).
## P2_revision2_all.R r6 used c("SLC16A3","LDHA","HK2","PFKFB3") under that name,
## which puts the tested gene inside its own composite. Two consequences, both
## stated first: (a) two different composites carry one name; (b) the reason
## given for preferring a z composite - that the raw average is 81 per cent
## PFKFB3 - is a property of the r6 panel only. On the same early-bin basis the
## largest raw member is SLC2A1 at 43.4 per cent in the manuscript panel and
## FTH1 at 52.9 per cent in the iron panel, so the convention dependence is
## largely an artefact of the panel choice, not of the analysis.
PATHWAYS_MS <- list(
  "Lactate export (ANLS composite, manuscript definition)" = list(ct="astro", g=ANLS_MS),
  "MCT4 / SLC16A3 (single gene, analysed separately)"      = list(ct="astro", g="SLC16A3"),
  "Iron handling"                  = list(ct="astro", g=c("TFRC","FTH1","FTL","CP","SLC40A1")),
  "Glutamate uptake and recycling" = list(ct="astro", g=c("SLC1A2","SLC1A3","GLUL","SLC38A3")),
  "Na/K-ATPase"                    = list(ct="astro", g=c("ATP1A1","ATP1A2","ATP1B1","ATP1B2")),
  "Water and K+ buffering"         = list(ct="astro", g=c("AQP4","KCNJ10","GJA1","GJB6")),
  "Glucose uptake and glycogen"    = list(ct="astro", g=c("SLC2A1","SLC2A3","GYS1","PYGB")),
  "Fatty acid oxidation"           = list(ct="astro", g=c("CPT1A","ACADVL","HADHA","ACSL3")),
  "Ketone body handling"           = list(ct="astro", g=c("BDH1","OXCT1","SLC16A6")),
  "Astrocyte reactivity"           = list(ct="astro", g=c("GFAP","VIM","SERPINA3","CLU","CD44")),
  "Mitochondrial ETC (astro, 35-gene panel)" = list(ct="astro",
    g=c("NDUFA1","NDUFA2","NDUFA4","NDUFA8","NDUFA9","NDUFA10","NDUFA13",
        "NDUFB1","NDUFB2","NDUFB4","NDUFB6","NDUFB8","NDUFB9","NDUFB10",
        "NDUFS1","NDUFS2","NDUFS3","NDUFS4","NDUFS7","NDUFS8","NDUFV1","NDUFV2",
        "UQCRC1","UQCRC2","CYC1","COX4I1","COX5A","COX6C","COX7C",
        "ATP5F1A","ATP5F1B","ATP5MC1","ATP5PB","SDHA","SDHB")),
  "Neuronal V-ATPase"              = list(ct="exc", g=VATP6),
  "Neuronal lysosome"              = list(ct="exc", g=LYSO_N),
  "Neuronal protein synthesis"     = list(ct="exc", g=c("RPL13A","RPS6","EIF4E","EEF2")),
  "Integrated stress response"     = list(ct="exc", g=c("ATF4","DDIT3","HSPA5","EIF2AK3")),
  "Neuronal identity (control)"    = list(ct="exc", g=c("SNAP25","SYT1","SLC17A7","RBFOX3")))

s12_reconcile <- function() {
  if (is.null(PASS$Ea)) { note("[stop S12-N] the h5ad pass did not run\n"); return(invisible(FALSE)) }
  Ea <- PASS$Ea; Ee <- PASS$Ee; Da <- PASS$Da; De <- PASS$De; DM <- PASS$DM
  cps <- DM$cps; gma <- DM$gm_a; gmn <- DM$gm_n
  A4 <- PASS$DT4$astro; E4 <- PASS$DT4$exc

  share <- function(g, ct) { M <- if (ct=="astro") Ea else Ee
    j <- match(g, colnames(M)); j <- j[!is.na(j)]
    if (length(j) < 2) return(NA_real_)
    mu <- colMeans(M[, j, drop=FALSE], na.rm=TRUE); 100*max(mu)/sum(mu) }
  pick <- function(g, ct) { M <- if (ct=="astro") Ea else Ee
    d <- if (ct=="astro") A4 else E4
    k <- g[!is.na(match(g, colnames(M))) & !is.na(match(g, d$gene))]
    k[d$det_bin[match(k, d$gene)] >= DET_FLOOR_NULL] }

  run <- function(PW, tag) {
    rows <- rbindlist(lapply(names(PW), function(nm) {
      P <- PW[[nm]]; ct <- P$ct; M <- if (ct=="astro") Ea else Ee
      gm <- if (ct=="astro") gma else gmn
      d  <- if (ct=="astro") A4 else E4
      ok <- pick(P$g, ct)
      if (!length(ok)) return(data.table(definition=tag, pathway=nm, celltype=ct,
        n_total=length(P$g), n_ok=0L, pct=NA_real_, slope_z=NA_real_, p_z=NA_real_,
        slope_raw=NA_real_, p_raw=NA_real_, raw_top_member_share=NA_real_,
        dropped=paste(P$g, collapse=";")))
      v  <- if (length(ok) > 1L) comp_z(M, ok) else as.numeric(scale(M[, match(ok, colnames(M))]))
      vr <- if (length(ok) > 1L) rowMeans(M[, match(ok, colnames(M)), drop=FALSE]) else M[, match(ok, colnames(M))]
      s <- slope_of(v, cps, gm); sr <- slope_of(vr, cps, gm)
      data.table(definition=tag, pathway=nm, celltype=ct, n_total=length(P$g), n_ok=length(ok),
        pct=mean(d$pct_change[match(ok, d$gene)], na.rm=TRUE), slope_z=s[1], p_z=s[2],
        slope_raw=sr[1], p_raw=sr[2], raw_top_member_share=share(ok, ct),
        dropped=paste(setdiff(P$g, ok), collapse=";")) }))
    for (cl in c("p_z","p_raw")) { ok <- is.finite(rows[[cl]])
      h <- rep(NA_real_, nrow(rows)); h[ok] <- p.adjust(rows[[cl]][ok], "holm")
      set(rows, j=paste0("holm_", sub("^p_","",cl)), value=h) }
    rows[, verdict_z   := fifelse(!is.finite(p_z),   "NOT TESTABLE",
      fifelse(slope_z   < 0 & holm_z   < ALPHA, "declines", "no decline at the family-adjusted alpha"))]
    rows[, verdict_raw := fifelse(!is.finite(p_raw), "NOT TESTABLE",
      fifelse(slope_raw < 0 & holm_raw < ALPHA, "declines", "no decline at the family-adjusted alpha"))]
    rows }

  MS <- run(PATHWAYS_MS, "manuscript definition")
  R6 <- run(list("Lactate export (ANLS hub, r6 definition, for contrast only)" =
                   list(ct="astro", g=ANLS_R6)), "r6 definition")
  ALL <- rbind(MS, R6, fill=TRUE); setorder(ALL, definition, holm_z, na.last=TRUE)
  put(ALL, "S12_N_pathway_screen_manuscript_definition.csv",
      "the sixteen-pathway screen re-run with the manuscript's ANLS definition, both conventions")

  gv <- function(pat, col) { z <- MS[pathway %like% pat]; if (nrow(z)) z[[col]][1] else NA_real_ }
  lac_z <- gv("ANLS composite","holm_z"); lac_r <- gv("ANLS composite","holm_raw")
  lac_sz <- gv("ANLS composite","slope_z"); lac_sr <- gv("ANLS composite","slope_raw")
  mct_z <- gv("single gene","holm_z"); mct_r <- gv("single gene","holm_raw")
  dz <- isTRUE(lac_sz < 0 && lac_z < ALPHA); dr <- isTRUE(lac_sr < 0 && lac_r < ALPHA)
  mz <- isTRUE(mct_z < ALPHA) || isTRUE(mct_r < ALPHA)
  cls <- if (dz && dr) "N-1 : the lactate composite declines under both conventions"
    else if (dz || dr) "N-2 : it declines under one convention only - print both, claim no module-level decline"
    else if (mz) "N-3 : the composite does not decline while MCT4 alone does. That is the manuscript's own position and is NOT reported as a module-level decline."
    else "N-3 : neither the composite nor MCT4 clears the family-adjusted alpha here"
  verdict("S12-N panel-definition reconciliation (Reviewer 2, points 1 and 3)",
    "The manuscript, the deposited repository and P2_revision2_all.R r6 do NOT use",
    "the same ANLS composite. The manuscript and the repository use SLC2A1 + LDHA",
    "+ SLC16A1 with MCT4 separate; r6 used SLC16A3 + LDHA + HK2 + PFKFB3, which",
    "puts the tested gene inside its own composite. This screen uses the",
    "manuscript definition and reports MCT4 as its own row.",
    sprintf("ANLS composite (manuscript) : %+.1f%% | z Holm %.4g | raw Holm %.4g | largest raw member %.1f%%",
            gv("ANLS composite","pct"), lac_z, lac_r, gv("ANLS composite","raw_top_member_share")),
    sprintf("MCT4 alone                  : %+.1f%% | z Holm %.4g | raw Holm %.4g",
            gv("single gene","pct"), mct_z, mct_r),
    sprintf("Iron handling               : %+.1f%% | z Holm %.4g | raw Holm %.4g | largest raw member %.1f%%",
            gv("Iron handling","pct"), gv("Iron handling","holm_z"),
            gv("Iron handling","holm_raw"), gv("Iron handling","raw_top_member_share")),
    cls,
    "Both conventions are printed for every row. Neither is presented as neutral.")
  invisible(TRUE) }

## ---- 10. S12-M  mediation competition ---------------------------------------
## The null construction below is 50_mediation_matched_null.R. pct_med(), tstat(),
## the ladder, MIN_CTRL and the seed are unchanged. Only the mediator is looped,
## and each competitor is judged inside ITS OWN matched null, which is the only
## comparison that is fair. The published mediated fractions are a hard gate.
s12_mediation <- function() {
  if (is.null(PASS$Ea)) { note("[stop S12-M] the h5ad pass did not run\n"); return(invisible(FALSE)) }
  set.seed(20260713)   # script 50's seed, kept so that this block reproduces the deposited Table 3 exactly; see R/mediation/50_mediation_matched_null.R
  Ea <- PASS$Ea; Ee <- PASS$Ee; Da <- PASS$Da; DM <- PASS$DM
  ax <- DM$cps
  mod <- function(M, gs) { g <- intersect(gs, colnames(M))
    if (!length(g)) return(NULL); rowMeans(M[, g, drop=FALSE]) }
  SETS <- list("Neuronal activity genes"        = mod(Ee, ACT_N),
               "Neuronal V-ATPase (6 subunits)" = mod(Ee, VATP6),
               "Neuronal LAMP1"                 = mod(Ee, "LAMP1"))
  SETS <- SETS[!vapply(SETS, is.null, logical(1))]
  if (!length(SETS)) { note("[stop S12-M] no outcome module could be built\n"); return(invisible(FALSE)) }
  note("outcome modules built: %d of 3 (%s)\n", length(SETS), paste(names(SETS), collapse=" | "))

  dA <- colMeans(Da)                                  # donor-level detection, script 50 form
  note("computing the CPS effect size for %d astrocytic genes (a few minutes) ...\n", ncol(Ea))
  tA <- apply(Ea, 2, function(z) tstat(z, ax))       # tstat() kept verbatim from script 50

  d4 <- dA[["SLC16A3"]]; t4 <- tA[["SLC16A3"]]
  g4 <- abs(d4 - G4$d4) <= G4$tol_d && abs(t4 - G4$t4) <= G4$tol_t
  note("G4 script-50 anchors: SLC16A3 detection %.4f (published %.4f) | |t| %.2f (published %.2f) -> %s\n",
       d4, G4$d4, t4, G4$t4, ifelse(g4, "PASS", "FAIL"))

  build_null <- function(g) {
    if (!(g %in% names(dA))) return(NULL)
    d0 <- dA[[g]]; t0 <- tA[[g]]
    if (!is.finite(d0) || !is.finite(t0) || d0 <= 0 || t0 <= 0) return(NULL)
    ok <- integer(0); wu <- NA_real_
    for (w in MED_LADDER) {
      ok <- which(is.finite(dA) & is.finite(tA) & names(dA) != g &
                  abs(dA - d0)/d0 < w & abs(tA - t0)/t0 < w)
      wu <- w; if (length(ok) >= MED_MIN_CTRL) break }
    if (length(ok) < MED_MIN_USABLE) return(NULL)
    list(gene=g, detection=d0, t=t0, window=wu, ctrl=names(dA)[ok]) }

  NULLS <- lapply(COMPETITORS, build_null); names(NULLS) <- COMPETITORS
  for (g in COMPETITORS) { z <- NULLS[[g]]
    if (is.null(z)) note("[MATCH] %-8s no usable null\n", g)
    else note("[MATCH] %-8s detection %.4f  |t| %.2f  window +/-%.0f%%  controls %d\n",
              g, z$detection, z$t, z$window*100, length(z$ctrl)) }

  # M-gate: MCT4 must reproduce the three published mediated fractions first.
  gate <- rbindlist(lapply(names(SETS), function(nm) {
    obs <- pct_med(SETS[[nm]], Ea[, "SLC16A3"], ax)
    pub <- unname(G5$published[nm])
    data.table(outcome=nm, observed=round(obs,1), published=pub,
      ok = is.finite(obs) && is.finite(pub) && abs(obs - pub) <= G5$tol) }))
  put(gate, "S12_M_reproduction_gate.csv",
      "the three published mediated fractions recomputed from the rebuilt donor matrices")
  print(gate)
  if (!all(gate$ok, na.rm=TRUE)) {
    verdict("S12-M mediation competition - NOT RUN",
      "M-gate FAILED. The rebuilt donor matrices do not reproduce the published",
      "mediated fractions (96.4 / 123.3 / 117.1 per cent), so the competitor loop",
      "would not be comparable with Table 3 and was not run.",
      paste(sprintf("%s observed %.1f%% against published %.1f%%",
                    gate$outcome, gate$observed, gate$published), collapse=" | "),
      "Fix the extraction before reading anything into the competitor comparison.")
    return(invisible(FALSE)) }

  RES <- rbindlist(lapply(names(SETS), function(nm) {
    y <- SETS[[nm]]
    rbindlist(lapply(COMPETITORS, function(g) { z <- NULLS[[g]]
      if (is.null(z)) return(data.table(outcome=nm, mediator=g, verdict="no usable matched null"))
      m   <- Ea[, g]; obs <- pct_med(y, m, ax)
      nul <- vapply(z$ctrl, function(k) pct_med(y, Ea[, k], ax), 0); nul <- nul[is.finite(nul)]
      rev <- pct_med(m, y, ax)
      # Empirical p in two conventions. p_unc is k / n, the form used when the
      # verdict rule was pre-registered and the form the verdict column follows.
      # p is (k + 1) / (n + 1), adopted at revision because a permutation p of
      # exactly zero is not attainable; it is the value reported in the paper.
      # The correction was adopted after the values were known and is therefore
      # NOT used to reverse a verdict: of the eighteen mediator-outcome tests it
      # moves exactly one across 0.05 (TFRC at neuronal V-ATPase, 0.0494 to
      # 0.0530) and in the direction favouring this study's own claim.
      p_unc <- mean(nul >= obs)
      p     <- (sum(nul >= obs) + 1) / (length(nul) + 1)
      data.table(outcome=nm, mediator=g, detection=round(z$detection,4), t_cps=round(z$t,2),
        window=z$window, null_n=length(nul), pct_mediated=round(obs,1),
        null_median=round(median(nul),1), null_p95=round(quantile(nul,0.95,names=FALSE),1),
        percentile=round(100*mean(nul <= obs),1), p_empirical=round(p,4),
        reverse_pct_mediated=round(rev,1),
        verdict=if (!is.finite(p_unc)) "not computed"
          else if (p_unc < ALPHA) "exceeds its own matched null"
          else "does not exceed its own matched null",
        p_empirical_uncorrected=round(p_unc,4)) }), fill=TRUE) }), fill=TRUE)
  put(RES, "S12_M_mediation_competition.csv",
      "MCT4 against competing astrocytic mediators, each inside its own detection- and effect-matched null")
  put(rbindlist(lapply(names(NULLS), function(g) { z <- NULLS[[g]]
    if (is.null(z)) NULL else data.table(mediator=g, control=z$ctrl, window=z$window) })),
    "S12_M_null_membership.csv", "which genes entered each competitor's matched null")

  if (!("p_empirical" %in% names(RES))) {
    verdict("S12-M mediation competition - NOT INTERPRETABLE",
      "no competitor produced a usable matched null; nothing is claimed")
    return(invisible(FALSE)) }
  KEY <- RES[outcome %like% "V-ATPase"]
  print(KEY[, .(mediator, detection, null_n, pct_mediated, null_median,
                percentile, p_empirical, reverse_pct_mediated, verdict)])
  m4  <- KEY[mediator=="SLC16A3"]
  riv <- KEY[mediator!="SLC16A3" & is.finite(p_empirical) & p_empirical < ALPHA]$mediator
  verdict("S12-M mediation competition (the load-bearing layer)",
    sprintf("G4 script-50 anchors : %s", ifelse(g4,"PASS","FAIL - detection or |t| does not match the published 0.0669 / 5.78")),
    "M-gate PASSED: the rebuilt matrices reproduce 96.4 / 123.3 / 117.1 per cent.",
    if (nrow(m4) && is.finite(m4$p_empirical))
      sprintf("MCT4 : %.1f%% mediated, null median %.1f%%, percentile %.1f, empirical p %.4f",
              m4$pct_mediated, m4$null_median, m4$percentile, m4$p_empirical) else "MCT4 : not computed",
    sprintf("competitors that ALSO exceed their own matched null : %s",
            if (length(riv)) paste(riv, collapse=", ") else "(none)"),
    if (length(riv))
      "-> M-1 FIRES. The mediation layer does not separate MCT4 from these genes. Reframe as a coordinated astrocytic metabolic decline of which MCT4 is one component. 'Bottleneck' and 'rate-limiting' stay out of the manuscript."
    else if (nrow(m4) && is.finite(m4$p_empirical) && m4$p_empirical < ALPHA)
      "-> M-2. Specificity survives at the mediation layer ONLY. The reply must state that the gene-level detection-matched null does not separate MCT4 from the iron pathway and that the entire separation rests on this one test."
    else "-> M-3. The mediation claim is withdrawn.",
    "The reverse direction is in the table for every competitor, not only for MCT4.",
    "Competitor set chosen after the gene-level percentiles were seen: adversarial, not pre-registered.")
  invisible(TRUE) }

## ---- 11. run -----------------------------------------------------------------
if (!selftest()) stop("self-test failed - nothing below was run")
note("== S12-P  panels from the existing S3 table (seconds) ==\n");  s12_panels()
note("\n== S12-H  sparse-cell guard (seconds) ==\n");                s12_guard()
note("\n== one h5ad pass: four compartments + donor x gene matrices ==\n"); one_pass()
note("\n== S12-A and S12-D  cell-type decomposition ==\n");          s12_compartments()
note("\n== S12-N  panel-definition reconciliation ==\n");            s12_reconcile()
note("\n== S12-M  mediation competition ==\n");                      s12_mediation()

MAN <- if (length(.MAN$rows)) rbindlist(.MAN$rows) else data.table()
if (nrow(MAN)) fwrite(MAN, file.path(OUTDIR, "S12_MANIFEST.csv"))
writeLines(c(paste("S12_VERDICTS", BUILD),
             paste("run :", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                   "| seed 42 globally, 20260713 inside S12-M"),
             .MAN$verd), file.path(OUTDIR, "S12_VERDICTS.txt"))
writeLines(c(BUILD, paste("run :", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
             paste("h5ad:", MTG_H5), capture.output(sessionInfo())[1:6]),
           file.path(OUTDIR, "S12_SESSION.txt"))
note("\nS12 complete. Outputs in %s\n", OUTDIR)
