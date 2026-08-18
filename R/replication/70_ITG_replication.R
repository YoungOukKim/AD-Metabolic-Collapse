# =============================================================================
# 70_ITG_replication.R
#   External replication of the astrocytic MCT4 decline in an independent cohort
#   (MGH-AbbVie AD Progression Atlas, Serrano-Pozo et al., Nat Neurosci 2024).
#
# WHY THIS EXISTS
#   P2 states that "confirmation in a second intermediate-vulnerability region is
#   needed". The inferior temporal gyrus (ITG) of this atlas is that region: it is
#   temporal association neocortex adjacent to the MTG, staged along the same tau
#   progression, with 31 donors. The data are astrocyte- and microglia-enriched
#   (neurons and oligodendrocytes were removed by FANS), so:
#
#     testable here : the astrocytic MCT4 decline, its detection-matched
#                     specificity, and the ambient control against microglia
#     NOT testable  : the MCT4 x neuronal V-ATPase coupling (no neurons)
#
# WHAT IS FIXED BEFORE THE DATA IS TOUCHED
#   Axis. The download metadata carries Pathology Stage (1-4), Braak, CERAD, Sex,
#   Donor ID and APOE. It does NOT carry the pTau/tau ratio used in the paper, so
#   the pre-registered axis is Pathology Stage: 4 ordered levels, 7/8/8/8 donors,
#   complete for all 31. Braak is not used as primary because three donors are
#   coded "unk" or "0".
#
#   Power. n = 31. With 80% power this design detects |rho| >= 0.48 on a
#   continuous axis and somewhat larger on a 4-level ordinal one. SEA-AD observed
#   rho = -0.502. A null therefore excludes a SEA-AD-sized effect but does NOT
#   exclude a smaller one, so the result is reported as an interval, never as a
#   bare p-value, and "not replicated" is only written if the interval excludes
#   the SEA-AD estimate.
#
#   APOE. 21 donors are 3/3, 7 are 3/4 or 4/4, 3 carry an e2 allele. Seven
#   carriers is below this project's MIN_CARRIERS of 10, so no APOE
#   stratification is attempted here and none will be reported.
#
# Usage
#   Paths are RELATIVE to the repository root. Place the atlas downloads in
#   data/raw/astMet/ or override with ASTMET_DIR. The raw matrices are not
#   redistributable; obtain them from https://ad-progression-atlas.partners.org
#
#   Rscript 70_ITG_replication.R selftest   logic only, no data
#   Rscript 70_ITG_replication.R            full run
# =============================================================================
suppressMessages({ library(data.table) })
options(stringsAsFactors = FALSE)
set.seed(42)
SELFTEST_ONLY <- identical(commandArgs(TRUE)[1], "selftest")

## ---- config -----------------------------------------------------------------
DIR    <- Sys.getenv("ASTMET_DIR",  unset = "data/raw/astMet")   # atlas downloads
OUTDIR <- Sys.getenv("ITG_OUT_DIR", unset = "output/tables/ITG_replication")
PFX    <- Sys.getenv("ASTMET_PREFIX", unset = "2026-07-19")
REGION <- Sys.getenv("ASTMET_REGION", unset = "ITG")   # ITG (pre-registered) or EC, PFC, V1, V2
CHUNK  <- 2e7          # MTX triplets read per block; lower if memory-limited
N_BG   <- 4000L        # background genes for the detection-matched null
N_BOOT <- 5000L

## ---- pre-registered thresholds and rules ------------------------------------
MIN_CELL      <- 10L    # cells per donor for that donor to enter the pseudobulk
MIN_DONOR     <- 20L    # donors required for any test
DET_MIN       <- 0.05   # detection floor, as in P2
## The source paper's own QC, quoted from its Methods. The download is the
## PRE-filter matrix; the published analyses use these thresholds, and the paper
## states the resulting counts for BA20 (= ITG), which makes them a technical gate.
QC_MIN_GENES  <- 2000L  # nuclei with more than this many genes detected
QC_MAX_UMI    <- 25000  # and fewer than this many UMIs after gene filtering
PAPER_PC_GENES<- 18283L # "18,283 protein-coding genes (noncoding genes were excluded)"
QC_GENE_MINCT <- 100    # genes with sum of counts above this
QC_GENE_FRAC  <- 0.30   # in at least this fraction of samples
PAPER_ITG_NUC <- 70984L # "n = 70,984 nuclei for BA20"
PAPER_ITG_GEN <- 11148L # "n = 11,148 for BA20"
QC_TOL        <- 0.10   # the reader must land within this fraction of both
BAND_LO       <- 0.75   # detection-matched band, multiples of MCT4's detection
BAND_HI       <- 1.25
ALPHA         <- 0.05
SEAAD_RHO     <- -0.502 # the SEA-AD MTG estimate this run is trying to replicate
MDE           <- -0.48  # smallest effect this n can detect at 80% power
NULL_PCTL     <- 5      # matched-null percentile a replication must reach
## G1 gate: the reader must reproduce a result the source paper reports.
## Temporal gene set 3 (proteostasis / heat-shock / metallothionein) peaks at
## pathology stage 3 and returns towards baseline at stage 4; the authors
## validated HSPB1 by RNAScope with exactly that shape. A composite is used
## rather than one gene so the gate does not hinge on a single measurement.
## CORRECTED. The first version of this gate was curated by hand from the paper's
## prose and contained MT2A, which the authors assign to temporal gene set 1
## (highest at the earliest stage, then falling) -- the opposite shape to the one
## being tested. It is removed, MT1G and MT1M are added, and the proteostasis
## members are taken in full as the paper lists them for set 3.
##
## Better still, the authors published the actual gene sets as Supplementary
## Data 5 of the paper. If that file is placed in ASTMET_DIR the gate uses it and
## no judgement of ours enters the gate at all. G1_GENES is only the fallback.
G1_SUPP <- Sys.getenv("G1_SUPPDATA5", unset = file.path(DIR, "41593_2024_1791_MOESM7_ESM.xlsx"))
G1_GENES <- c(
  ## set 3, proteostasis / heat-shock response
  "AHSA1","BAG3","CCT4","CRYAB","DNAJA1","DNAJB1","DNAJB6","HSP90AA1","HSP90AB1",
  "HSPA1A","HSPA1B","HSPA4","HSPA4L","HSPA8","HSPA9","HSPB1","HSPD1","HSPE1",
  "HSPH1","MIB1","OTUD7B","SQSTM1","ST13","UBB","UBC","UBE2B","UBE2E2",
  ## set 3, antioxidant defence and metallothioneins  (MT2A deliberately absent)
  "NFE2L2","NXN","PRDX1","SOD1","SOD2","MT1E","MT1F","MT1G","MT1M","MT1X","MT3",
  ## set 3, energy metabolism and ETC
  "ENO1","GAPDH","LDHB","PGK1","ATP5F1E","COX6A1","COX6C","COX7C",
  "NDUFA4","NDUFB2","NDUFB4","NDUFC1","UQCRB")
G1_WRONG <- c("HSPB1","HSPA1A","HSPA1B","CRYAB","DNAJB1","HSPH1","HSP90AA1",
              "MT1X","MT2A","MT3","MT1E","MT1F")   # the first version, reported for the record
G1_MIN_RATIO <- 1.05    # stage-3 composite must exceed stage-1/2 by this factor
                        # AND must fall from stage 3 to stage 4

TARGET   <- "SLC16A3"
ANCHOR2  <- c("SLC1A2","SLC1A3","AQP4","GFAP")   # identity genes, reported alongside

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
die <- function(...) { cat(sprintf("\n*** ABORT: %s\n", sprintf(...))); quit(status = 1) }

## =============================================================================
## SECTION A -- statistics
## =============================================================================

## verbatim from R/01_mtg_b55a_robust.R -- do not edit
pcor_sp <- function(x, y, z) {
  rz <- rank(z)
  ex <- residuals(lm(rank(x) ~ rz))
  ey <- residuals(lm(rank(y) ~ rz))
  cor(ex, ey)
}
## Fisher interval for a Spearman/partial estimate on n donors, k covariates.
rho_ci <- function(r, n, k = 1L, conf = 0.95) {
  if (!is.finite(r) || n - k - 3L <= 0L) return(c(NA_real_, NA_real_))
  se <- 1 / sqrt(n - k - 3); z <- atanh(r); cr <- qnorm(1 - (1 - conf) / 2)
  tanh(c(z - cr * se, z + cr * se))
}
rho_p <- function(r, n, k = 1L) {
  df <- n - k - 2L
  if (!is.finite(r) || df <= 0L) return(NA_real_)
  2 * pt(-abs(r * sqrt(df / (1 - r^2))), df = df)
}
boot_rho <- function(x, y, z, B = N_BOOT) {
  ok <- complete.cases(x, y, z); x <- x[ok]; y <- y[ok]; z <- z[ok]; n <- length(x)
  v <- vapply(seq_len(B), function(b) { i <- sample.int(n, n, replace = TRUE)
    suppressWarnings(tryCatch(pcor_sp(x[i], y[i], z[i]), error = function(e) NA_real_)) }, numeric(1))
  unname(quantile(v, c(0.025, 0.975), na.rm = TRUE))
}

## =============================================================================
## SECTION B -- self-test. Runs every time.
## =============================================================================
selftest <- function() {
  ok <- TRUE
  chk <- function(c, m) { cat(sprintf("  [%s] %s\n", ifelse(c, "PASS", "FAIL"), m)); ok <<- ok && isTRUE(c) }
  cat("=== SELF-TEST ===\n"); set.seed(42)   # project seed; self-test assertions are structural or distribution-based and do not depend on the draw

  ## B1 partial Spearman removes a purely shared covariate
  n <- 31
  pr <- replicate(1000, { z <- rnorm(n); x <- z + rnorm(n, 0, .2); y <- z + rnorm(n, 0, .2)
                         c(raw = cor(x, y), partial = pcor_sp(x, y, z)) })
  chk(median(abs(pr["partial", ])) < 0.25 && median(pr["raw", ]) > 0.8,
      sprintf("shared covariate removed: median |partial| %.2f vs raw %.2f (300 draws)",
              median(abs(pr["partial", ])), median(pr["raw", ])))
  ## The entire verdict rests on the confidence interval, so the interval must be
  ## calibrated FOR THIS ESTIMATOR at this n. Under the null the Fisher-z of a
  ## partial correlation has sd = 1/sqrt(n-k-3); check the simulated width against
  ## that, rather than against a threshold chosen to make the test pass.
  zs <- atanh(pr["partial", ]); theo <- 1 / sqrt(n - 1 - 3)
  chk(abs(sd(zs) - theo) < 0.05,
      sprintf("CI is calibrated for pcor_sp at n=%d: simulated sd(z) %.3f vs theory %.3f", n, sd(zs), theo))

  ## B2 the interval must be wide enough to make an underpowered null visible
  ci <- rho_ci(-0.20, 31, 1)
  chk(ci[2] > 0, sprintf("a weak estimate at n=31 gives an interval crossing zero [%.2f, %.2f]", ci[1], ci[2]))
  ci2 <- rho_ci(SEAAD_RHO, 31, 1)
  chk(ci2[2] < 0, sprintf("a SEA-AD-sized estimate does not cross zero [%.2f, %.2f]", ci2[1], ci2[2]))

  ## B3 MTX triplet accumulation must equal a dense reference
  ng <- 40L; nc <- 60L
  M <- matrix(rpois(ng * nc, 0.6), ng, nc)
  tri <- which(M > 0, arr.ind = TRUE); val <- M[M > 0]
  don <- rep(1:5, length.out = nc)
  ref <- t(sapply(1:5, function(d) rowSums(M[, don == d, drop = FALSE])))
  acc <- matrix(0, 5, ng)
  lin <- don[tri[, 2]] + (tri[, 1] - 1L) * 5L
  s <- rowsum(val, lin, reorder = FALSE); i <- as.integer(rownames(s))
  acc[i] <- acc[i] + s[, 1]
  chk(max(abs(acc - ref)) < 1e-9, sprintf("triplet accumulation equals the dense sum (max diff %.1e)", max(abs(acc - ref))))

  ## B4 blocked reading of the same triplets must give the same answer
  ord <- sample(nrow(tri)); tri <- tri[ord, ]; val <- val[ord]
  acc2 <- matrix(0, 5, ng)
  for (a in seq(1, nrow(tri), by = 137)) {
    b <- min(a + 136, nrow(tri))
    l <- don[tri[a:b, 2]] + (tri[a:b, 1] - 1L) * 5L
    s <- rowsum(val[a:b], l, reorder = FALSE); i <- as.integer(rownames(s))
    acc2[i] <- acc2[i] + s[, 1] }
  chk(max(abs(acc2 - ref)) < 1e-9, "blocked reading equals single-pass reading")

  ## B5 the G1 gate must accept the published shape and reject a flat one
  shape_ok <- function(v) (v[3] / mean(v[1:2]) >= G1_MIN_RATIO) && (v[4] < v[3])
  chk(shape_ok(c(1.0, 1.0, 1.4, 1.1)), "G1 accepts a stage-3 peak that falls at stage 4")
  chk(!shape_ok(c(1.0, 1.0, 1.0, 1.0)), "G1 rejects a flat trajectory")
  chk(!shape_ok(c(1.0, 1.1, 1.3, 1.5)), "G1 rejects a monotonic rise")

  ## B6 the verdict rule must be decidable in all three regimes
  verdict <- function(r, ci) {
    if (!is.finite(r)) return("NOT TESTABLE")
    if (r <= MDE && ci[2] < 0) return("REPLICATED")
    if (r < 0 && ci[2] >= 0)   return("DIRECTION CONSISTENT, UNDERPOWERED")
    if (ci[1] > MDE)           return("NOT REPLICATED")
    "INCONCLUSIVE" }
  chk(verdict(-0.55, rho_ci(-0.55, 31)) == "REPLICATED", "verdict: strong negative -> REPLICATED")
  chk(verdict(-0.20, rho_ci(-0.20, 31)) == "DIRECTION CONSISTENT, UNDERPOWERED",
      "verdict: weak negative -> underpowered, not a failure")
  chk(verdict(+0.10, rho_ci(+0.10, 31)) == "NOT REPLICATED", "verdict: positive -> NOT REPLICATED")

  cat(sprintf("=== SELF-TEST %s ===\n\n", ifelse(ok, "PASSED", "FAILED")))
  if (!ok) die("self-test failed; nothing was read")
  invisible(TRUE)
}
selftest()
if (SELFTEST_ONLY) { cat("selftest mode: stopping before touching data.\n"); quit(status = 0) }

## =============================================================================
## SECTION C -- read one cell type, streaming the MTX
## =============================================================================
fp <- function(...) file.path(DIR, paste0(PFX, "_", paste0(..., collapse = "")))

read_meta <- function(ct) {
  ## the download ships microglia metadata as CSV and astrocyte metadata as XLSX
  cands <- c(fp(ct, "_", REGION, "_metadata.csv"), fp(ct, "_", REGION, "_metadata.xlsx"),
             fp(ct, "_metadata.csv"), fp(ct, "_metadata.xlsx"))
  f <- cands[file.exists(cands)][1]
  if (is.na(f)) stop(sprintf("no metadata file for %s. Looked for:\n    %s", ct,
                             paste(cands, collapse = "\n    ")))
  if (grepl("\\.xlsx$", f)) {
    if (!requireNamespace("readxl", quietly = TRUE))
      die("%s is XLSX and package 'readxl' is not installed.\n  Either install.packages('readxl') or open it in Excel and save as\n  %s", f, sub("\\.xlsx$", ".csv", f))
    d <- as.data.table(readxl::read_excel(f))
  } else d <- fread(f)
  setnames(d, 1, "barcode")
  cat(sprintf("[meta  ] %s: %s rows | %s\n", basename(f), format(nrow(d), big.mark = ","),
              paste(setdiff(names(d), "barcode"), collapse = ", ")))
  d
}

load_ct <- function(ct) {
  cells <- fread(fp(ct, "_", REGION, "_cell_annotation.txt"), header = FALSE)[[1]]
  genes <- fread(fp(ct, "_", REGION, "_row_annotation.txt"),  header = FALSE)[[1]]
  mtx   <- fp(ct, "_", REGION, "_matrix.mtx")
  if (!file.exists(mtx)) die("matrix not found: %s", mtx)
  meta  <- read_meta(ct)

  ## ---- MTX header: skip comments, read dims -------------------------------
  con <- file(mtx, "r"); hdr <- 0L
  repeat { l <- readLines(con, 1L); hdr <- hdr + 1L; if (!startsWith(l, "%")) break }
  close(con)
  dims <- as.numeric(strsplit(trimws(l), "\\s+")[[1]])
  nr <- dims[1]; nc <- dims[2]; nnz <- dims[3]
  cat(sprintf("[mtx   ] %s: %s x %s, %s non-zeros (header %d lines)\n", basename(mtx),
              format(nr, big.mark=","), format(nc, big.mark=","), format(nnz, big.mark=","), hdr))

  ## ---- orientation: rows must be genes, columns cells ---------------------
  if (nr == length(genes) && nc == length(cells))      transpose <- FALSE
  else if (nr == length(cells) && nc == length(genes)) transpose <- TRUE
  else die("matrix %g x %g matches neither %d genes x %d cells nor its transpose",
           nr, nc, length(genes), length(cells))
  cat(sprintf("[mtx   ] orientation: %s\n", if (transpose) "cells x genes -> transposed" else "genes x cells"))

  ## ---- donor key ----------------------------------------------------------
  ## Astrocytes and microglia were sorted from the SAME sequencing libraries, so
  ## the library ID in a barcode suffix identifies the donor for either cell type.
  ## If this cell type has no metadata sheet for this region, borrow the companion
  ## cell type's sheet: it supplies donor, stage, Braak and APOE, all of which are
  ## constant within a library. Verified: 31/31 ITG libraries shared, 100% of
  ## astrocyte cells assigned.
  lib_of0 <- function(b) sub("^[ACGTN]+-1_", "", b)
  ## The library ID lives INSIDE the barcode on both sheets, not in a column of
  ## its own, so materialise it before testing coverage or joining. Omitting this
  ## made the companion sheet score 0 and the swap silently never happened.
  add_lib <- function(md) { if (is.null(md) || !nrow(md)) return(md)
    md <- as.data.table(md); md[, .lib := lib_of0(as.character(barcode))]; md[] }
  meta <- add_lib(meta)
  meta_covers <- function(md) {
    if (is.null(md) || !nrow(md)) return(0)
    ul <- unique(lib_of0(cells))
    max(c(mean(cells %in% as.character(md$barcode)),
          vapply(names(md), function(cn) mean(ul %in% unique(as.character(md[[cn]]))), numeric(1))),
        na.rm = TRUE) }
  if (meta_covers(meta) < 0.9) {
    other <- setdiff(c("Astrocytes", "Microglia"), ct)
    alt <- add_lib(tryCatch(read_meta(other), error = function(e) NULL))
    if (!is.null(alt) && meta_covers(alt) >= 0.9) {
      cat(sprintf("[meta  ] this sheet covers %.0f%% of %s %s; the %s sheet covers %.0f%% -- using it for the\n",
                  100*meta_covers(meta), ct, REGION, other, 100*meta_covers(alt)))
      cat("         library-to-donor mapping (both cell types come from the same libraries).\n")
      meta <- alt
    } }

  ## The two downloads use different schemas: the microglia CSV is the curated
  ## per-cell table (Donor ID / Pathology Stage), the astrocyte XLSX is the raw
  ## Seurat metadata (Donor.ID.y / Path..Group. / Ptau.Total.Tau..A.U..). Detect
  ## rather than assume, and prefer the continuous pTau/tau axis when it is there.
  pick1 <- function(pats) { for (p in pats) { h <- grep(p, names(meta), value = TRUE, ignore.case = TRUE)
                                              if (length(h)) return(h[1]) }; NA_character_ }
  dcol <- pick1(c("^Donor.?ID", "^Donor"))
  scol <- pick1(c("^Pathology.?Stage", "^Path\\.{2}Group", "^Path.?Group", "^Path\\."))
  pcol <- pick1(c("^Ptau.*Total.*Tau", "pTau.*tau", "^Ptau"))
  rcol <- pick1(c("^Unified_region$", "^Brain.?region$", "^Region$"))
  cat(sprintf("[cols  ] donor=%s | stage=%s | pTau/tau=%s | region=%s\n",
              dcol, scol, ifelse(is.na(pcol), "-", pcol), ifelse(is.na(rcol), "-", rcol)))
  if (is.na(dcol)) die("no donor column found. Columns present:\n    %s",
                       paste(names(meta), collapse = ", "))
  if (is.na(scol) && is.na(pcol))
    die("neither a pathology-stage nor a pTau/tau column was found. Columns present:\n    %s",
        paste(names(meta), collapse = ", "))
  ## keep only this region when the table spans several
  if (!is.na(rcol)) {
    reg <- toupper(trimws(as.character(meta[[rcol]])))
    cat(sprintf("[region] %s values: %s\n", rcol,
                paste(sprintf("%s(%d)", names(table(reg)), as.integer(table(reg))), collapse = " ")))
    want <- switch(REGION, ITG = "ITG|BA20|INFERIOR", EC = "^EC$|ENTORHIN",
                   PFC = "PFC|BA46|PREFRONT", V1 = "^V1$|BA17", V2 = "^V2$|BA18|BA19", REGION)
    if (!any(grepl(want, reg)) && meta_covers(meta) < 0.9)
      die(paste("this metadata sheet covers %s, not %s. The atlas download page serves the",
                "Metadata file for whichever region is selected, so re-download it with %s active",
                "under the %s tab. Alternatively set ASTMET_REGION to the region you have."),
          paste(names(table(reg)), collapse = "/"), REGION, REGION, ct)
    sel <- grepl(want, reg)
    if (any(sel) && !all(sel)) { before <- nrow(meta); meta <- meta[sel]
      cat(sprintf("[region] filtered to ITG: %s of %s rows\n",
                  format(nrow(meta), big.mark=","), format(before, big.mark=","))) } }
  meta[[dcol]] <- as.character(meta[[dcol]])

  ## ---- join ---------------------------------------------------------------
  ## Only donor identity and the progression axis are needed, and both are
  ## constant within a sequencing library. So if per-cell barcodes do not match
  ## between the matrix and the metadata sheet, fall back to joining on the
  ## library ID carried in the barcode suffix (AAACCC...-1_6289-MW-0005).
  lib_of <- function(b) sub("^[ACGTN]+-1_", "", b)
  m <- meta[match(cells, meta$barcode)]
  hit <- sum(!is.na(m[[dcol]]))
  how <- "barcode"
  if (hit < 0.5 * length(cells)) {
    clib <- lib_of(cells); ulib <- unique(clib)
    score <- vapply(names(meta), function(cn) {
      v <- unique(as.character(meta[[cn]])); mean(ulib %in% v) }, numeric(1))
    best <- names(score)[which.max(score)]
    cat(sprintf("[join  ] barcode match %.1f%%; best library-ID column is '%s' covering %.1f%% of libraries\n",
                100*hit/length(cells), best, 100*max(score)))
    if (max(score) >= 0.9) {
      key <- as.character(meta[[best]])
      lut <- meta[!duplicated(key)]; rownames(lut) <- NULL
      m <- lut[match(clib, as.character(lut[[best]]))]
      hit <- sum(!is.na(m[[dcol]])); how <- sprintf("library ID via '%s'", best)
    } else {
      cat("\n  first 3 matrix barcodes : ", paste(head(cells, 3), collapse = " | "), "\n")
      cat("  their library suffixes  : ", paste(head(clib, 3), collapse = " | "), "\n")
      cat("  first 3 metadata keys   : ", paste(head(as.character(meta$barcode), 3), collapse = " | "), "\n\n")
      die(paste("cannot join %s ITG to this metadata sheet: %.1f%% by barcode, %.1f%% by library ID.",
                "This sheet has %s rows against %s ITG cells. On the atlas download page, select the",
                "Astrocytes tab and take the metadata that accompanies ITG (the microglia equivalent is",
                "%s_Microglia_ITG_metadata.csv), then rerun."),
          ct, 100*hit/length(cells), 100*max(score),
          format(nrow(meta), big.mark=","), format(length(cells), big.mark=","), PFX) } }
  cat(sprintf("[join  ] %s of %s cells matched on %s (%.1f%%)\n",
              format(hit, big.mark=","), format(length(cells), big.mark=","), how, 100*hit/length(cells)))
  if (hit < 0.5 * length(cells)) die("join still below 50%% after fallback")
  don_f  <- factor(m[[dcol]]); donors <- levels(don_f); nD <- nlevels(don_f)
  dcode  <- as.integer(don_f)
  one_per_donor <- function(col) vapply(donors, function(d) {
      v <- suppressWarnings(as.numeric(m[[col]][which(m[[dcol]] == d)]))
      v <- v[is.finite(v)]; if (!length(v)) NA_real_ else median(v) }, numeric(1))
  axis_name <- NA_character_; axis_d <- NULL
  if (!is.na(pcol)) { a <- one_per_donor(pcol)
    if (sum(is.finite(a)) >= MIN_DONOR && length(unique(a[is.finite(a)])) > 4) {
      axis_d <- a; axis_name <- pcol } }
  if (is.null(axis_d) && !is.na(scol)) { a <- one_per_donor(scol)
    if (sum(is.finite(a)) >= MIN_DONOR) { axis_d <- a; axis_name <- scol } }
  if (is.null(axis_d)) die("no usable progression axis at donor level")
  stage_d <- if (!is.na(scol)) one_per_donor(scol) else rep(NA_real_, nD)
  nG <- length(genes)
  cat(sprintf("[axis  ] %s: %s | %d donors with a value | range %.3f - %.3f | %d distinct\n",
              ct, axis_name, sum(is.finite(axis_d)), min(axis_d, na.rm=TRUE), max(axis_d, na.rm=TRUE),
              length(unique(axis_d[is.finite(axis_d)]))))
  if (all(is.finite(stage_d))) cat(sprintf("[units ] stage counts %s\n",
      paste(sprintf("%d:%d", 1:4, tabulate(round(stage_d), 4)), collapse = " ")))

  ## ---- stream the triplets ------------------------------------------------
  CELLG <- integer(length(cells)); CELLU <- numeric(length(cells))
  GSUM  <- numeric(nG); GSAMP <- matrix(0, nD, nG)
  SUM <- matrix(0, nD, nG); NZ <- matrix(0, nD, nG); TOT <- numeric(nD); CNT <- integer(nD)
  seen <- logical(length(cells))
  done <- 0; t0 <- Sys.time()
  ## One connection, read forward. fread(skip=) would rescan from byte zero on
  ## every chunk, which is quadratic across 334 million triplets.
  con <- file(mtx, "r"); invisible(readLines(con, hdr))
  on.exit(try(close(con), silent = TRUE), add = TRUE)
  while (done < nnz) {
    k <- min(CHUNK, nnz - done)
    raw <- scan(con, what = list(integer(), integer(), double()), nmax = k,
                quiet = TRUE, multi.line = FALSE)
    if (!length(raw[[1]])) break
    tr <- list(i = raw[[1]], j = raw[[2]], x = raw[[3]]); k <- length(tr$i)
    gi <- if (transpose) tr$j else tr$i
    ci <- if (transpose) tr$i else tr$j
    d  <- dcode[ci]
    keep <- !is.na(d)
    if (any(keep)) {
      gi <- gi[keep]; ci <- ci[keep]; d <- d[keep]; xv <- tr$x[keep]
      cg <- rowsum(rep(1, length(ci)), ci, reorder = FALSE); ic <- as.integer(rownames(cg))
      CELLG[ic] <- CELLG[ic] + cg[, 1]
      c2 <- rowsum(xv, ci, reorder = FALSE); iu <- as.integer(rownames(c2)); CELLU[iu] <- CELLU[iu] + c2[, 1]
      gs <- rowsum(xv, gi, reorder = FALSE); ig <- as.integer(rownames(gs)); GSUM[ig] <- GSUM[ig] + gs[, 1]
      lin <- d + (gi - 1L) * nD
      ga <- rowsum(xv, lin, reorder = FALSE); ia <- as.integer(rownames(ga)); GSAMP[ia] <- GSAMP[ia] + ga[, 1] }
    done <- done + k
    cat(sprintf("         %s / %s triplets  %.0fs\r", format(done, big.mark=","),
                format(nnz, big.mark=","), as.numeric(difftime(Sys.time(), t0, units="secs"))))
    rm(raw, tr, gi, ci, d, keep); gc(FALSE)
  }
  close(con)
  cat(sprintf("\n[read  ] %s pass 1 (%.1f min)\n", ct, as.numeric(difftime(Sys.time(), t0, units="mins"))))
  ## Step 1 of the paper's filter, which the first version of this script omitted:
  ## "18,283 protein-coding genes (that is, noncoding genes were excluded)". The
  ## download ships symbols only, with no biotype column, so use an annotation
  ## package when one is installed and fall back to a symbol heuristic otherwise.
  ## Either way the surviving count is reported against the paper's 18,283.
  pc_by_pkg <- function(sym) {
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) return(NULL)
    ok <- tryCatch(suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
            keys = sym, keytype = "SYMBOL", column = "ENTREZID")), error = function(e) NULL)
    if (is.null(ok)) NULL else !is.na(ok) }
  NONCODING <- paste(c("^(AC|AL|AP|AD|BX|CR|CU|FP|FO|Z)[0-9]{5,6}\\.[0-9]+$",
                       "^LINC[0-9]+$", "-AS[0-9]*$", "-DT$", "-OT[0-9]*$", "-IT[0-9]*$",
                       "^MIR[0-9]", "^MIRLET", "^SNOR", "^SCARNA", "^RNU[0-9]", "^RNA5S",
                       "^RNA5-8S", "^RN7S", "^RNY[0-9]", "^VTRNA", "^TRNA", "^LNCSRLR",
                       "^Y_RNA$", "^Metazoa_SRP$", "^[0-9]+$"), collapse = "|")
  pc <- pc_by_pkg(genes)
  how_pc <- "org.Hs.eg.db"
  if (is.null(pc)) { pc <- !grepl(NONCODING, genes); how_pc <- "symbol heuristic" }
  cat(sprintf("[pc    ] protein-coding filter (%s): %s of %s genes kept (paper states %s)\n",
              how_pc, format(sum(pc), big.mark=","), format(nG, big.mark=","),
              format(PAPER_PC_GENES, big.mark=",")))
  keep_gene <- pc & GSUM > QC_GENE_MINCT & colMeans(GSAMP > QC_GENE_MINCT) >= QC_GENE_FRAC
  ## recount detected genes per nucleus over the surviving gene set
  CELLG2 <- integer(length(cells)); done <- 0
  con <- file(mtx, "r"); invisible(readLines(con, hdr))
  while (done < nnz) {
    raw <- scan(con, what = list(integer(), integer(), double()), nmax = min(CHUNK, nnz - done), quiet = TRUE)
    if (!length(raw[[1]])) break
    g2 <- if (transpose) raw[[2]] else raw[[1]]; c2i <- if (transpose) raw[[1]] else raw[[2]]
    done <- done + length(raw[[1]])
    ok2 <- keep_gene[g2]
    if (any(ok2)) { cg2 <- rowsum(rep(1, sum(ok2)), c2i[ok2], reorder = FALSE)
                    i2 <- as.integer(rownames(cg2)); CELLG2[i2] <- CELLG2[i2] + cg2[, 1] }
    rm(raw, g2, c2i, ok2); gc(FALSE) }
  close(con)

  ## the >2,000-genes threshold is applied after gene filtering in the paper, so
  ## recount detected genes per nucleus over the surviving genes only
  keep_cell <- CELLG2 > QC_MIN_GENES & CELLU < QC_MAX_UMI & !is.na(dcode)
  cat(sprintf("[qc    ] %s: nuclei %s -> %s | genes %s -> %s\n", ct,
              format(sum(!is.na(dcode)), big.mark=","), format(sum(keep_cell), big.mark=","),
              format(nG, big.mark=","), format(sum(keep_gene), big.mark=",")))
  gmap2 <- rep(NA_integer_, nG); gmap2[which(keep_gene)] <- seq_len(sum(keep_gene))
  nG2 <- sum(keep_gene); genes2 <- genes[keep_gene]
  SUM <- matrix(0, nD, nG2); NZ <- matrix(0, nD, nG2); TOT <- numeric(nD); CNT <- integer(nD)
  seen <- logical(length(cells)); done <- 0; t1 <- Sys.time()
  con <- file(mtx, "r"); invisible(readLines(con, hdr))
  while (done < nnz) {
    raw <- scan(con, what = list(integer(), integer(), double()), nmax = min(CHUNK, nnz - done), quiet = TRUE)
    if (!length(raw[[1]])) break
    gi <- if (transpose) raw[[2]] else raw[[1]]; ci <- if (transpose) raw[[1]] else raw[[2]]
    xv <- raw[[3]]; done <- done + length(raw[[1]])
    ok <- keep_cell[ci] & !is.na(gmap2[gi])
    if (any(ok)) {
      ci <- ci[ok]; gp <- gmap2[gi[ok]]; xv <- xv[ok]; d <- dcode[ci]
      lin <- d + (gp - 1L) * nD
      s <- rowsum(xv, lin, reorder = FALSE); ii <- as.integer(rownames(s)); SUM[ii] <- SUM[ii] + s[, 1]
      z <- rowsum(rep(1, length(lin)), lin, reorder = FALSE); jj <- as.integer(rownames(z)); NZ[jj] <- NZ[jj] + z[, 1]
      td <- rowsum(xv, d, reorder = FALSE); tt <- as.integer(rownames(td)); TOT[tt] <- TOT[tt] + td[, 1]
      u <- !seen[ci]
      if (any(u)) { cu <- ci[u]; cu <- cu[!duplicated(cu)]; seen[cu] <- TRUE
        cc <- rowsum(rep(1, length(cu)), dcode[cu], reorder = FALSE)
        kk <- as.integer(rownames(cc)); CNT[kk] <- CNT[kk] + cc[, 1] } }
    cat(sprintf("         pass 2: %s / %s  %.0fs\r", format(done, big.mark=","),
                format(nnz, big.mark=","), as.numeric(difftime(Sys.time(), t1, units="secs"))))
    rm(raw, gi, ci, xv); gc(FALSE) }
  close(con)
  cat(sprintf("\n[read  ] %s pass 2 (%.1f min)\n", ct, as.numeric(difftime(Sys.time(), t1, units="mins"))))
  genes <- genes2; nG <- nG2
  qc <- list(n_nuclei = sum(keep_cell), n_genes = nG2, n_pc = sum(pc), how_pc = how_pc,
             mean_genes_per_nucleus = mean(CELLG2[keep_cell]))
  keepd <- CNT >= MIN_CELL
  if (sum(keepd) < MIN_DONOR) die("only %d donors reached %d cells in %s", sum(keepd), MIN_CELL, ct)
  list(genes = genes, qc = qc, donors = donors[keepd], stage = axis_d[keepd],
       axis_name = axis_name, stage_ord = stage_d[keepd],
       expr = SUM[keepd, , drop = FALSE] / CNT[keepd],       # per-cell mean, per donor
       det  = NZ[keepd, , drop = FALSE] / CNT[keepd],        # detection rate, per donor
       gmean = (TOT[keepd] / CNT[keepd]) / nG,               # global factor, as in P2
       ncell = CNT[keepd])
}

AST <- load_ct("Astrocytes")
MIC <- tryCatch(load_ct("Microglia"), error = function(e) { cat("[note  ] microglia not loaded: ", conditionMessage(e), "\n"); NULL })

## ---- is the matrix counts or already normalised? ----------------------------
is_counts <- all(abs(AST$expr[AST$expr > 0][1:min(1e5, sum(AST$expr > 0))] -
                     round(AST$expr[AST$expr > 0][1:min(1e5, sum(AST$expr > 0))])) < 1e-8)
cat(sprintf("\n[X     ] astrocyte matrix looks like %s\n",
            if (is_counts) "raw counts (per-donor means are non-integer, expected)" else "normalised values"))

## =============================================================================
## G1 GATE -- reproduce a result the source paper reports
## =============================================================================
cat("\n=== G1 GATE (technical): does the reader reproduce the paper's own QC counts? ===\n")
cat(sprintf("  paper, BA20 (= ITG): %s nuclei, %s genes, ~2,500 genes per nucleus\n",
            format(PAPER_ITG_NUC, big.mark=","), format(PAPER_ITG_GEN, big.mark=",")))
cat(sprintf("  this reader        : %s nuclei, %s genes, %.0f genes per nucleus\n",
            format(AST$qc$n_nuclei, big.mark=","), format(AST$qc$n_genes, big.mark=","),
            AST$qc$mean_genes_per_nucleus))
g1tech <- isTRUE(REGION != "ITG" ||
                 (abs(AST$qc$n_nuclei/PAPER_ITG_NUC - 1) <= QC_TOL &&
                  abs(AST$qc$n_genes /PAPER_ITG_GEN - 1) <= QC_TOL))
cat(sprintf("  deviation: nuclei %+.1f%% | genes %+.1f%% | tolerance %.0f%%\n  G1-technical >>> %s\n",
            100*(AST$qc$n_nuclei/PAPER_ITG_NUC - 1), 100*(AST$qc$n_genes/PAPER_ITG_GEN - 1), 100*QC_TOL,
            ifelse(g1tech, "PASS", "FAIL -- the extraction does not match the published cohort")))
fwrite(data.table(metric=c("protein_coding_genes","genes_after_count_filter","nuclei","genes_per_nucleus"),
                  paper=c(PAPER_PC_GENES,PAPER_ITG_GEN,PAPER_ITG_NUC,2500),
                  reader=c(AST$qc$n_pc,AST$qc$n_genes,AST$qc$n_nuclei,round(AST$qc$mean_genes_per_nucleus)),
                  method=c(AST$qc$how_pc,"","",""), pass=g1tech), file.path(OUTDIR,"G1_technical.csv"))

cat("\n=== G1b (reported, NOT gating): temporal set-3 stage-3 peak ===\n")
cat("  The paper derives its temporal gene sets ACROSS all five regions with region as a\n")
cat("  latent variable, and states the response may be attenuated or already exhausted in\n")
cat("  severely affected areas; ITG is the second most affected. A single-region test of a\n")
cat("  cross-region construct is therefore reported, not used as a gate.\n")
so <- round(AST$stage_ord)
if (!all(is.finite(so))) die("the G1 gate needs the 4-level pathology stage, which is absent")

## the authors' own list, if the supplementary file is present
supp_set3 <- NULL
if (file.exists(G1_SUPP) && requireNamespace("readxl", quietly = TRUE)) {
  sd5 <- tryCatch(as.data.table(readxl::read_excel(G1_SUPP)), error = function(e) NULL)
  if (!is.null(sd5)) {
    gcol <- names(sd5)[which(vapply(sd5, function(v) mean(as.character(v) %in% AST$genes), numeric(1)) > 0.5)][1]
    scol2 <- grep("set|cluster|group", names(sd5), value = TRUE, ignore.case = TRUE)[1]
    if (!is.na(gcol) && !is.na(scol2)) {
      hit <- grepl("3", as.character(sd5[[scol2]]))
      if (sum(hit) >= 20) { supp_set3 <- unique(as.character(sd5[[gcol]][hit]))
        cat(sprintf("  using the authors' published gene set: %d genes from %s\n",
                    length(supp_set3), basename(G1_SUPP))) } } } }
if (is.null(supp_set3))
  cat(sprintf("  Supplementary Data 5 not usable at %s; falling back to the list transcribed\n  from the paper's Methods. Place that file there to remove our judgement from the gate.\n",
              basename(G1_SUPP)))

run_gate <- function(gset, label) {
  gi <- match(intersect(gset, AST$genes), AST$genes)
  if (length(gi) < 6) return(data.table(gate = label, n_genes = length(gi), peak = NA, fall = NA, pass = NA))
  comp <- rowMeans(scale(AST$expr[, gi, drop = FALSE]))
  sm <- vapply(1:4, function(s) mean(comp[so == s], na.rm = TRUE), numeric(1))
  r <- sm - min(sm) + 1
  pk <- r[3] / mean(r[1:2]) >= G1_MIN_RATIO; fl <- r[4] < r[3]
  cat(sprintf("  %-22s %d genes | stage means %s | peak %s | falls %s\n", label, length(gi),
              paste(sprintf("%+.3f", sm), collapse = " "), ifelse(pk,"yes","NO"), ifelse(fl,"yes","NO")))
  data.table(gate = label, n_genes = length(gi), s1 = round(sm[1],4), s2 = round(sm[2],4),
             s3 = round(sm[3],4), s4 = round(sm[4],4), peak = pk, fall = fl, pass = isTRUE(pk && fl))
}
GT <- rbindlist(list(
  if (!is.null(supp_set3)) run_gate(supp_set3, "authors' Supp Data 5") else NULL,
  run_gate(G1_GENES, "paper Methods, corrected"),
  run_gate(G1_WRONG, "first version (MT2A in)")), fill = TRUE)
print(GT, row.names = FALSE); fwrite(GT, file.path(OUTDIR, "G1_gate.csv"))
## the primary gate is the authors' list when available, otherwise the corrected one
g1bio <- isTRUE(GT$pass[1]); g1ok <- g1tech
cat(sprintf("  G1b >>> %s  (list: %s)\n",
            ifelse(g1bio, "stage-3 peak present in ITG", "stage-3 peak NOT present in ITG alone"), GT$gate[1]))

## =============================================================================
## R1 -- the primary test
## =============================================================================
cat(sprintf("\n=== R1. ASTROCYTIC MCT4 vs %s ===\n", AST$axis_name))
ti <- match(TARGET, AST$genes); if (is.na(ti)) die("%s absent from the row annotation", TARGET)
det_t <- mean(AST$det[, ti])  # detection uses all donors
keepa <- is.finite(AST$stage)
x <- AST$expr[keepa, ti]; y <- as.numeric(AST$stage[keepa]); z <- AST$gmean[keepa]; n <- length(x)
rho <- pcor_sp(x, y, z); ci <- rho_ci(rho, n, 1L); p <- rho_p(rho, n, 1L); bci <- boot_rho(x, y, z)
cat(sprintf("  detection rate %.4f %s\n", det_t,
            if (det_t < DET_MIN) sprintf("*** BELOW the %.2f floor -- see reading rule", DET_MIN) else ""))
cat(sprintf("  rho = %+.3f   p = %.3g   Fisher 95%% CI [%+.3f, %+.3f]   bootstrap [%+.3f, %+.3f]   n = %d\n",
            rho, p, ci[1], ci[2], bci[1], bci[2], n))
cat(sprintf("  SEA-AD MTG reference: rho = %+.3f (n = 84)\n", SEAAD_RHO))
cat("\n  identity genes, for context (a global collapse would move these too):\n")
IDT <- rbindlist(lapply(intersect(ANCHOR2, AST$genes), function(g) {
  j <- match(g, AST$genes); r <- pcor_sp(AST$expr[keepa, j], y, z)
  c2 <- rho_ci(r, n, 1L)
  data.table(gene = g, rho = round(r, 3), lo = round(c2[1], 3), hi = round(c2[2], 3),
             detection = round(mean(AST$det[keepa, j]), 4)) }))
print(IDT, row.names = FALSE); fwrite(IDT, file.path(OUTDIR, "R1_identity_genes.csv"))

## =============================================================================
## R2 -- detection-matched null
## =============================================================================
cat("\n=== R2. DETECTION-MATCHED NULL ===\n")
gdet <- colMeans(AST$det[keepa, , drop = FALSE])
pool <- which(gdet >= DET_MIN & apply(AST$expr[keepa, , drop = FALSE], 2, sd) > 0)
band <- pool[gdet[pool] >= det_t * BAND_LO & gdet[pool] <= det_t * BAND_HI]
band <- setdiff(band, ti)
if (length(band) < 100) {
  cat(sprintf("  only %d matched genes; the null is not reported.\n", length(band)))
  pctl <- NA_real_
} else {
  rb <- vapply(band, function(j) suppressWarnings(pcor_sp(AST$expr[keepa, j], y, z)), numeric(1))
  rb <- rb[is.finite(rb)]
  pctl <- 100 * mean(rb <= rho)
  cat(sprintf("  %d genes matched on detection (band %.3f-%.3f)\n", length(rb), det_t*BAND_LO, det_t*BAND_HI))
  cat(sprintf("  MCT4 rho %+.3f at the %.2fth percentile | band median %+.3f\n", rho, pctl, median(rb)))
  fwrite(data.table(gene = c(TARGET, AST$genes[band][is.finite(vapply(band, function(j)
           suppressWarnings(pcor_sp(AST$expr[keepa, j], y, z)), numeric(1)))]),
         rho = c(rho, rb)), file.path(OUTDIR, "R2_detection_null.csv"))
}

## =============================================================================
## R3 -- ambient control against microglia
## =============================================================================
cat(sprintf("\n=== R3. AMBIENT CONTROL (astrocytes vs microglia, %s) ===\n", REGION))
if (is.null(MIC)) cat("  microglia not available; ambient control not run.\n") else {
  mj <- match(TARGET, MIC$genes)
  sh <- intersect(AST$donors, MIC$donors)
  if (is.na(mj) || length(sh) < MIN_DONOR) cat("  insufficient overlap; ambient control not run.\n") else {
    ia <- match(sh, AST$donors); im <- match(sh, MIC$donors)
    ra <- pcor_sp(AST$expr[ia, ti], as.numeric(AST$stage[ia]), AST$gmean[ia])
    rm_ <- pcor_sp(MIC$expr[im, mj], as.numeric(MIC$stage[im]), MIC$gmean[im])
    AM <- data.table(cell_type = c("astrocyte","microglia"),
                     mean_expr = round(c(mean(AST$expr[ia, ti]), mean(MIC$expr[im, mj])), 4),
                     detection = round(c(mean(AST$det[ia, ti]), mean(MIC$det[im, mj])), 4),
                     rho_vs_stage = round(c(ra, rm_), 3), n_donors = length(sh))
    print(AM, row.names = FALSE); fwrite(AM, file.path(OUTDIR, "R3_ambient.csv"))
    cat("  Ambient bleed-through would move the highest-expressing population most.\n")
  }
}

## =============================================================================
## VERDICT -- rule fixed before the run
## =============================================================================
cat("\n=== VERDICT ===\n")
call_verdict <- function(g1ok, rho, ci, pctl) {
  if (!isTRUE(g1ok))                       return("UNTRUSTED -- G1 gate failed")
  if (!is.finite(rho))                     return("NOT TESTABLE")
  if (rho <= MDE && ci[2] < 0 && (is.na(pctl) || pctl <= NULL_PCTL))
                                           return("REPLICATED")
  if (rho <= MDE && ci[2] < 0)             return("REPLICATED IN DIRECTION AND SIZE, BUT NOT SPECIFIC TO ITS DETECTION RATE")
  if (rho < 0 && ci[2] >= 0)               return("DIRECTION CONSISTENT, UNDERPOWERED -- NOT a failure to replicate")
  if (ci[1] > MDE)                         return("NOT REPLICATED")
  "INCONCLUSIVE"
}
verdict <- call_verdict(g1ok, rho, ci, pctl)
if (!g1ok) verdict <- sprintf("%s | had the gate passed, the pre-set rule would return: %s",
                              verdict, call_verdict(TRUE, rho, ci, pctl))
cat(sprintf("  %s\n", verdict))
cat(sprintf("\n  Reading rule, fixed in advance: n = 31 detects |rho| >= %.2f at 80%% power.\n", abs(MDE)))
cat("  A null therefore excludes a SEA-AD-sized effect but not a smaller one, and is\n")
cat("  reported as an interval. 'Not replicated' is written only when the interval\n")
cat("  excludes the SEA-AD estimate. Neurons are absent from this dataset, so the\n")
cat("  cross-cellular coupling is not tested here and no claim about it is made.\n")
fwrite(data.table(target = TARGET, n_donors = n, rho = round(rho, 4),
                  ci_lo = round(ci[1], 4), ci_hi = round(ci[2], 4), p = signif(p, 3),
                  boot_lo = round(bci[1], 4), boot_hi = round(bci[2], 4),
                  detection = round(det_t, 4), null_percentile = round(pctl, 2),
                  seaad_rho = SEAAD_RHO, g1_pass = g1ok, verdict = verdict),
       file.path(OUTDIR, "VERDICT.csv"))
cat(sprintf("\nWritten to %s\n", OUTDIR))
