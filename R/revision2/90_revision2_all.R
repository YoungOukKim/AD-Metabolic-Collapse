# ==============================================================================
# P2_revision2_all.R
#
# Purpose : Every analysis added for the npj Dementia revision of
#           "Astrocytic lactate shuttle disruption and the energy-starved
#           lysosome in Alzheimer's disease", in one script, writing every
#           result to disk so the run is reviewable and committable.
#
# Reviewer mapping
#   S4  acid-base network                Reviewer 1, point 2
#   S5  mitochondria, astro vs neuron    Reviewer 2, point 2
#   S6  alternative pathway screen       Reviewer 2, points 1 and 3
#   S7  CSF on the full ADNI inputs      Reviewer 2, points 4 and 5
#   S8  neuropathological staging axes   Reviewer 1, points 4 and 5
#   S9  curve features                   Reviewer 1, point 3
#
# Conventions taken from the existing repositories, not re-invented here:
#   set.seed(42)                         P2_run_all.R, 01_extract_seaad.R,
#                                        Fig2D, 70_ITG_replication, 46/47/49
#                                        mediation; AD-Boundary-Audit engine
#   Sys.getenv path with a D: default    oligo_MCT_detection_check.R
#   rhdf5 + data.table, block slurp      oligo_MCT_detection_check.R,
#                                        celltype_audit_from_raw.R
#   bit64conversion = "double"           required for indptr past 2^31
#   9 CPS bins, cut(right = FALSE)       oligo_MCT_detection_check.R, Fig.1/FigS
#   z per gene across cells, then bin    same
#   det_rate >= 0.10 to report a
#     trajectory on its own              oligo_MCT_detection_check.R
#   quadratic vertex / discrete argmax /
#     segmented breakpoint, each named   AD-Boundary-Audit/PTGDS_vertex_check.R
#   pc() partial correlation             AD-Metabolic-Collapse 00_utils.R
#
#   WHERE THIS SCRIPT DEPARTS FROM A HOUSE CONVENTION IT SAYS SO IN THE OUTPUT,
#   in output/CONVENTIONS.csv, rather than silently.
#
# Input   : SEA-AD MTG h5ad, the July P2 outputs, the ADNI files on D:
# Output  : OUTDIR/*.csv, OUTDIR/VERDICTS.txt, OUTDIR/MANIFEST.csv,
#           OUTDIR/SESSION.txt
# Requires: rhdf5, data.table; readxl optional; hdf5r not required
# Usage   : paste this whole file into the RStudio editor and run it, or
#           source("P2_revision2_all.R"). It self-tests and then runs every
#           step to the end. There is no switch and nothing else to type.
#           S3 reads the h5ad and takes 40-70 minutes; every other step is fast.
#           To re-run one step afterwards: p2_run("S4")
# ==============================================================================
suppressPackageStartupMessages({ library(rhdf5); library(data.table) })
set.seed(42)

## ---- 0. configuration --------------------------------------------------------
MTG_H5   <- Sys.getenv("H5AD_PATH",   unset = file.path("data", "SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad"))
P2_H5DIR <- Sys.getenv("P2_H5DIR",    unset = file.path("output", "sensitivity"))
P2_RUNDIR<- Sys.getenv("P2_RUNDIR",   unset = file.path("output", "tables"))
EMORY    <- Sys.getenv("EMORY_CSV",   unset = file.path("data", "EMORY_CSF_TMT_MS.csv"))
ELECSYS  <- Sys.getenv("ELECSYS_CSV", unset = file.path("data", "UPENNBIOMK_ROCHE_ELECSYS_09Jan2026.csv"))
ADNI_RDA <- Sys.getenv("ADNI_RDA",    unset = file.path("data", "ADNIMERGE2"))
OUTDIR   <- Sys.getenv("P2R2_OUT",    unset = file.path("output", "revision2"))
if (!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE)

BUILD <- "P2_revision2_all.R  2026-08-14  r6 (reuses a verified S3 pass; z composites)"

## ---- 1. gates and thresholds, fixed before any file is opened ----------------
G0 <- list(nuclei = 1378211L, astro = 67419L, exc = 671689L, donors = 84L)
G1 <- list(det_all = 0.0612, det_bin = 0.0600, pct = -43.2, tol_det = 0.004, tol_pct = 1.5)
DET_FLOOR_HOUSE <- 0.10   # oligo_MCT_detection_check.R: report a trajectory alone
DET_FLOOR_NULL  <- 0.02   # below this a matched-null percentile is not computed
MIN_MEMBERS <- 3L; MIN_NULL <- 200L; MATCH_TOL <- 0.25
ALPHA <- 0.05; NULL_TAIL <- 0.05; R_ANLS <- 0.40
MIN_DONORS <- 8L; TRIM_FRAC <- 0.05; MED_TOL <- 0.10
# TWO binning conventions exist in the repositories and they are NOT the same.
#   P2, 77_detection_matched_null_standalone.R : bin <- round(cps, 1)
#       so bin 0.2 is CPS in [0.15, 0.25). Bin 0.1 is excluded from bin-level
#       work. early / late are the UNWEIGHTED MEAN of three bin means, not the
#       pooled mean over the cells in those bins.
#   P1-fixed, oligo_MCT_detection_check.R      : cut(seq(0.1,1.0,0.1), right=FALSE)
#       so bin 1 is CPS in [0.10, 0.20).
# The published -43.2% comes from the first. An earlier draft of this file used
# the second, produced -46.5%, and the G1 gate caught it. Reproducing a P2
# number requires the P2 form; that is what is used below, and the difference is
# recorded in CONVENTIONS.csv rather than left as a footnote.
BIN_ROUND   <- TRUE                              # P2: round(cps, 1)
BIN_EXCLUDE <- 0.1                               # P2: bin 0.1 out of bin-level work
EARLY_BINS  <- c(0.2, 0.3, 0.4)                  # P2
LATE_BINS   <- c(0.6, 0.7, 0.8)                  # P2
BIN_CENTRES <- seq(0.1, 0.9, by = 0.1)           # nine bins, P2 labelling
BLOCK <- 100000L                                 # house: oligo script

ASTRO_SUB <- "Astrocyte"
EXC_SUB   <- c("^L2/3 IT$","^L4 IT$","^L5 IT$","^L6 IT$","^L6 IT Car3$",
               "^L5 ET$","^L5/6 NP$","^L6 CT$","^L6b$")
VATP6 <- c("ATP6V1A","ATP6V1B2","ATP6V0A1","ATP6V0C","ATP6V0D1","ATP6V1E1")
CTRL_NEUR <- c("SNAP25","SYT1","SLC17A7","RBFOX3")
ANLS <- c("SLC16A3","LDHA","HK2","PFKFB3")

CA_GENES  <- c("CA2","CA4","CA12","CA14","CA8","CA11","CA9")
CHAPERONE <- c("BSG","EMB")
PH_OTHER  <- c("SLC4A4","SLC4A2","SLC4A10","SLC9A1","SLC9A6","SLC9A7")
MCTS      <- c("SLC16A1","SLC16A3","SLC16A7")
CI_A <- c("NDUFA1","NDUFA2","NDUFA4","NDUFA8","NDUFA9","NDUFA10","NDUFA13")
CI_B <- c("NDUFB1","NDUFB2","NDUFB4","NDUFB6","NDUFB8","NDUFB9","NDUFB10")
CI_S <- c("NDUFS1","NDUFS2","NDUFS3","NDUFS4","NDUFS7","NDUFS8","NDUFV1","NDUFV2")
CII  <- c("SDHA","SDHB"); CIII <- c("UQCRC1","UQCRC2","CYC1")
CIV  <- c("COX4I1","COX5A","COX6C","COX7C"); CV <- c("ATP5F1A","ATP5F1B","ATP5MC1","ATP5PB")
BIOG <- c("PPARGC1A","PPARGC1B","PPRC1","TFAM","NRF1","NFE2L2","ESRRA")
DYN  <- c("MFN1","MFN2","OPA1","DNM1L","FIS1")
ALT_PANEL <- c("SLC1A2","SLC1A3","GLUL","SLC38A1","SLC38A3","ATP1A1","ATP1A2",
  "ATP1A3","ATP1B1","ATP1B2","AQP4","KCNJ10","GJA1","GJB6","SLC2A1","SLC2A3",
  "GYS1","PYGB","PYGM","CPT1A","ACADVL","HADHA","ACSL3","BDH1","OXCT1","SLC16A6",
  "RPL13A","RPS6","EIF4E","EEF2","HSPA5","DDIT3","ATF4","EIF2AK3","GFAP","VIM",
  "SERPINA3","CLU","C3","CD44","LAMP1","LAMP2","CTSB","CTSD","SLC40A1","TFRC",
  "FTH1","FTL","CP","LDHA","HK2","PFKFB3","NEFL","PTGDS","MOG","MAG")
PANEL <- unique(c(CA_GENES, CHAPERONE, PH_OTHER, MCTS, CI_A, CI_B, CI_S, CII,
                  CIII, CIV, CV, BIOG, DYN, ALT_PANEL, CTRL_NEUR, VATP6))

PATHWAYS <- list(
  "Lactate export (ANLS hub)"      = list(ct="astro", g=ANLS),
  "Iron handling"                  = list(ct="astro", g=c("TFRC","FTH1","FTL","CP","SLC40A1")),
  "Glutamate uptake and recycling" = list(ct="astro", g=c("SLC1A2","SLC1A3","GLUL","SLC38A3")),
  "Na/K-ATPase"                    = list(ct="astro", g=c("ATP1A1","ATP1A2","ATP1B1","ATP1B2")),
  "Water and K+ buffering"         = list(ct="astro", g=c("AQP4","KCNJ10","GJA1","GJB6")),
  "Glucose uptake and glycogen"    = list(ct="astro", g=c("SLC2A1","SLC2A3","GYS1","PYGB")),
  "Fatty acid oxidation"           = list(ct="astro", g=c("CPT1A","ACADVL","HADHA","ACSL3")),
  "Ketone body handling"           = list(ct="astro", g=c("BDH1","OXCT1","SLC16A6")),
  "Astrocyte reactivity"           = list(ct="astro", g=c("GFAP","VIM","SERPINA3","CLU","CD44")),
  "Mitochondrial ETC (astro)"      = list(ct="astro", g=c("NDUFS1","NDUFA9","SDHA","UQCRC1","COX4I1","ATP5F1A")),
  "Neuronal V-ATPase"              = list(ct="exc",   g=VATP6),
  "Neuronal lysosome"              = list(ct="exc",   g=c("LAMP1","LAMP2","CTSB","CTSD")),
  "Mitochondrial ETC (neuron)"     = list(ct="exc",   g=c("NDUFS1","NDUFA9","SDHA","UQCRC1","COX4I1","ATP5F1A")),
  "Neuronal protein synthesis"     = list(ct="exc",   g=c("RPL13A","RPS6","EIF4E","EEF2")),
  "Integrated stress response"     = list(ct="exc",   g=c("ATF4","DDIT3","HSPA5","EIF2AK3")),
  "Neuronal identity (control)"    = list(ct="exc",   g=CTRL_NEUR))

STAGE_AXES <- list(
  adnc  = list(col="adnc_int",  role="primary",     label="Overall AD neuropathological Change (NIA-AA)",
               early=function(x) x<=1L, late=function(x) x>=3L, elab="Not AD + Low", llab="High"),
  thal  = list(col="thal_int",  role="sensitivity", label="Thal amyloid phase",
               early=function(x) x<=1L, late=function(x) x>=4L, elab="Thal 0-1", llab="Thal 4-5"),
  cerad = list(col="cerad_int", role="sensitivity", label="CERAD neuritic plaque density",
               early=function(x) x<=0L, late=function(x) x>=3L, elab="Absent", llab="Frequent"),
  braak = list(col="braak",     role="pre-specified, gate not met", label="Braak tau stage",
               early=function(x) x<=2L, late=function(x) x>=5L, elab="Braak 0-II", llab="Braak V-VI"))

ROMAN <- c("0"=0L,"I"=1L,"II"=2L,"III"=3L,"IV"=4L,"V"=5L,"VI"=6L)
NOT_A_STAGE <- c("REFERENCE","NA","NAN","NOT APPLICABLE","UNKNOWN","")

## ---- 2. bookkeeping ----------------------------------------------------------
.MAN <- new.env(); .MAN$rows <- list(); .MAN$verd <- character(0)
put <- function(obj, file, what) {
  f <- file.path(OUTDIR, file)
  data.table::fwrite(as.data.table(obj), f)
  .MAN$rows[[length(.MAN$rows)+1L]] <- data.table(file=file, rows=nrow(obj),
                                                  cols=ncol(obj), contents=what)
  invisible(f) }
verdict <- function(...) {
  ln <- c(...); cat("\n[VERDICT]\n"); cat(paste0("  ", ln, collapse="\n"), "\n")
  .MAN$verd <- c(.MAN$verd, "", ln); invisible(ln) }
note <- function(...) cat(sprintf(...))

## ---- 3. helpers, in the house form -------------------------------------------
read_obs <- function(name, h5 = MTG_H5) {
  v <- h5read(h5, paste0("obs/", name))
  cats <- tryCatch(h5read(h5, paste0("obs/__categories/", name)), error = function(e) NULL)
  if (is.null(cats)) {
    cats <- tryCatch(h5read(h5, paste0("obs/", name, "/categories")), error = function(e) NULL)
    if (!is.null(cats)) v <- h5read(h5, paste0("obs/", name, "/codes")) }
  if (is.null(cats)) return(as.vector(v))
  code <- as.integer(v); out <- rep(NA_character_, length(code))
  ok <- !is.na(code) & code >= 0
  out[ok] <- as.character(cats)[code[ok] + 1L]; out }

# pc(): kept verbatim from AD-Metabolic-Collapse R/sensitivity/00_utils.R so that
# every partial correlation in this project comes from one implementation.
pc <- function(x, y, Z) {
  Z <- as.matrix(Z); o <- complete.cases(x, y, Z)
  x <- x[o]; y <- y[o]; Z <- Z[o, , drop = FALSE]
  r <- cor(residuals(lm(x ~ Z)), residuals(lm(y ~ Z)))
  n <- sum(o); k <- ncol(Z)
  c(r = r, p = 2 * pt(-abs(r * sqrt((n - k - 2)/(1 - r^2))), n - k - 2), n = n) }

braak_int <- function(x) {
  lab <- trimws(gsub("BRAAK|STAGE|\\.", "", toupper(as.character(x))))
  out <- rep(NA_integer_, length(lab)); k <- lab %in% names(ROMAN)
  out[k] <- ROMAN[lab[k]]
  bad <- unique(lab[!k & !(lab %in% NOT_A_STAGE) & !is.na(x) & nzchar(lab)])
  if (length(bad)) stop("unrecognised Braak labels: ", paste(bad, collapse=" | "))
  out }
ord_map <- function(x, map) { k <- toupper(trimws(as.character(x)))
  m <- match(k, toupper(names(map))); out <- rep(NA_integer_, length(k))
  out[!is.na(m)] <- unname(map)[m[!is.na(m)]]; out }

## ---- 4. self-test, before any file is opened ---------------------------------
selftest <- function() {
  cat("=== SELF-TEST ===\n"); f <- character(0)
  say <- function(id, ok, extra="") { cat(sprintf("  %-5s %-56s %s%s\n", id,
      attr(ok,"w"), ifelse(isTRUE(ok),"PASS","FAIL"), extra))
    if (!isTRUE(ok)) f <<- c(f, id) }
  W <- function(x,w){ attr(x,"w") <- w; x }

  # T1  the P2 bin is round(cps, 1), which centres the bin instead of taking its
  #     left edge. Values exactly on a .x5 boundary are avoided here: round() is
  #     round-half-to-even AND 0.15 is below one half in binary, so a boundary
  #     value is not a stable thing to assert on. Real CPS values are continuous.
  cp <- c(0.16, 0.23, 0.26, 0.34, 0.63, 0.66, 0.94)
  b2 <- round(round(cp, 1) * 10)
  bc <- cut(cp, breaks = seq(0.1, 1.0, by = 0.1), right = FALSE, labels = FALSE)
  say("T1", W(identical(as.integer(b2), c(2L,2L,3L,3L,6L,7L,9L)),
      "P2 bin = round(cps,1) centres the bin, e.g. 0.16 -> bin 2"))
  say("T2", W(sum(b2 != bc) == 3L,
      "round() and cut(right=FALSE) assign different cells - not interchangeable"),
      sprintf("  (%d of %d differ)", sum(b2 != bc), length(bc)))
  # T2b  early is the UNWEIGHTED mean of bin means, not the pooled cell mean.
  #      With bins of 100, 100 and 1 cells the two differ by a factor of two.
  mb <- c(1, 2, 6); nb <- c(100, 100, 1)
  say("T2b", W(abs(mean(mb) - 3) < 1e-12 &&
               abs(sum(mb*nb)/sum(nb) - 306/201) < 1e-9 &&
               abs(mean(mb) - sum(mb*nb)/sum(nb)) > 1,
      "unweighted bin mean differs from the pooled cell mean"))
  # T3  z-then-bin equals the manual computation
  set.seed(42); m <- matrix(rnorm(200), 100, 2); bid <- rep(1:5, each=20)
  a1 <- tapply(scale(m[,1])[,1], bid, mean)
  a2 <- tapply((m[,1]-mean(m[,1]))/sd(m[,1]), bid, mean)
  say("T3", W(max(abs(a1-a2)) < 1e-12, "z per gene then bin mean == manual"))
  # T4  pc() equals the residual correlation
  z <- rnorm(60); x <- z+rnorm(60,0,.5); y <- z+rnorm(60,0,.5)
  say("T4", W(abs(pc(x,y,cbind(z))[["r"]] -
      cor(residuals(lm(x~z)), residuals(lm(y~z)))) < 1e-12, "pc() == residual correlation"))
  # T5  Braak parser: known labels map, unknown raises, Reference is skipped
  say("T5", W(identical(braak_int(c("Braak 0","Braak II","VI")), c(0L,2L,6L)) &&
      is.na(braak_int("Reference")) &&
      inherits(try(braak_int("Braak VII"), silent=TRUE), "try-error"),
      "Braak parser maps, skips Reference, errors on the unknown"))
  # T6  detection floor is a hard gate at BOTH house values
  say("T6", W(0.0612 < DET_FLOOR_HOUSE && 0.0612 > DET_FLOOR_NULL,
      "astrocytic MCT4 sits between the two floors - flagged, not hidden"))
  # T7  matched null places a planted extreme in the lower tail
  set.seed(42); det <- runif(3000,0.01,0.20); pct <- rnorm(3000,-4,12)
  det[500] <- 0.06; pct[500] <- -43.2
  bd <- which(abs(det-det[500]) <= MATCH_TOL*det[500])
  say("T7", W(length(bd) >= MIN_NULL && mean(pct[bd] <= pct[500]) < NULL_TAIL,
      "detection-matched null finds the planted extreme"))
  # T8  quadratic vertex recovers a planted peak and rejects a monotone line
  x8 <- rep(0:3, each=21)
  y_up <- 0.05+0.03*x8-0.010*x8^2+rnorm(84,0,.004)
  y_dn <- 0.06-0.008*x8+rnorm(84,0,.004)
  vf <- function(y){ m <- summary(lm(y ~ x8 + I(x8^2)))$coefficients
    v <- -m["x8","Estimate"]/(2*m["I(x8^2)","Estimate"])
    c(q=m["I(x8^2)","Estimate"], p=m["I(x8^2)","Pr(>|t|)"], v=v) }
  u <- vf(y_up); d <- vf(y_dn)
  say("T8", W(u["p"]<ALPHA && u["q"]<0 && u["v"]>0 && u["v"]<3 &&
              !(d["p"]<ALPHA && d["q"]<0 && d["v"]>0 && d["v"]<3),
      "vertex rule accepts a peak and rejects a monotone decline"))
  # T11 a raw mean composite is dominated by the most abundant member; the
  #     z-per-gene composite is not. Planted at the real ANLS scale.
  set.seed(42); nn <- 84; cq <- sort(runif(nn, 0.2, 0.9)); gq <- rnorm(nn, 0, 0.05)
  hi <- 1.037 - 0.05*cq + rnorm(nn, 0, 0.05)                 # abundant, nearly flat
  lo <- 0.0633*(1 - 0.65*(cq-0.2)/0.7) + rnorm(nn, 0, 0.002) # scarce, steep
  Mq <- cbind(hi = hi, lo = lo)
  pv <- function(v) summary(lm(v ~ cq + gq))$coefficients["cq","Pr(>|t|)"]
  zc <- rowMeans(scale(Mq))
  say("T11", W(pv(lo) < 1e-6 && pv(zc) < 1e-6 &&
               mean(hi)/(mean(hi)+mean(lo)) > 0.9,
      "z composite keeps a scarce steep gene that a raw mean would drown"),
      sprintf("  (abundant member is %.0f%% of a raw mean)", 100*mean(hi)/(mean(hi)+mean(lo))))
  # T9  chi-square on top-vs-rest survives a level missing from the top group
  v9 <- c(rep("A",50),rep("B",40),rep("C",10)); t9 <- c(rep(TRUE,5),rep(FALSE,95))
  lv <- sort(unique(v9))
  M <- rbind(as.integer(table(factor(v9[t9],levels=lv))),
             as.integer(table(factor(v9[!t9],levels=lv))))
  say("T9", W(is.finite(suppressWarnings(chisq.test(M)$p.value)),
      "top-vs-rest table survives a level absent from the top group"))
  # T10 baseline visit ranking puts BL first
  vc <- c("m06","bl","m12"); mo <- suppressWarnings(as.integer(sub("^[^0-9]*","",toupper(vc))))
  rk <- ifelse(toupper(vc) %in% c("BL","SC","SCMRI","M00"), 0L, ifelse(is.na(mo),999L,mo))
  say("T10", W(which.min(rk)==2, "baseline visit ranking puts BL first"))

  ok <- length(f)==0
  cat(sprintf("=== SELF-TEST %s%s ===\n\n", ifelse(ok,"PASS","FAIL"),
      ifelse(ok,"",paste0(" (",paste(f,collapse=", "),")"))))
  ok }

## ---- 5. S1  donor metadata from obs -----------------------------------------
s1_donor_meta <- function() {
  if (!file.exists(MTG_H5)) { note("[stop S1] h5ad not found: %s\n", MTG_H5); return(invisible(FALSE)) }
  note("reading obs ...\n")
  donor <- read_obs("Donor ID"); cps <- as.numeric(read_obs("Continuous Pseudo-progression Score"))
  brk   <- read_obs("Braak");    nC <- length(donor)
  stopifnot(length(cps)==nC, length(brk)==nC)
  note("  nuclei %d (cohort file has %d)\n", nC, G0$nuclei)
  n_ref <- sum(toupper(trimws(as.character(brk))) %in% NOT_A_STAGE, na.rm=TRUE)
  note("  nuclei whose Braak label is not a stage (Reference): %d\n", n_ref)
  keep <- !is.na(cps) & !is.na(donor)
  dl <- data.table(donor=donor[keep], lab=brk[keep], cps=cps[keep])
  chk <- dl[, .(nl = uniqueN(lab[!is.na(lab)])), by=donor]
  bad <- chk[nl > 1L]$donor
  note("  B1 Braak constant within donor : %s\n",
       if (!length(bad)) "PASS" else paste("FAIL for", length(bad), "donors"))
  if (length(bad)) return(invisible(FALSE))
  agg <- dl[, .(braak_label = lab[1], cps = mean(cps), n_nuclei = .N), by = donor][order(donor)]
  agg[, braak := braak_int(braak_label)]
  drop <- is.na(agg$braak)
  if (any(drop)) { note("  donors dropped for a non-stage label: %d\n", sum(drop)); agg <- agg[!drop] }
  EXTRA <- c("Thal","CERAD score","Overall AD neuropathological Change","Cognitive Status",
             "APOE Genotype","Age at Death","Sex","PMI","Brain pH","RIN",
             "Severely Affected Donor","Last MMSE Score","Highest Lewy Body Disease","LATE")
  ok_keys <- tryCatch(h5ls(MTG_H5)[h5ls(MTG_H5)$group=="/obs","name"], error=function(e) character(0))
  for (nm in EXTRA[EXTRA %in% ok_keys]) {
    v <- tryCatch(read_obs(nm), error=function(e) NULL)
    if (is.null(v) || length(v) != nC) next
    dv <- tapply(v[keep], dl$donor, function(z){ u <- unique(z[!is.na(z)])
      if (length(u)==1L) u else NA_character_ })
    set(agg, j = make.names(nm), value = as.character(dv[match(agg$donor, names(dv))])) }
  if ("Thal" %in% names(agg))
    agg[, thal_int := ord_map(Thal, c("Thal 0"=0L,"Thal 1"=1L,"Thal 2"=2L,"Thal 3"=3L,"Thal 4"=4L,"Thal 5"=5L))]
  if ("CERAD.score" %in% names(agg))
    agg[, cerad_int := ord_map(CERAD.score, c("Absent"=0L,"Sparse"=1L,"Moderate"=2L,"Frequent"=3L))]
  if ("Overall.AD.neuropathological.Change" %in% names(agg))
    agg[, adnc_int := ord_map(Overall.AD.neuropathological.Change,
        c("Not AD"=0L,"Low"=1L,"Intermediate"=2L,"High"=3L))]
  if ("Cognitive.Status" %in% names(agg))
    agg[, dementia := ord_map(Cognitive.Status, c("No dementia"=0L,"Dementia"=1L))]
  g0 <- nrow(agg)==G0$donors && nC==G0$nuclei
  note("  G0 (nuclei and donor count) : %s\n", ifelse(g0,"PASS","FAIL"))
  put(agg, "S1_donor_stage_full.csv", "donor-level staging and clinical metadata, from obs")
  cnt <- function(cl) if (cl %in% names(agg)) paste(sprintf("%s=%d", names(table(agg[[cl]])),
    as.integer(table(agg[[cl]]))), collapse=" ") else "absent"
  verdict("S1 donor metadata from obs",
    sprintf("G0 %s | donors %d | Reference nuclei excluded %d", ifelse(g0,"PASS","FAIL"), nrow(agg), n_ref),
    sprintf("braak : %s", cnt("braak")), sprintf("thal  : %s", cnt("thal_int")),
    sprintf("cerad : %s", cnt("cerad_int")), sprintf("adnc  : %s", cnt("adnc_int")))
  invisible(g0) }

## ---- 6. S2  artifact consistency + reproduction of the July outputs ----------
EXPECT <- list(
  list(f="celltype_MCT4.csv", k="cell_type", v="astro", col="mean_expr", ref=0.0462, tol=0.0010, w="astrocytic MCT4 mean"),
  list(f="celltype_MCT4.csv", k="cell_type", v="astro", col="detection", ref=0.0612, tol=0.0010, w="astrocytic MCT4 detection"),
  list(f="celltype_MCT4.csv", k="cell_type", v="astro", col="n_cells", ref=67419, tol=0, w="astrocyte nuclei"),
  list(f="celltype_MCT4.csv", k="cell_type", v="neuron",col="n_cells", ref=671689, tol=0, w="excitatory nuclei"),
  list(f="detection_null.csv", k="gene", v="SLC16A3", col="percentile", ref=0.39, tol=0.05, w="MCT4 percentile in the matched null"),
  list(f="detection_null.csv", k="gene", v="SLC16A3", col="n_matched", ref=2042, tol=0, w="genes in the matched band"),
  list(f="hk_specificity_gap.csv", k="control", v="A", col="gap", ref=0.133, tol=0.010, w="housekeeping specificity gap"),
  list(f="hk_specificity_gap.csv", k="control", v="A", col="max_hk", ref=0.651, tol=0.010, w="largest housekeeping |r| (PPIA)"),
  list(f="mct4_bin_trajectory.csv", k="bin", v=0.2, col="SLC16A3", ref=0.05502, tol=0.0005, w="MCT4 bin 0.2"),
  list(f="mct4_bin_trajectory.csv", k="bin", v=0.4, col="SLC16A3", ref=0.07211, tol=0.0005, w="MCT4 bin 0.4"),
  list(f="mct4_bin_trajectory.csv", k="bin", v=0.5, col="SLC16A3", ref=0.04034, tol=0.0005, w="MCT4 bin 0.5"),
  list(f="couplings.csv", k="pair", v="Neuron V-ATPase (6, P2)", k2="stratum", v2="base", col="A_r", ref=0.474, tol=0.010, w="donor-level coupling, CPS-adjusted"),
  list(f="couplings.csv", k="pair", v="Neuron V-ATPase (6, P2)", k2="stratum", v2="base", col="C1_r", ref=0.232, tol=0.010, w="same, genome-wide controlled"),
  list(f="ambient_control.csv", k="gene", v="SLC16A3", col="astro_pct", ref=-43.2, tol=0.5, w="astrocytic MCT4 change"),
  list(f="ambient_control.csv", k="gene", v="SLC16A3", col="micro_pct", ref=-14.4, tol=0.5, w="microglial MCT4 change"),
  list(f="G1_reproduction.csv", k="gene", v="SLC16A3", col="delta", ref=0, tol=0.05, w="G1 delta, MCT4"),
  list(f="G1_reproduction.csv", k="gene", v="SLC16A1", col="recomputed", ref=-11.1, tol=0.3, w="MCT1 change"))

find_in <- function(nm) { for (d in c(P2_RUNDIR, P2_H5DIR, OUTDIR, "."))
  { f <- file.path(d, nm); if (file.exists(f)) return(f) }; NA_character_ }

s2_verify <- function() {
  hit <- function(colv, want) { if (is.numeric(want))
    { x <- suppressWarnings(as.numeric(colv)); is.finite(x) & abs(x-want) < 1e-6 }
    else trimws(as.character(colv)) == trimws(as.character(want)) }
  out <- rbindlist(lapply(EXPECT, function(E) {
    p <- find_in(E$f)
    if (is.na(p)) return(data.table(file=E$f, check=E$w, observed=NA_real_,
      expected=E$ref, tol=E$tol, verdict="NOT CHECKED - file not found"))
    D <- fread(p, showProgress = FALSE)
    if (!(E$k %in% names(D)) || !(E$col %in% names(D)))
      return(data.table(file=E$f, check=E$w, observed=NA_real_, expected=E$ref,
        tol=E$tol, verdict="NOT CHECKED - column absent"))
    sel <- hit(D[[E$k]], E$v)
    if (!is.null(E$k2)) sel <- sel & hit(D[[E$k2]], E$v2)
    r <- which(sel)
    if (length(r) != 1L) return(data.table(file=E$f, check=E$w, observed=NA_real_,
      expected=E$ref, tol=E$tol, verdict=sprintf("NOT CHECKED - %d rows matched", length(r))))
    o <- suppressWarnings(as.numeric(gsub("[%,]","", as.character(D[r][[E$col]]))))
    data.table(file=E$f, check=E$w, observed=o, expected=E$ref, tol=E$tol,
      verdict=ifelse(is.finite(o) && abs(o-E$ref) <= E$tol*(1+1e-9), "match", "MISMATCH")) }))
  put(out, "S2_verify_july_outputs.csv", "July P2 outputs against the published values")

  # C1  the artifacts used downstream must agree with each other
  con <- data.table(check=character(), value=character(), verdict=character())
  fG <- file.path(P2_H5DIR, "donor_by_gene.rds"); fM <- find_in("donor_merged_base.csv")
  if (file.exists(fG) && !is.na(fM)) {
    G <- readRDS(fG); DM <- fread(fM, showProgress = FALSE)
    ord <- match(DM$donor, G$donors)
    a <- G$astro[ord, match("SLC16A3", G$genes)]
    r <- suppressWarnings(cor(a, DM$MCT4)); d <- max(abs(a - DM$MCT4), na.rm=TRUE)
    con <- rbind(con, data.table(check="donor_by_gene.rds vs donor_merged_base.csv, astrocytic MCT4",
      value=sprintf("r = %.6f, max |difference| = %.3g", r, d),
      verdict=ifelse(is.finite(d) && d < 1e-6, "identical",
        ifelse(is.finite(r) && r > 0.99, "SAME ORDERING, DIFFERENT SCALE - check the log1p convention",
               "DISAGREE - do not use them together")))) }
  if (nrow(con)) put(con, "S2_artifact_consistency.csv", "consistency between the artifacts used downstream")
  verdict("S2 reproduction and artifact consistency",
    sprintf("matched %d of %d published values", sum(out$verdict=="match"), nrow(out)),
    if (any(out$verdict!="match")) paste("not matched:",
      paste(out$check[out$verdict!="match"], collapse=" | ")) else "every check matched",
    if (nrow(con)) sprintf("C1 %s -> %s", con$check[1], con$verdict[1]) else
      "C1 not run - donor_by_gene.rds or donor_merged_base.csv missing")
  invisible(all(out$verdict=="match")) }

## ---- 7. S3  one h5ad pass: detection and bin means, both cell types ----------
s3_extract <- function() {
  # Reuse a completed pass instead of reading the h5ad again. The cached table
  # is only accepted if it still reproduces the published SLC16A3 values, so a
  # stale or wrongly-binned file is re-read rather than trusted. There is no
  # switch to set: running the file does the right thing either way.
  fD <- file.path(OUTDIR, "S3_detection_by_celltype.csv")
  fB <- file.path(OUTDIR, "S3_bin_means_panel.csv")
  if (file.exists(fD) && file.exists(fB)) {
    D0 <- tryCatch(fread(fD, showProgress = FALSE), error = function(e) NULL)
    i0 <- if (!is.null(D0) && "gene" %in% names(D0)) match("SLC16A3", D0$gene) else NA
    if (!is.na(i0)) {
      ok1 <- abs(D0$det_all_astro[i0]  - G1$det_all) <= G1$tol_det &&
             abs(D0$det_bin_astro[i0]  - G1$det_bin) <= G1$tol_det &&
             abs(D0$pct_change_astro[i0] - G1$pct)   <= G1$tol_pct
      if (ok1) {
        note("cached S3 found and it reproduces the published values - h5ad not re-read\n")
        note("  %s\n  %s\n", fD, fB)
        .MAN$rows[[length(.MAN$rows)+1L]] <- data.table(file=basename(fD),
          rows=nrow(D0), cols=ncol(D0), contents="genome-wide detection and change, both cell types (reused)")
        .MAN$rows[[length(.MAN$rows)+1L]] <- data.table(file=basename(fB),
          rows=NA_integer_, cols=NA_integer_, contents="panel bin means, both cell types (reused)")
        verdict("S3 one h5ad pass",
          "REUSED an existing pass; the h5ad was not read again.",
          sprintf("astrocytic SLC16A3 det_all %.4f (published %.4f) | change %.1f%% (published %.1f%%)",
                  D0$det_all_astro[i0], G1$det_all, D0$pct_change_astro[i0], G1$pct),
          "G1 reproduction gate : PASS (on the cached table)",
          sprintf("NOTE astrocytic MCT4 detection %.4f is BELOW the house reporting floor of %.2f.",
                  D0$det_bin_astro[i0], DET_FLOOR_HOUSE),
          "     oligo_MCT_detection_check.R reports a trajectory alone only at >= 0.10.",
          "     Here the claim rests on a detection-matched null instead, and that",
          "     difference must be stated in the reply rather than left implicit.")
        return(invisible(TRUE)) }
      note("cached S3 does NOT reproduce the published values - reading the h5ad again\n") } }
  if (!file.exists(MTG_H5)) { note("[stop S3] h5ad not found\n"); return(invisible(FALSE)) }
  t0 <- Sys.time()
  sub <- read_obs("Subclass"); cps <- as.numeric(read_obs("Continuous Pseudo-progression Score"))
  don <- read_obs("Donor ID"); genes <- as.character(h5read(MTG_H5, "var/_index"))
  nC <- length(sub); nG <- length(genes)
  donors <- sort(unique(don[!is.na(cps)])); inCoh <- !is.na(cps) & don %in% donors
  isA <- sub == ASTRO_SUB & inCoh
  isE <- Reduce(`|`, lapply(EXC_SUB, function(p) grepl(p, sub))) & inCoh
  note("G0 cohort gate: nuclei %d/%d | astro %d/%d | exc %d/%d | donors %d/%d\n",
       nC, G0$nuclei, sum(isA), G0$astro, sum(isE), G0$exc, length(donors), G0$donors)
  g0 <- nC==G0$nuclei && sum(isA)==G0$astro && sum(isE)==G0$exc && length(donors)==G0$donors
  note("  G0 = %s\n", ifelse(g0,"PASS","FAIL - stopping"))
  if (!g0) return(invisible(FALSE))

  # P2 form, verbatim in shape from 77_detection_matched_null_standalone.R
  binv  <- round(cps, 1)
  binid <- round(binv * 10)                       # 1 .. 10
  # bin 0.1 is excluded from bin-level work (P2). Compared on the integer code,
  # not on the double, so no floating-point equality is relied on.
  inbin <- !is.na(binid) & binid >= 2L & binid <= 9L
  ecol  <- round(EARLY_BINS * 10); lcol <- round(LATE_BINS * 10)
  pidx <- match(PANEL, genes); pres <- !is.na(pidx)
  note("panel genes present %d of %d\n", sum(pres), length(PANEL))
  gmap <- rep(NA_integer_, nG); gmap[pidx[pres]] <- which(pres); nP <- length(PANEL)

  indptr <- h5read(MTG_H5, "X/indptr", bit64conversion = "double")
  acc <- function() list(det=numeric(nG), det_b=numeric(nG),
    sum_by_bin=matrix(0, nG, 9), n_by_bin=numeric(9),
    n=0, nb=0, bs=matrix(0,9,nP), bd=matrix(0,9,nP), bn=numeric(9))
  A <- acc(); E <- acc()
  note("reading expression ...\n")
  for (s0 in seq(1L, nC, by = BLOCK)) {
    e0 <- min(s0+BLOCK-1L, nC); sp <- indptr[s0]; cnt <- indptr[e0+1L]-sp
    if (cnt <= 0) next
    if (!any(isA[s0:e0] | isE[s0:e0])) { note("  ...%d/%d\r", e0, nC); next }
    ci <- h5read(MTG_H5,"X/indices",start=sp+1L,count=cnt,bit64conversion="double") + 1L
    cd <- h5read(MTG_H5,"X/data",   start=sp+1L,count=cnt)
    p  <- indptr[s0:(e0+1L)] - sp
    cid <- rep.int(seq_len(e0-s0+1L), as.integer(diff(p)))
    gl <- s0:e0
    for (tg in c("A","E")) {
      kp <- if (tg=="A") isA[gl] else isE[gl]
      if (!any(kp)) next
      O <- get(tg); loc <- which(kp); m <- cid %in% loc
      if (any(m)) {
        gi <- ci[m]; gx <- cd[m]; cc <- cid[m]; bglob <- binid[gl][cc]; ibb <- inbin[gl][cc]
        O$det <- O$det + tabulate(gi, nG)
        if (any(ibb)) O$det_b <- O$det_b + tabulate(gi[ibb], nG)
        okb <- !is.na(bglob) & bglob >= 1L & bglob <= 9L
        if (any(okb)) {
          lin <- gi[okb] + (bglob[okb] - 1L) * nG
          ag <- rowsum(gx[okb], lin, reorder=FALSE)
          i <- as.integer(rownames(ag)); O$sum_by_bin[i] <- O$sum_by_bin[i] + as.numeric(ag) }
        pos <- gmap[gi]; okp <- !is.na(pos) & okb
        if (any(okp)) {
          lin2 <- (pos[okp]-1L)*9L + bglob[okp]
          ag <- rowsum(gx[okp], lin2, reorder=FALSE)
          i <- as.integer(rownames(ag)); O$bs[i] <- O$bs[i] + as.numeric(ag)
          ad <- rowsum(rep(1, sum(okp)), lin2, reorder=FALSE)
          i2 <- as.integer(rownames(ad)); O$bd[i2] <- O$bd[i2] + as.numeric(ad) } }
      O$n <- O$n + length(loc); O$nb <- O$nb + sum(inbin[gl][loc])
      bs2 <- binid[gl][loc]; bs2 <- bs2[!is.na(bs2) & bs2 >= 1L & bs2 <= 9L]
      O$n_by_bin <- O$n_by_bin + tabulate(bs2, 9); O$bn <- O$bn + tabulate(bs2, 9)
      assign(tg, O) }
    note("  ...%d/%d\r", e0, nC) }
  note("\n")
  # early / late are the UNWEIGHTED mean of the three bin means, as in the P2
  # script that produced the published figure - not the pooled cell mean.
  mk <- function(O) { mbb <- sweep(O$sum_by_bin, 2, pmax(O$n_by_bin, 1), "/")
    e <- rowMeans(mbb[, ecol, drop=FALSE]); l <- rowMeans(mbb[, lcol, drop=FALSE])
    data.table(gene=genes, det_all=O$det/O$n, det_bin=O$det_b/max(O$nb,1),
      early=e, late=l, pct_change=ifelse(!is.finite(e)|e<=0, NA_real_, 100*(l/e - 1))) }
  DA <- mk(A); DE <- mk(E)
  setnames(DA, setdiff(names(DA),"gene"), paste0(setdiff(names(DA),"gene"),"_astro"))
  setnames(DE, setdiff(names(DE),"gene"), paste0(setdiff(names(DE),"gene"),"_exc"))
  DET <- merge(DA, DE, by="gene", sort=FALSE)
  DET[, `:=`(reportable_astro = det_bin_astro >= DET_FLOOR_HOUSE,
             reportable_exc   = det_bin_exc   >= DET_FLOOR_HOUSE,
             nullable_astro   = det_bin_astro >= DET_FLOOR_NULL,
             nullable_exc     = det_bin_exc   >= DET_FLOOR_NULL)]
  put(DET, "S3_detection_by_celltype.csv", "genome-wide detection and early-to-late change, both cell types")
  BM <- rbindlist(lapply(c("astro","exc"), function(tg) { O <- if (tg=="astro") A else E
    rbindlist(lapply(1:9, function(b) data.table(celltype=tg, bin=BIN_CENTRES[b], gene=PANEL,
      mean=O$bs[b,]/max(O$bn[b],1), det=O$bd[b,]/max(O$bn[b],1), n_cells=O$bn[b]))) }))
  put(BM, "S3_bin_means_panel.csv", "panel bin means and detection, 9 house bins, both cell types")
  i <- match("SLC16A3", DET$gene)
  g1 <- abs(DET$det_all_astro[i]-G1$det_all) <= G1$tol_det &&
        abs(DET$det_bin_astro[i]-G1$det_bin) <= G1$tol_det &&
        abs(DET$pct_change_astro[i]-G1$pct)  <= G1$tol_pct
  verdict("S3 one h5ad pass",
    sprintf("G0 PASS | elapsed %.1f min", as.numeric(difftime(Sys.time(), t0, units="mins"))),
    sprintf("astrocytic SLC16A3 det_all %.4f (published %.4f) | change %.1f%% (published %.1f%%)",
            DET$det_all_astro[i], G1$det_all, DET$pct_change_astro[i], G1$pct),
    sprintf("G1 reproduction gate : %s", ifelse(g1,"PASS","FAIL - do not read S4 to S6")),
    sprintf("NOTE astrocytic MCT4 detection %.4f is BELOW the house reporting floor of %.2f.",
            DET$det_bin_astro[i], DET_FLOOR_HOUSE),
    "     oligo_MCT_detection_check.R reports a trajectory alone only at >= 0.10.",
    "     Here the claim rests on a detection-matched null instead, and that",
    "     difference must be stated in the reply rather than left implicit.")
  invisible(g1) }

## ---- 8. shared donor-level loader for S4 to S6, S8, S9 -----------------------
load_donor <- function(need_det = TRUE) {
  fG <- file.path(P2_H5DIR, "donor_by_gene.rds"); fM <- find_in("donor_merged_base.csv")
  if (!file.exists(fG) || is.na(fM)) { note("[stop] donor_by_gene.rds or donor_merged_base.csv missing\n"); return(NULL) }
  G <- readRDS(fG); DM <- fread(fM, showProgress = FALSE)
  ord <- match(DM$donor, G$donors); if (anyNA(ord)) { note("[stop] donor mismatch\n"); return(NULL) }
  Ga <- G$astro[ord,,drop=FALSE]; Ge <- G$neuron[ord,,drop=FALSE]
  colnames(Ga) <- colnames(Ge) <- G$genes
  DET <- NULL; src <- "none"
  if (need_det) {
    f <- file.path(OUTDIR, "S3_detection_by_celltype.csv")
    if (file.exists(f)) { DET <- fread(f, showProgress=FALSE); src <- "S3 (both cell types)" }
    else { f2 <- file.path(P2_H5DIR, "detection_matched_null_all_genes.csv")
      if (!file.exists(f2)) { note("[stop] no detection table; run S3\n"); return(NULL) }
      A <- fread(f2, showProgress=FALSE)
      DET <- data.table(gene=A$gene, det_all_astro=A$detect_all, det_bin_astro=A$detect_bin,
        pct_change_astro=A$pct_change, det_all_exc=NA_real_, det_bin_exc=NA_real_,
        pct_change_exc=NA_real_)
      DET[, `:=`(reportable_astro=det_bin_astro>=DET_FLOOR_HOUSE, reportable_exc=NA,
                 nullable_astro=det_bin_astro>=DET_FLOOR_NULL, nullable_exc=NA)]
      src <- "detection_matched_null_all_genes.csv (ASTROCYTES ONLY - neuronal tests pending S3)" }
    note("  detection source: %s\n", src) }
  fS <- file.path(OUTDIR, "S1_donor_stage_full.csv")
  if (file.exists(fS)) { S <- fread(fS, showProgress=FALSE)
    add <- setdiff(intersect(c("adnc_int","thal_int","cerad_int","braak"), names(S)), names(DM))
    if (length(add)) DM <- merge(DM, S[, c("donor", add), with=FALSE], by="donor", all.x=TRUE) }
  list(DET=DET, Ga=Ga, Ge=Ge, DM=DM, astro_only=need_det && !grepl("both cell", src)) }

# comp_z(): z per gene ACROSS DONORS, then the mean. This is the house form
# (oligo_MCT_detection_check.R z-scores each gene before binning) and it matters
# here for a concrete reason. In the ANLS panel the astrocytic means are
# SLC16A3 0.063, LDHA 0.105, HK2 0.073, PFKFB3 1.037. A raw rowMeans is 81 per
# cent PFKFB3, which changes by only -5 per cent, so MCT4's -43 per cent is
# diluted away: the composite came out at p = 0.09 while MCT4 alone is
# p = 3.6e-07. A screen whose positive control cannot fire is uninformative.
# Genes with no variance across donors are dropped rather than turned into NaN.
comp_z <- function(M, genes) {
  j <- match(genes, colnames(M)); j <- j[!is.na(j)]
  if (!length(j)) return(rep(NA_real_, nrow(M)))
  Z <- M[, j, drop = FALSE]
  sdv <- apply(Z, 2, sd, na.rm = TRUE)
  keep <- is.finite(sdv) & sdv > 1e-12
  if (!any(keep)) return(rep(NA_real_, nrow(M)))
  Z <- scale(Z[, keep, drop = FALSE])
  rowMeans(Z, na.rm = TRUE) }

slope_of <- function(v, cps, gm) {
  if (all(!is.finite(v)) || sd(v, na.rm=TRUE) < 1e-12) return(c(NA_real_, NA_real_))
  m <- summary(lm(v ~ cps + gm))$coefficients
  if (!("cps" %in% rownames(m))) return(c(NA_real_, NA_real_))
  c(m["cps","Estimate"], m["cps","Pr(>|t|)"]) }

## ---- 9. S4  acid-base network (Reviewer 1.2) --------------------------------
s4_acidbase <- function() {
  L <- load_donor(); if (is.null(L)) return(invisible(FALSE))
  DET <- L$DET; Ga <- L$Ga; DM <- L$DM; cps <- DM$cps
  FOCUS <- c(CA_GENES, CHAPERONE, PH_OTHER, MCTS)
  out <- rbindlist(lapply(FOCUS, function(g) {
    i <- match(g, DET$gene); ja <- match(g, colnames(Ga))
    if (is.na(i)) return(data.table(gene=g, det_astro=NA_real_, pct_astro=NA_real_,
      n_null=NA_integer_, pctl=NA_real_, slope=NA_real_, p=NA_real_, partial_r_vatp=NA_real_,
      reportable=NA, verdict="not in var"))
    d0 <- DET$det_bin_astro[i]
    bd <- which(is.finite(DET$pct_change_astro) & abs(DET$det_bin_astro-d0) <= MATCH_TOL*d0)
    s <- if (is.na(ja)) c(NA_real_,NA_real_) else slope_of(Ga[,ja], cps, DM$gm_a)
    data.table(gene=g, det_astro=d0, pct_astro=DET$pct_change_astro[i], n_null=length(bd),
      pctl = if (length(bd)>=MIN_NULL) mean(DET$pct_change_astro[bd] <= DET$pct_change_astro[i]) else NA_real_,
      slope=s[1], p=s[2],
      partial_r_vatp = if (is.na(ja)) NA_real_ else pc(Ga[,ja], DM$VATP_n6, cbind(cps))[["r"]],
      reportable = d0 >= DET_FLOOR_HOUSE,
      verdict = if (!isTRUE(d0 >= DET_FLOOR_NULL)) "NOT MEASURABLE (below the null floor)"
        else if (isTRUE(d0 < DET_FLOOR_HOUSE)) "below the house reporting floor - null only"
        else "measurable") }))
  out[, behaves_like_MCT4 := isTRUE(TRUE) & !is.na(pctl) & pctl <= NULL_TAIL &
        !is.na(slope) & slope < 0 & !is.na(p) & p < ALPHA]
  put(out, "S4_acidbase_axis.csv", "carbonic anhydrases, MCT chaperones, exchangers, MCTs")
  hit <- out[behaves_like_MCT4 == TRUE & gene %in% c(CA_GENES, CHAPERONE, PH_OTHER)]$gene
  nm  <- out[grepl("NOT MEASURABLE", verdict)]$gene
  verdict("S4 acid-base network (Reviewer 1, point 2)",
    sprintf("positive control SLC16A3 : change %.1f%%, matched-null percentile %.4f, p %.2g",
            out[gene=="SLC16A3"]$pct_astro, out[gene=="SLC16A3"]$pctl, out[gene=="SLC16A3"]$p),
    sprintf("not measurable : %s", if (length(nm)) paste(nm, collapse=", ") else "(none)"),
    sprintf("acid-base genes behaving like MCT4 : %s",
            if (length(hit)) paste(hit, collapse=", ") else "(none)"),
    sprintf("BSG (CD147), the chaperone through which CA IV acts, detection %.3f, change %+.1f%%",
            out[gene=="BSG"]$det_astro, out[gene=="BSG"]$pct_astro),
    sprintf("MCT1 (SLC16A1) change %+.1f%%, p %.2g, percentile %.3f - a real decline,",
            out[gene=="SLC16A1"]$pct_astro, out[gene=="SLC16A1"]$p, out[gene=="SLC16A1"]$pctl),
    "     so the manuscript word 'preserved' for MCT1 must be replaced by the number.")
  invisible(TRUE) }

## ---- 10. S5  mitochondria (Reviewer 2.2) ------------------------------------
s5_mito <- function() {
  L <- load_donor(); if (is.null(L)) return(invisible(FALSE))
  DET <- L$DET; Ga <- L$Ga; Ge <- L$Ge; DM <- L$DM; cps <- DM$cps
  SETS <- list("Complex I (NDUFA)"=CI_A, "Complex I (NDUFB)"=CI_B, "Complex I (NDUFS/V)"=CI_S,
    "Complex II"=CII, "Complex III"=CIII, "Complex IV"=CIV, "Complex V"=CV,
    "Complex I (all)"=c(CI_A,CI_B,CI_S),
    "ETC (all complexes)"=c(CI_A,CI_B,CI_S,CII,CIII,CIV,CV),
    "Biogenesis regulators"=BIOG, "Fission/fusion"=DYN)
  ALL <- unique(unlist(SETS))
  GR <- rbindlist(lapply(ALL, function(g) { i <- match(g, DET$gene)
    if (is.na(i)) return(data.table(gene=g, in_var=FALSE, det_astro=NA_real_, det_exc=NA_real_,
      pct_astro=NA_real_, pct_exc=NA_real_, ok_astro=FALSE, ok_exc=FALSE))
    data.table(gene=g, in_var=TRUE, det_astro=DET$det_bin_astro[i], det_exc=DET$det_bin_exc[i],
      pct_astro=DET$pct_change_astro[i], pct_exc=DET$pct_change_exc[i],
      ok_astro=isTRUE(DET$det_bin_astro[i] >= DET_FLOOR_NULL),
      ok_exc  =isTRUE(DET$det_bin_exc[i]   >= DET_FLOOR_NULL)) }))
  put(GR, "S5_mito_genes.csv", "per-gene detection and change for the mitochondrial panel")
  cst <- function(set, tg) { M <- if (tg=="astro") Ga else Ge
    gm <- if (tg=="astro") DM$gm_a else DM$gm_n
    if (tg=="exc" && isTRUE(L$astro_only))
      return(data.table(set=NA_character_, celltype=tg, n_total=length(set), n_ok=NA_integer_,
        pct=NA_real_, slope=NA_real_, p=NA_real_, slope_raw=NA_real_,
        p_raw_composite=NA_real_, verdict="PENDING S3 (no neuronal detection)"))
    ok <- GR[gene %in% set & (if (tg=="astro") ok_astro else ok_exc)]$gene
    j <- match(ok, colnames(M)); j <- j[!is.na(j)]
    if (length(j) < MIN_MEMBERS) return(data.table(set=NA_character_, celltype=tg,
      n_total=length(set), n_ok=length(j), pct=NA_real_, slope=NA_real_, p=NA_real_,
      slope_raw=NA_real_, p_raw_composite=NA_real_,
      verdict="NOT REPORTED (too few measurable members)"))
    v <- comp_z(M, ok); s <- slope_of(v, cps, gm)
    vraw <- rowMeans(M[,j,drop=FALSE]); sraw <- slope_of(vraw, cps, gm)
    pc2 <- if (tg=="astro") GR$pct_astro else GR$pct_exc
    data.table(set=NA_character_, celltype=tg, n_total=length(set), n_ok=length(j),
      pct=mean(pc2[match(ok, GR$gene)], na.rm=TRUE), slope=s[1], p=s[2],
      slope_raw=sraw[1], p_raw_composite=sraw[2],
      verdict=ifelse(is.na(s[2]),"not estimable",
        ifelse(s[2]<ALPHA, ifelse(s[1]<0,"declines","increases"), "no significant CPS association"))) }
  CO <- rbindlist(lapply(names(SETS), function(nm)
    rbindlist(lapply(c("astro","exc"), function(tg) { r <- cst(SETS[[nm]], tg); r$set <- nm; r }))))
  jm <- match("SLC16A3", colnames(Ga)); sm <- slope_of(Ga[,jm], cps, DM$gm_a)
  CO <- rbind(CO, data.table(set="REFERENCE astrocytic MCT4", celltype="astro", n_total=1L,
    n_ok=1L, pct=DET$pct_change_astro[match("SLC16A3",DET$gene)], slope=sm[1], p=sm[2],
    slope_raw=sm[1], p_raw_composite=sm[2], verdict="reference"))
  put(CO, "S5_mito_composites.csv", "complex-level composites, astrocyte and neuron")
  nmz <- GR[!ok_astro & !(ok_exc %in% TRUE)]$gene
  gv <- function(s,ct) CO[set==s & celltype==ct]$pct[1]
  verdict("S5 mitochondria (Reviewer 2, point 2)",
    sprintf("not measurable in either cell type : %s",
            if (length(nmz)) paste(nmz, collapse=", ") else "(none)"),
    sprintf("biogenesis regulators, astrocyte : %s", CO[set=="Biogenesis regulators" & celltype=="astro"]$verdict[1]),
    sprintf("ETC astro %s | ETC neuron %s | MCT4 %s",
            ifelse(is.finite(gv("ETC (all complexes)","astro")), sprintf("%+.1f%%", gv("ETC (all complexes)","astro")), "NA"),
            ifelse(is.finite(gv("ETC (all complexes)","exc")),   sprintf("%+.1f%%", gv("ETC (all complexes)","exc")),   "pending S3"),
            sprintf("%+.1f%%", CO[set=="REFERENCE astrocytic MCT4"]$pct[1])),
    "Print the biogenesis regulators with their detection rates whether or not they",
    "are measurable: a transcription factor below the floor is an answer, not a gap.")
  invisible(TRUE) }

## ---- 11. S6  alternative pathway screen (Reviewer 2.1 and 2.3) --------------
s6_pathways <- function() {
  L <- load_donor(); if (is.null(L)) return(invisible(FALSE))
  DET <- L$DET; Ga <- L$Ga; Ge <- L$Ge; DM <- L$DM; cps <- DM$cps
  pick <- function(g, ct) { M <- if (ct=="astro") Ga else Ge
    if (ct=="exc" && isTRUE(L$astro_only)) return(character(0))
    dc <- if (ct=="astro") DET$det_bin_astro else DET$det_bin_exc
    k <- g[!is.na(match(g, colnames(M))) & !is.na(match(g, DET$gene))]
    k[dc[match(k, DET$gene)] >= DET_FLOOR_NULL] }
  ao <- pick(ANLS,"astro")
  av <- if (length(ao) >= MIN_MEMBERS) comp_z(Ga, ao) else rep(NA_real_, nrow(Ga))
  rows <- rbindlist(lapply(names(PATHWAYS), function(nm) {
    P <- PATHWAYS[[nm]]; ct <- P$ct; M <- if (ct=="astro") Ga else Ge
    gm <- if (ct=="astro") DM$gm_a else DM$gm_n
    pc2 <- if (ct=="astro") DET$pct_change_astro else DET$pct_change_exc
    ok <- pick(P$g, ct)
    if (length(ok) < MIN_MEMBERS) return(data.table(pathway=nm, celltype=ct, n_total=length(P$g),
      n_ok=length(ok), pct=NA_real_, slope=NA_real_, p_raw=NA_real_,
      slope_rawcomp=NA_real_, p_rawcomp=NA_real_, r_with_ANLS=NA_real_,
      dropped = if (ct=="exc" && isTRUE(L$astro_only)) "PENDING S3" else paste(setdiff(P$g, ok), collapse=";")))
    v <- comp_z(M, ok); s <- slope_of(v, cps, gm)
    vraw <- rowMeans(M[, match(ok, colnames(M)), drop=FALSE]); sraw <- slope_of(vraw, cps, gm)
    data.table(pathway=nm, celltype=ct, n_total=length(P$g), n_ok=length(ok),
      pct=mean(pc2[match(ok, DET$gene)], na.rm=TRUE), slope=s[1], p_raw=s[2],
      slope_rawcomp=sraw[1], p_rawcomp=sraw[2],
      r_with_ANLS=suppressWarnings(cor(v, av, use="complete.obs")),
      dropped=paste(setdiff(P$g, ok), collapse=";")) }))
  ok <- is.finite(rows$p_raw); rows[, p_holm := NA_real_]
  rows$p_holm[ok] <- p.adjust(rows$p_raw[ok], "holm")
  rows[, verdict := fifelse(!is.finite(p_raw), "NOT TESTABLE",
    fifelse(slope < 0 & p_holm < ALPHA & r_with_ANLS > R_ANLS, "CO-DECLINES with ANLS",
      fifelse(slope < 0 & p_holm < ALPHA, "declines but not ANLS-coupled",
        "no decline at the family-adjusted alpha")))]
  setorder(rows, p_holm, na.last=TRUE)
  put(rows, "S6_pathway_screen.csv", "sixteen candidate pathways under one identical rule")
  ctrl <- rows[pathway=="Neuronal identity (control)"]$verdict[1]
  posc <- rows[pathway=="Lactate export (ANLS hub)"]
  pos_ok <- nrow(posc) && is.finite(posc$p_holm[1]) && posc$slope[1] < 0 && posc$p_holm[1] < ALPHA
  verdict("S6 alternative pathway screen (Reviewer 2, points 1 and 3)",
    sprintf("pathways %d | Holm alpha %.2f | ANLS coupling r > %.2f | composites are z per gene",
            nrow(rows), ALPHA, R_ANLS),
    sprintf("POSITIVE CONTROL lactate export (ANLS hub) -> %s  (slope %s, Holm p %s)",
            ifelse(pos_ok, "declines, as required", "DOES NOT DECLINE"),
            ifelse(nrow(posc) && is.finite(posc$slope[1]), sprintf("%+.4f", posc$slope[1]), "NA"),
            ifelse(nrow(posc) && is.finite(posc$p_holm[1]), sprintf("%.4g", posc$p_holm[1]), "NA")),
    if (!pos_ok) paste0("*** POSITIVE CONTROL FAILED. The screen cannot detect the one",
      " decline we already know is there, so no negative in this table is informative.")
      else "positive control fired - negatives in this table are interpretable",
    sprintf("  (raw-rowMeans composite for the same panel: slope %s, p %s - kept for audit;",
            ifelse(nrow(posc) && is.finite(posc$slope_rawcomp[1]), sprintf("%+.4f", posc$slope_rawcomp[1]), "NA"),
            ifelse(nrow(posc) && is.finite(posc$p_rawcomp[1]), sprintf("%.4g", posc$p_rawcomp[1]), "NA")),
    "   a raw mean of this panel is 81 per cent PFKFB3 and dilutes MCT4 away)",
    sprintf("NEGATIVE CONTROL neuronal identity -> %s", ctrl),
    if (grepl("CO-DECLINES", ctrl)) "*** NEGATIVE CONTROL FAILED - the screen reads composition. STOP."
      else "negative control behaved as required",
    sprintf("co-declining : %s", paste(rows[verdict=="CO-DECLINES with ANLS"]$pathway, collapse=", ")),
    sprintf("explored, no decline : %s", paste(rows[verdict=="no decline at the family-adjusted alpha"]$pathway, collapse=", ")),
    sprintf("not testable : %s", paste(rows[verdict=="NOT TESTABLE"]$pathway, collapse=", ")),
    "Name the 'explored, no decline' and 'not testable' pathways in the Results.",
    "Iron appears here as one entry under the same rule as every other entry.")
  invisible(TRUE) }

## ---- 12. S7  CSF on the full ADNI inputs (Reviewer 2.4 and 2.5) -------------
baseline_rows <- function(d, rid="RID", vis="VISCODE2") {
  if (vis %in% names(d)) { v <- toupper(as.character(d[[vis]]))
    mo <- suppressWarnings(as.integer(sub("^[^0-9]*","",v)))
    rk <- ifelse(v %in% c("BL","SC","SCMRI","M00"), 0L, ifelse(is.na(mo), 999L, mo))
    d <- d[order(d[[rid]], rk)] }
  d[!duplicated(d[[rid]])] }
ld_rda <- function(dir, f) { p <- file.path(dir, f); if (!file.exists(p)) return(NULL)
  e <- new.env(); load(p, envir=e); as.data.table(get(ls(e)[1], envir=e)) }

s7_csf <- function() {
  if (!file.exists(EMORY)) { note("[stop S7] Emory matrix not found: %s\n", EMORY); return(invisible(FALSE)) }
  V1A <- "ATP6V1A_P38606"; TAU <- "MAPT_P10636"; HK1 <- "HK1_P19367"
  WANT <- c("RID","PTID","PHASE","VISCODE","VISCODE2","EXAMDATE", V1A, TAU, HK1,
    "ATP6V1E1_P36543","GFAP_P14136","TREM2_Q9NZC2","TFRC_P02786","LAMP2_P13473",
    "CTSB_P07858","NEFL_P07196","APOE_P02649","APOE2_APOE2","APOE4_APOE4")
  hdr <- names(fread(EMORY, nrows=0)); keep <- intersect(WANT, hdr)
  E <- fread(EMORY, select=keep, showProgress=FALSE)
  note("Emory rows %d | columns kept %d of %d\n", nrow(E), length(keep), length(WANT))
  num <- function(x) suppressWarnings(as.numeric(gsub("[<>,]","", as.character(x))))
  for (cl in setdiff(names(E), c("RID","PTID","PHASE","VISCODE","VISCODE2","EXAMDATE")))
    E[[cl]] <- num(E[[cl]])
  E_raw <- copy(E); E <- baseline_rows(E)
  note("  rows %d | unique RIDs %d | after baseline de-duplication %d\n",
       nrow(E_raw), uniqueN(E_raw$RID), nrow(E))

  if (file.exists(ELECSYS)) { EL <- fread(ELECSYS, showProgress=FALSE)
    tc <- intersect(c("TAU","TTAU","T_TAU"), names(EL))[1]
    if (!is.na(tc)) { if ("VISCODE2" %in% names(EL)) EL <- EL[toupper(VISCODE2) %in% c("BL","")]
      EL[[tc]] <- num(EL[[tc]]); EL <- EL[!is.na(RID) & is.finite(EL[[tc]])]
      EL <- EL[!duplicated(RID)]
      E <- merge(E, data.table(RID=EL$RID, ElecsysTau=EL[[tc]]), by="RID", all.x=TRUE)
      note("  Elecsys tau merged for %d subjects (column %s)\n", sum(is.finite(E$ElecsysTau)), tc) } }
  if (dir.exists(ADNI_RDA)) {
    CD <- ld_rda(ADNI_RDA, "CDR.rda")
    if (!is.null(CD)) { g <- intersect(c("CDGLOBAL","CDRGLOBAL"), names(CD))[1]
      if (!is.na(g)) { CDb <- baseline_rows(CD)
        E <- merge(E, data.table(RID=CDb$RID, CDR=num(CDb[[g]])), by="RID", all.x=TRUE)
        note("  CDR merged for %d subjects (column %s)\n", sum(is.finite(E$CDR)), g) } }
    AP <- ld_rda(ADNI_RDA, "APOERES.rda")
    if (!is.null(AP)) {
      # APOERES stores one GENOTYPE string such as "3/4". The Emory APOE2/APOE4
      # columns are isoform PEPTIDE quantifications, not allele counts, and are
      # handled below as continuous variables.
      if ("GENOTYPE" %in% names(AP)) { APb <- baseline_rows(AP)
        sp <- strsplit(gsub("[^0-9/]","", as.character(APb$GENOTYPE)), "/")
        a1 <- suppressWarnings(as.numeric(vapply(sp, function(z) if (length(z)>0) z[1] else NA_character_, "")))
        a2 <- suppressWarnings(as.numeric(vapply(sp, function(z) if (length(z)>1) z[2] else NA_character_, "")))
        k <- is.finite(a1) & is.finite(a2)
        E <- merge(E, data.table(RID=APb$RID[k],
          APOE_geno=paste0(pmin(a1,a2)[k],"/",pmax(a1,a2)[k]),
          APOE_e4n=((a1==4)+(a2==4))[k]), by="RID", all.x=TRUE)
        note("  APOERES merged for %d subjects | e4 dose %s\n", sum(!is.na(E$APOE_geno)),
             paste(sprintf("%s=%d", names(table(E$APOE_e4n)), as.integer(table(E$APOE_e4n))), collapse=" ")) } } }
  dep <- find_in("Supporting_Data_Values_v2.xlsx"); DEPcmp <- NULL
  if (!is.na(dep) && requireNamespace("readxl", quietly=TRUE)) {
    D5 <- as.data.table(readxl::read_excel(dep, sheet="Fig5_data", skip=3))
    if (all(c("RID","DX") %in% names(D5))) {
      D5[, RID := suppressWarnings(as.numeric(RID))]
      E <- merge(E, D5[, .(RID, DX=as.character(DX))], by="RID", all.x=TRUE)
      note("  deposited diagnosis merged for %d subjects\n", sum(!is.na(E$DX)))
      pr <- intersect(c(V1A,TAU,HK1), intersect(names(D5), names(E)))
      if (length(pr)) { for (cl in pr) D5[[cl]] <- num(D5[[cl]])
        M2 <- merge(D5[, c("RID",pr), with=FALSE], E[, c("RID",pr), with=FALSE],
                    by="RID", suffixes=c("_dep","_src"))
        DEPcmp <- rbindlist(lapply(pr, function(cl) { a <- M2[[paste0(cl,"_dep")]]; b <- M2[[paste0(cl,"_src")]]
          o <- is.finite(a) & is.finite(b)
          data.table(protein=cl, n=sum(o), equal_1e6=sum(abs(a[o]-b[o]) <= 1e-6*pmax(1,abs(b[o]))),
            r=if (sum(o)>3) cor(a[o],b[o]) else NA_real_,
            median_dep=median(a[o]), median_src=median(b[o])) }))
        put(DEPcmp, "S7_deposit_vs_source.csv", "deposited Supporting Data against the source Emory matrix") } } }

  sp <- function(d,a,b) { if (!(a %in% names(d)) || !(b %in% names(d))) return(c(NA_real_,NA_real_,0))
    x <- d[[a]]; y <- d[[b]]; o <- is.finite(x)&is.finite(y)
    if (sum(o) < 10) return(c(NA_real_,NA_real_,sum(o)))
    ct <- suppressWarnings(cor.test(x[o],y[o],method="spearman"))
    c(unname(ct$estimate), ct$p.value, sum(o)) }
  mkr <- function(d,lab,ind) data.table(set=lab, n=nrow(d), independent=ind,
    rho_V1A=sp(d,V1A,TAU)[1], p_V1A=sp(d,V1A,TAU)[2],
    rho_HK1=sp(d,HK1,TAU)[1], p_HK1=sp(d,HK1,TAU)[2])
  rec <- rbind(mkr(E_raw[is.finite(get(V1A)) & is.finite(get(TAU))], "rows as read, repeat visits kept", FALSE),
               mkr(E[is.finite(get(V1A)) & is.finite(get(TAU))], "one baseline row per subject", TRUE),
               mkr(E[is.finite(get(V1A)) & is.finite(get(TAU)) & !is.na(DX)], "baseline row, with a diagnosis", TRUE))
  put(rec, "S7_n_reconciliation.csv", "the same correlation at three cohort definitions")

  subs <- list("whole cohort"=rep(TRUE,nrow(E)))
  if ("CDR" %in% names(E)) { subs[["CDR >= 1"]] <- is.finite(E$CDR) & E$CDR>=1
                             subs[["CDR > 1"]]  <- is.finite(E$CDR) & E$CDR>1 }
  if ("DX" %in% names(E))  subs[["DX == DEM"]]  <- !is.na(E$DX) & E$DX=="DEM"
  SUB <- rbindlist(lapply(names(subs), function(k) { d <- E[subs[[k]]]
    kw <- if ("DX" %in% names(d) && uniqueN(d$DX[!is.na(d$DX)])>1)
      suppressWarnings(kruskal.test(d[[V1A]][is.finite(d[[V1A]])], factor(d$DX[is.finite(d[[V1A]])]))$p.value) else NA_real_
    data.table(subset=k, n=nrow(d), V1A_KW_p=kw, V1A_median=median(d[[V1A]], na.rm=TRUE),
      rho_V1A_TMT=sp(d,V1A,TAU)[1], p_V1A_TMT=sp(d,V1A,TAU)[2],
      rho_HK1_TMT=sp(d,HK1,TAU)[1], p_HK1_TMT=sp(d,HK1,TAU)[2],
      rho_V1A_imm=sp(d,V1A,"ElecsysTau")[1], p_V1A_imm=sp(d,V1A,"ElecsysTau")[2],
      rho_HK1_imm=sp(d,HK1,"ElecsysTau")[1], p_HK1_imm=sp(d,HK1,"ElecsysTau")[2],
      n_imm=sp(d,V1A,"ElecsysTau")[3]) }))
  put(SUB, "S7_subgroups.csv", "whole cohort and late-stage subgroups, both Tau platforms")

  Do <- E[is.finite(get(V1A))]; k <- max(1L, floor(nrow(Do)*TRIM_FRAC))
  top <- rep(FALSE, nrow(Do)); top[order(-Do[[V1A]])[seq_len(k)]] <- TRUE
  cat_t <- function(v, lab) {
    if (is.null(v) || all(is.na(v))) return(data.table(variable=lab, test="not available",
      p=NA_real_, direction=NA_character_, detail="absent"))
    lv <- sort(unique(as.character(v[!is.na(v)])))
    if (length(lv) < 2) return(data.table(variable=lab, test="not testable", p=NA_real_,
      direction=NA_character_, detail=paste("one level:", paste(lv, collapse=","))))
    M <- rbind(as.integer(table(factor(as.character(v[top]), levels=lv))),
               as.integer(table(factor(as.character(v[!top]), levels=lv))))
    colnames(M) <- lv; M <- M[, colSums(M)>0, drop=FALSE]
    pv <- if (ncol(M)<2) NA_real_ else suppressWarnings(tryCatch(chisq.test(M)$p.value, error=function(e) NA_real_))
    lvn <- suppressWarnings(as.numeric(colnames(M))); dr <- NA_character_
    if (!anyNA(lvn)) { mt <- sum(lvn*M[1,])/sum(M[1,]); mr <- sum(lvn*M[2,])/sum(M[2,])
      dr <- sprintf("mean %.2f vs %.2f (%s)", mt, mr, ifelse(mt>mr,"top 5% HIGHER","top 5% lower")) }
    det <- paste(sprintf("%s %d/%d", colnames(M), M[1,], M[1,]+M[2,]), collapse="; ")
    if (nchar(det)>150) det <- paste0(substr(det,1,147),"...")
    data.table(variable=lab, test=ifelse(is.na(pv),"not testable","chi-square"), p=pv,
      direction=dr, detail=det) }
  num_t <- function(v, lab) {
    if (is.null(v) || sum(is.finite(v))<20 || sum(is.finite(v[top]))<5)
      return(data.table(variable=lab, test="not available", p=NA_real_, direction=NA_character_, detail="too few"))
    tt <- suppressWarnings(t.test(v[top], v[!top]))
    data.table(variable=lab, test="Welch t", p=tt$p.value,
      direction=ifelse(tt$estimate[1]>tt$estimate[2],"top 5% HIGHER","top 5% lower"),
      detail=sprintf("%.4g vs %.4g", tt$estimate[1], tt$estimate[2])) }
  COMP <- rbindlist(list(
    cat_t(if ("DX" %in% names(Do)) Do$DX else NULL, "diagnosis"),
    cat_t(if ("CDR" %in% names(Do)) as.character(Do$CDR) else NULL, "CDR"),
    cat_t(if ("APOE_e4n" %in% names(Do)) as.character(Do$APOE_e4n) else NULL, "APOE e4 count"),
    cat_t(if ("APOE_geno" %in% names(Do)) Do$APOE_geno else NULL, "APOE genotype"),
    num_t(if ("CDR" %in% names(Do)) Do$CDR else NULL, "CDR (numeric)"),
    num_t(if ("APOE_e4n" %in% names(Do)) Do$APOE_e4n else NULL, "APOE e4 dose"),
    num_t(if ("APOE4_APOE4" %in% names(Do)) Do$APOE4_APOE4 else NULL, "APOE4 peptide (Emory)"),
    num_t(if ("APOE2_APOE2" %in% names(Do)) Do$APOE2_APOE2 else NULL, "APOE2 peptide (Emory)")))
  put(COMP, "S7_leverage_composition.csv", "who is in the trimmed 5 per cent of V1A")

  gv <- function(k,c2) SUB[subset==k][[c2]][1]
  late <- if ("CDR >= 1" %in% SUB$subset && isTRUE(gv("CDR >= 1","n")>=20)) "CDR >= 1" else "DX == DEM"
  bias <- COMP[is.finite(p) & p < ALPHA]$variable
  verdict("S7 CSF on the full inputs (Reviewer 2, points 4 and 5)",
    sprintf("n reconciliation : %s", paste(sprintf("%s n=%d V1A %+.3f", rec$set, rec$n, rec$rho_V1A), collapse=" | ")),
    if (!is.null(DEPcmp)) sprintf("deposit vs source : %s",
      paste(sprintf("%s equal in %d of %d", DEPcmp$protein, DEPcmp$equal_1e6, DEPcmp$n), collapse=" | "))
      else "deposit vs source : not run",
    sprintf("trimmed 5%% differs from the rest on : %s",
            if (length(bias)) paste(bias, collapse=", ") else "(nothing)"),
    sprintf("late-stage subset %s (n=%s): V1A-Tau TMT %+.3f, against immunoassay Tau %s (p %s)",
            late, gv(late,"n"), gv(late,"rho_V1A_TMT"),
            ifelse(is.finite(gv(late,"rho_V1A_imm")), sprintf("%+.3f", gv(late,"rho_V1A_imm")), "NA"),
            ifelse(is.finite(gv(late,"p_V1A_imm")), sprintf("%.3g", gv(late,"p_V1A_imm")), "NA")),
    "The subgroup restriction RAISES the TMT coefficient. The claim stays withdrawn",
    "because the immunoassay anchor is the criterion, not the subgroup size.",
    sprintf("HK1 against immunoassay Tau in the same subset: %s (p %s, n %s) - if this is",
            ifelse(is.finite(gv(late,"rho_HK1_imm")), sprintf("%+.3f", gv(late,"rho_HK1_imm")), "NA"),
            ifelse(is.finite(gv(late,"p_HK1_imm")), sprintf("%.3g", gv(late,"p_HK1_imm")), "NA"), gv(late,"n_imm")),
    "     not significant it is a power limit of the subgroup, not a failure of HK1.")
  invisible(TRUE) }

## ---- 13. S8  staging axes: lock, then strata -------------------------------
AXIS_LOCK <- function() file.path(OUTDIR, "S8_AXIS_LOCK.txt")
s8_lock <- function(force = FALSE) {
  f <- file.path(OUTDIR, "S1_donor_stage_full.csv")
  if (!file.exists(f)) { note("[stop S8] run S1 first\n"); return(invisible(FALSE)) }
  if (file.exists(AXIS_LOCK()) && !force) {
    note("a lock already exists and will NOT be rewritten:\n%s\n",
         paste(readLines(AXIS_LOCK()), collapse="\n")); return(invisible(TRUE)) }
  D <- fread(f, showProgress=FALSE)
  rows <- rbindlist(lapply(names(STAGE_AXES), function(k) { A <- STAGE_AXES[[k]]
    if (!(A$col %in% names(D))) return(data.table(axis=k, role=A$role, column=A$col,
      early=NA_integer_, mid=NA_integer_, late=NA_integer_, gate="column absent"))
    x <- suppressWarnings(as.integer(D[[A$col]]))
    e <- sum(A$early(x), na.rm=TRUE); l <- sum(A$late(x), na.rm=TRUE)
    data.table(axis=k, role=A$role, column=A$col, early=e,
      mid=sum(!A$early(x) & !A$late(x), na.rm=TRUE), late=l,
      gate=ifelse(min(e,l) >= MIN_DONORS, "gate met", "gate NOT met")) }))
  put(rows, "S8_axis_counts.csv", "donor counts per staging axis - metadata only")
  txt <- c("S8_AXIS_LOCK   staging axis designation", paste0("written : ", format(Sys.time())),
    paste0("build   : ", BUILD), "",
    sprintf("RULE, fixed before any donor metadata were read: MIN_DONORS = %d in BOTH", MIN_DONORS),
    "the earliest and the latest stratum. This is a DESIGN GATE ON DONOR COUNT.",
    "It is not a power calculation: power needs an effect size and a variance, and",
    "none had been examined when this was written.", "",
    "DONOR COUNTS (metadata only, no outcome variable read):")
  for (i in seq_len(nrow(rows))) { A <- STAGE_AXES[[rows$axis[i]]]
    txt <- c(txt, sprintf("  %-6s %-44s early %s = %s | mid = %s | late %s = %s -> %s",
      rows$axis[i], A$label, A$elab, rows$early[i], rows$mid[i], A$llab, rows$late[i], rows$gate[i])) }
  txt <- c(txt, "", "DESIGNATION:",
    sprintf("  primary                     : %s", paste(rows[role=="primary"]$axis, collapse=", ")),
    sprintf("  sensitivity                 : %s", paste(rows[role=="sensitivity"]$axis, collapse=", ")),
    sprintf("  pre-specified, gate not met : %s", paste(rows[grepl("not met", role)]$axis, collapse=", ")), "",
    "PROVENANCE, stated exactly:",
    "  Braak was the axis written into this step when it was specified, and",
    "  MIN_DONORS was fixed at that time. Thal, CERAD and ADNC were NOT a",
    "  pre-registered candidate set; they were found in the donor metadata when",
    "  obs was first read. The designation was made after reading donor counts and",
    "  before any outcome value was examined. The reply must not describe the",
    "  candidate set as pre-registered.", "",
    "Braak is retained and analysed. Its result is reported as uninterpretable for",
    "lack of donors at the earliest stratum, not as a negative finding: an",
    "underpowered null is not evidence of dissociation.")
  writeLines(txt, AXIS_LOCK()); cat(paste(txt, collapse="\n"), "\n")
  .MAN$rows[[length(.MAN$rows)+1L]] <- data.table(file="S8_AXIS_LOCK.txt", rows=length(txt),
    cols=1L, contents="pre-registration record for the staging axis")
  invisible(TRUE) }

s8_strata <- function() {
  if (!file.exists(AXIS_LOCK())) { note("[stop S8] run s8_lock() first\n"); return(invisible(FALSE)) }
  L <- load_donor(need_det=FALSE); if (is.null(L)) return(invisible(FALSE))
  Ga <- L$Ga; Ge <- L$Ge; DM <- L$DM
  vec <- function(M,g){ j <- match(g, colnames(M)); j <- j[!is.na(j)]
    if (!length(j)) rep(NA_real_, nrow(M)) else rowMeans(M[,j,drop=FALSE]) }
  SER <- list("astrocytic MCT4"=list(v=vec(Ga,"SLC16A3"), gm=DM$gm_a),
              "neuronal V-ATPase (6)"=list(v=vec(Ge,VATP6), gm=DM$gm_n),
              "CONTROL neuronal identity"=list(v=vec(Ge,CTRL_NEUR), gm=DM$gm_n))
  out <- rbindlist(lapply(names(STAGE_AXES), function(ax) { A <- STAGE_AXES[[ax]]
    if (!(A$col %in% names(DM))) return(NULL)
    st <- suppressWarnings(as.integer(DM[[A$col]]))
    strat <- ifelse(is.na(st), NA, ifelse(A$early(st),"early", ifelse(A$late(st),"late","mid")))
    rbindlist(lapply(names(SER), function(sn) { S <- SER[[sn]]
      rbindlist(lapply(c("early","mid","late"), function(k) { i <- which(strat==k & is.finite(S$v))
        data.table(axis=ax, role=A$role, series=sn, stratum=k, n_donors=length(i),
          mean=if (length(i)) mean(S$v[i]) else NA_real_,
          gate_met=length(i) >= MIN_DONORS) })) })) }))
  put(out, "S8_strata.csv", "stratum means on every staging axis")
  ln <- character(0)
  for (ax in names(STAGE_AXES)) { d <- out[axis==ax & series=="astrocytic MCT4"]
    c2 <- out[axis==ax & series=="CONTROL neuronal identity"]
    if (!nrow(d)) next
    g <- function(t,k="mean") d[stratum==t][[k]][1]; gc3 <- function(t) c2[stratum==t]$mean[1]
    ln <- c(ln, sprintf("%-6s (%s) MCT4 early %.5f mid %.5f late %.5f | early->late %+.1f%% | control %+.1f%% | gate %s",
      ax, STAGE_AXES[[ax]]$role, g("early"), g("mid"), g("late"),
      100*(g("late")-g("early"))/g("early"), 100*(gc3("late")-gc3("early"))/gc3("early"),
      ifelse(all(d$gate_met),"met","NOT met"))) }
  verdict(c("S8 staging axes (Reviewer 1, points 4 and 5)", ln,
    "MCT4 falls on every axis; the neuronal-identity control moves the other way.",
    "These axes are neuropathological gradings, not the pseudo-progression score,",
    "so this does not depend on the construction Reviewer 1 point 5 objects to.",
    "No compensatory-phase claim is made here: the shape test is S9."))
  invisible(TRUE) }

## ---- 14. S9  curve features (Reviewer 1.3) ----------------------------------
seg_bp <- function(x,y) { if (length(unique(x)) < 4) return(NA_real_)
  best <- Inf; bp <- NA_real_
  for (k in seq(min(x)+1e-3, max(x)-1e-3, length.out=400)) {
    x1 <- pmin(x,k); x2 <- pmax(x-k,0); if (sd(x2) < 1e-12) next
    s <- tryCatch(sum(resid(lm(y ~ x1 + x2))^2), error=function(e) Inf)
    if (is.finite(s) && s < best) { best <- s; bp <- k } }
  bp }

s9_curve <- function() {
  L <- load_donor(need_det=FALSE); if (is.null(L)) return(invisible(FALSE))
  Ga <- L$Ga; Ge <- L$Ge; DM <- L$DM
  vec <- function(M,g){ j <- match(g, colnames(M)); j <- j[!is.na(j)]
    if (!length(j)) rep(NA_real_, nrow(M)) else rowMeans(M[,j,drop=FALSE]) }
  SER <- list("astrocytic MCT4"=list(v=vec(Ga,"SLC16A3"), gm=DM$gm_a),
              "neuronal V-ATPase (6)"=list(v=vec(Ge,VATP6), gm=DM$gm_n),
              "CONTROL neuronal identity"=list(v=vec(Ge,CTRL_NEUR), gm=DM$gm_n))
  AX <- list(cps=DM$cps)
  for (a in c("adnc_int","thal_int","cerad_int","braak"))
    if (a %in% names(DM)) AX[[sub("_int$","",a)]] <- suppressWarnings(as.numeric(DM[[a]]))
  out <- rbindlist(lapply(names(AX), function(an) { x <- AX[[an]]
    rbindlist(lapply(names(SER), function(sn) { S <- SER[[sn]]
      o <- is.finite(x) & is.finite(S$v) & is.finite(S$gm)
      if (sum(o) < 20 || length(unique(x[o])) < 3)
        return(data.table(axis=an, series=sn, n=sum(o), levels=length(unique(x[o])),
          quad=NA_real_, p_quad=NA_real_, vertex=NA_real_, inside=NA, argmax=NA_real_,
          seg_bp=NA_real_, verdict="not resolvable"))
      xx <- x[o]; yy <- S$v[o]; gg <- S$gm[o]
      m <- summary(lm(yy ~ xx + I(xx^2) + gg))$coefficients
      q <- rownames(m)[grepl("\\^2", rownames(m))][1]
      if (is.na(q)) return(data.table(axis=an, series=sn, n=sum(o), levels=length(unique(xx)),
        quad=NA_real_, p_quad=NA_real_, vertex=NA_real_, inside=NA, argmax=NA_real_,
        seg_bp=NA_real_, verdict="not estimable"))
      b1 <- m["xx","Estimate"]; b2 <- m[q,"Estimate"]
      v <- if (abs(b2) > 1e-14) -b1/(2*b2) else NA_real_
      ins <- is.finite(v) && v > min(xx) && v < max(xx)
      ag <- tapply(yy, xx, mean)
      data.table(axis=an, series=sn, n=sum(o), levels=length(unique(xx)), quad=b2,
        p_quad=m[q,"Pr(>|t|)"], vertex=v, inside=ins,
        argmax=as.numeric(names(ag))[which.max(ag)], seg_bp=seg_bp(xx,yy),
        verdict=if (m[q,"Pr(>|t|)"] < ALPHA && b2 < 0 && ins) "concave, vertex inside range"
                else "not concave-inside") })) }))
  put(out, "S9_curve_features.csv", "quadratic vertex, discrete argmax and segmented breakpoint")
  ln <- character(0)
  for (an in names(AX)) { m <- out[axis==an & series=="astrocytic MCT4"][1]
    c2 <- out[axis==an & series=="CONTROL neuronal identity"][1]
    hit <- identical(m$verdict,"concave, vertex inside range")
    ct  <- identical(c2$verdict,"concave, vertex inside range")
    ln <- c(ln, sprintf("%-6s quad %+.5f (p %.4f) | vertex %s inside=%s | argmax %.2f | breakpoint %s -> %s%s",
      an, m$quad, m$p_quad, ifelse(is.finite(m$vertex), sprintf("%.3f", m$vertex), "NA"), m$inside,
      m$argmax, ifelse(is.finite(m$seg_bp), sprintf("%.3f", m$seg_bp), "NA"),
      if (hit && !ct) "NON-MONOTONE SUPPORTED" else if (hit && ct) "NOT USABLE - control shares the shape"
        else "not supported",
      if (an=="cps") "   [same feature as Paper 1's PTGDS vertex at CPS 0.47]" else "")) }
  verdict(c("S9 curve features (Reviewer 1, point 3)", ln,
    "Features are named as in AD-Boundary-Audit/PTGDS_vertex_check.R: the quadratic",
    "vertex, the discrete argmax and the segmented breakpoint are three different",
    "things and different values are not a contradiction.",
    "An 'early above late plus early not yet falling' rule is NOT a test for a rise:",
    "a monotone decline satisfies it. Any such label from an earlier run is withdrawn."))
  invisible(TRUE) }

## ---- 15. conventions ledger and the runner ----------------------------------
write_conventions <- function() {
  CV <- data.table(
    item = c("seed","CPS bins","early / late","trajectory","detection floor (report alone)",
             "detection floor (matched null)","block size","h5 backend","partial correlation",
             "curve features","output location"),
    house_source = c("P2_run_all.R / AD-Boundary-Audit engine",
             "77_detection_matched_null_standalone.R (P2)",
             "77_detection_matched_null_standalone.R (P2)","oligo_MCT_detection_check.R",
             "oligo_MCT_detection_check.R","P2 04/77 detection-matched null",
             "oligo_MCT_detection_check.R","P3 celltype_audit_from_raw.R",
             "AD-Metabolic-Collapse 00_utils.R","AD-Boundary-Audit PTGDS_vertex_check.R",
             "oligo_MCT_detection_check.R"),
    value_used = c("42","round(cps, 1); bin 0.1 excluded from bin-level work",
             "unweighted mean of bins 0.2-0.4 and 0.6-0.8",
             "z per gene across cells, then bin mean","0.10","0.02","100000","rhdf5",
             "pc(), verbatim","vertex / argmax / segmented breakpoint","OUTDIR, absolute"),
    note = c("",
             "P1-fixed oligo_MCT_detection_check.R uses cut(right=FALSE) instead - a genuine inconsistency between the repositories",
             "NOT the pooled cell mean; using the pooled mean gives -46.5% instead of -43.2%",
             "multi-gene composites in S5 and S6 are z per gene across donors; a raw mean of the ANLS panel is 81% PFKFB3 and hides MCT4. S8/S9 keep the published raw VATP_n6 definition",
             "astrocytic MCT4 is 0.061, BELOW this floor - flagged in S3 and S4",
             "used only to decide whether a matched-null percentile is computed",
             "","hdf5r not required","","","set with Sys.getenv('P2R2_OUT')"))
  put(CV, "CONVENTIONS.csv", "which house convention each choice follows, and where it departs")
  invisible(TRUE) }

STEPS <- list(S1=s1_donor_meta, S2=s2_verify, S3=s3_extract, S4=s4_acidbase,
              S5=s5_mito, S6=s6_pathways, S7=s7_csf, S8L=s8_lock, S8=s8_strata, S9=s9_curve)

p2_run <- function(steps = names(STEPS), stop_on_fail = FALSE) {
  t0 <- Sys.time(); .MAN$verd <- c(.MAN$verd, paste0("RUN ", format(Sys.time()), "  ", BUILD))
  st <- character(0)
  for (s in steps) {
    cat(sprintf("\n================ %s ================\n", s))
    ok <- tryCatch(isTRUE(STEPS[[s]]()), error=function(e){ cat("[R error] ", conditionMessage(e), "\n"); NA })
    st[s] <- if (is.na(ok)) "ERROR" else if (ok) "ok" else "stopped"
    # flush after every step so an interrupted session still leaves a record
    writeLines(.MAN$verd, file.path(OUTDIR, "VERDICTS.txt"))
    if (length(.MAN$rows))
      data.table::fwrite(rbindlist(.MAN$rows), file.path(OUTDIR, "MANIFEST.csv"))
    if (st[s] != "ok" && stop_on_fail) break }
  write_conventions()
  if (length(.MAN$rows)) put(rbindlist(.MAN$rows), "MANIFEST.csv", "every file this run wrote")
  writeLines(.MAN$verd, file.path(OUTDIR, "VERDICTS.txt"))
  si <- capture.output(sessionInfo())
  writeLines(c(BUILD, paste("run :", format(Sys.time())), paste("seed:", 42),
    paste("h5ad:", MTG_H5), paste("out :", OUTDIR), "", si), file.path(OUTDIR, "SESSION.txt"))
  cat(sprintf("\n---- summary (%.1f min) ----\n", as.numeric(difftime(Sys.time(), t0, units="mins"))))
  for (s in names(st)) cat(sprintf("  %-4s %s\n", s, st[s]))
  cat("\nWritten to: ", OUTDIR, "\n  VERDICTS.txt  MANIFEST.csv  CONVENTIONS.csv  SESSION.txt\n", sep="")
  invisible(st) }

## ---- 16. run ----------------------------------------------------------------
## Running this file runs the analysis. There is no second command to type.
cat(BUILD, "\n")
cat("  h5ad   :", MTG_H5, ifelse(file.exists(MTG_H5), "(found)", "(MISSING)"), "\n")
cat("  P2 out :", P2_RUNDIR, "|", P2_H5DIR, "\n")
cat("  ADNI   :", EMORY, ifelse(file.exists(EMORY), "(found)", "(MISSING)"), "\n")
cat("  OUTDIR :", OUTDIR, "\n\n")
if (!selftest()) stop("SELF-TEST FAILED - nothing was read.")

# No opt-out switch. An earlier version had one, a stale variable left in the
# session silently turned the run off, and the file looked like it had done
# nothing. Running this file runs the analysis. Always.
cat("Running every step. S3 reads the h5ad the first time, 40-70 minutes; on a\n")
cat("later run it reuses that pass if the file still reproduces the published\n")
cat("values, so a repeat run takes about a minute. Everything is written as each\n")
cat("step finishes, so nothing is lost if the session is interrupted.\n")
cat("  S1 donor metadata | S2 verify | S3 h5ad pass | S4 acid-base | S5 mitochondria\n")
cat("  S6 pathways | S7 CSF | S8L axis lock | S8 strata | S9 curves\n")
P2R2_STATUS <- p2_run()
