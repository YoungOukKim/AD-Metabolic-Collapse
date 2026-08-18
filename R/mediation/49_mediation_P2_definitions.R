# =============================================================================
# 49_mediation_P2_definitions.R
#
# WHY THIS EXISTS. I made two errors and I am correcting both.
#
#   ERROR 1 (withdrawn): I claimed P2's Methods contradicted its code on Bin 0.1.
#     It does not. Methods excludes Bin 0.1 from BIN-LEVEL analyses; 02_global's
#     CPS_MIN = 0.1 is a CELL filter for the DONOR-LEVEL analysis. Different
#     layers, mutually consistent. I misread it.
#
#   ERROR 2 (withdrawn): I claimed the V-ATPase composite was inconsistent
#     (6 subunits in 01_extract, 10 in 02_global). It is not: P2 deliberately
#     uses TEN subunits for ASTROCYTES (vatpase_genes) and SIX for NEURONS
#     (vatpase_n). That is a considered choice, not a bug.
#
#   BUT ERROR 2 CONTAMINATES MY OWN WORK. In script 42 I built the neuronal
#   V-ATPase composite from TEN subunits. So my headline mediation number --
#   "MCT4 mediates 106% of the CPS effect on neuronal V-ATPase" -- was computed
#   on a composite P2 does not use. If that number goes into the revision as-is,
#   the manuscript would contain a figure the repository cannot reproduce.
#
#   This script recomputes EVERYTHING with P2's own definitions, verbatim from
#   R/data_extraction/01_extract_seaad.R:
#
#     anls_genes    <- c("SLC2A1","LDHA","SLC16A1")                        # 3
#     vatpase_genes <- 10 subunits                                          # astro
#     vatpase_n     <- c("ATP6V1A","ATP6V1B2","ATP6V0A1",
#                        "ATP6V0C","ATP6V0D1","ATP6V1E1")                   # 6, neuron
#     ed_genes      <- c("SLC2A1","LDHA","SLC16A1","PKM","HK1")             # 5
#     Excitatory_Neuron := Subclass matching ^L[0-9]|IT$|ET$|CT$|NP$|L6b
#     Astrocyte        := Subclass ^Astro
#
# WHAT IT PRODUCES, ready to paste into the revision:
#   1. G1  P2's published numbers, reproduced from raw. If these do not match,
#          nothing below is trustworthy and the script stops.
#   2. Mediation, BOTH directions, on P2's composites. -> Table 3
#   3. The ANLS composite donor-level slope. P2 reports beta = -0.295, p = 0.016.
#      My earlier donor-level recomputation gave p = 0.132. Settle it from raw.
#   4. SENSITIVITY: the donor-level coupling (partial r = +0.466) with and
#      without the early donors (CPS < 0.2). This is a real robustness limit and
#      the revision must disclose it. It is NOT a Methods/code contradiction --
#      I withdraw that claim.
#   5. Detection-matched null for MCT4, and the ambient control
#      (astrocytic vs microglial MCT4), both as deposited numbers.
#
# set.seed(42).  Run with the working directory set to the repository root.
# =============================================================================

# Paths are RELATIVE to the repository root. Set the working directory there,
# or override with the environment variables SEAAD_H5AD / ROSMAP_ASTRO /
# ROSMAP_CLIN / P2_OUT_DIR. Raw data are not redistributable; see README.md.

suppressPackageStartupMessages({ library(rhdf5); library(Matrix); library(data.table) })
set.seed(42)

.chk <- c("mean","median","mad","sd","sum","rank","cor","lm","apply","sapply","tapply")
.bad <- Filter(function(f) { v <- get0(f, inherits=TRUE); !is.function(v) }, .chk)
if (length(.bad)) stop(sprintf("Masked base function(s): %s. rm(list=c(%s)) or restart R.",
                               paste(.bad, collapse=", "), paste(sprintf('"%s"', .bad), collapse=",")))

SEA <- Sys.getenv("SEAAD_H5AD", unset = "data/raw/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad")
OUT <- Sys.getenv("P2_OUT_DIR", unset = "output/tables"); dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
BS <- 5000L
BINS <- round(seq(0.2,0.9,0.1),1); KEY <- sprintf("%.1f", BINS); NB <- 8L
I_E <- 1:3; I_L <- 5:7

## ---- P2's definitions. Verbatim. Do not "improve" them. --------------------
ANLS_GENES    <- c("SLC2A1","LDHA","SLC16A1")
VATPASE_ASTRO <- c("ATP6V1A","ATP6V1B2","ATP6V0D1","ATP6V0A1",
                   "ATP6V1C1","ATP6V1E1","ATP6V1H","ATP6V0C","ATP6V0E1","ATP6V0B")   # 10
VATPASE_NEUR  <- c("ATP6V1A","ATP6V1B2","ATP6V0A1","ATP6V0C","ATP6V0D1","ATP6V1E1")  # 6
ED_GENES      <- c("SLC2A1","LDHA","SLC16A1","PKM","HK1")
NEUR_PAT      <- "^L[0-9]|IT$|ET$|CT$|NP$|L6b"
ASTR_PAT      <- "^Astro"

## ---- extra outcomes (mine, clearly labelled as not-P2) ---------------------
NAK_ASTRO  <- c("ATP1A2","ATP1B2","ATP1A1","ATP1B1","FXYD1","FXYD6")
GLYC_ASTRO <- c("HK1","HK2","PFKP","GAPDH","PKM","LDHA")
ACT_NEUR   <- c("FOS","NPAS4","ARC","EGR1","JUNB")
GLUT_ASTRO <- c("SLC1A2","SLC1A3","GLUL")
MITO_ASTRO <- c("NDUFS1","NDUFV1","NDUFA1","UQCRC1","COX4I1","ATP5F1A","ATP5F1B")

sink(file.path(OUT,"49_mediation_P2.txt"), split=TRUE)
cat("=== 49: mediation, recomputed with P2's OWN composite definitions ===\n\n")
cat("  neuronal V-ATPase = 6 subunits (P2's vatpase_n), NOT the 10 I used in 42.\n")
cat("  Every number below is reproducible from the deposited repository code.\n\n")

## ======================= reader =============================================
rd <- function(k){
  v <- tryCatch({ ca <- as.vector(h5read(SEA,paste0("obs/__categories/",k)))
                  co <- as.integer(h5read(SEA,paste0("obs/",k))); co[co<0] <- NA; ca[co+1L] },
                error=function(e) NULL)
  if (!is.null(v)) return(v); tryCatch(as.vector(h5read(SEA,paste0("obs/",k))), error=function(e) NULL) }
parse_braak <- function(x){ if (is.numeric(x)) return(as.integer(round(x)))
  s <- trimws(gsub("BRAAK","",toupper(trimws(as.character(x)))))
  r <- c("0"=0L,"I"=1L,"II"=2L,"III"=3L,"IV"=4L,"V"=5L,"VI"=6L); o <- r[s]
  n <- suppressWarnings(as.integer(s)); o[is.na(o)&!is.na(n)] <- n[is.na(o)&!is.na(n)]; as.integer(o) }

genes <- as.character(h5read(SEA,"var/_index")); nG <- length(genes)
ip  <- as.numeric(h5read(SEA,"X/indptr", bit64conversion="double")); nC <- length(ip)-1L
cps <- suppressWarnings(as.numeric(rd("Continuous Pseudo-progression Score")))
sc  <- as.character(rd("Subclass")); dn <- as.character(rd("Donor ID"))
brk <- parse_braak(rd("Braak"))
isA <- grepl(ASTR_PAT, sc, ignore.case=TRUE)
isN <- grepl(NEUR_PAT, sc)
isM <- grepl("^Micro", sc, ignore.case=TRUE)          # ambient control
cat(sprintf("[1] cells: astro %d | exc.neuron %d | microglia %d\n", sum(isA), sum(isN), sum(isM)))

bi  <- match(sprintf("%.1f", round(cps,1)), KEY)
don <- sort(unique(dn[is.finite(cps)])); nD <- length(don); di <- match(dn, don)
keep<- (isA | isN | isM) & is.finite(cps) & !is.na(di)

# donor x gene, per cell type. bin x gene for astro (to reproduce P2's bin stats).
mk <- function(n1,n2) Matrix(0,n1,n2,sparse=TRUE)
SdA<-mk(nD,nG); NdA<-mk(nD,nG); CdA<-numeric(nD)
SdN<-mk(nD,nG); CdN<-numeric(nD)
SdM<-mk(nD,nG); CdM<-numeric(nD)
SbA<-mk(NB,nG); NbA<-mk(NB,nG); CbA<-numeric(NB)
SbM<-mk(NB,nG); CbM<-numeric(NB)
t0 <- Sys.time()
for (s0 in seq(1L,nC,by=BS)) {
  e0 <- min(s0+BS-1L,nC); gi <- (s0:e0)[keep[s0:e0]]; if (!length(gi)) next
  sp <- ip[s0]; cnt <- ip[e0+1L]-sp; if (cnt<=0) next
  ci <- h5read(SEA,"X/indices",start=sp+1,count=cnt,bit64conversion="double")
  cd <- as.numeric(h5read(SEA,"X/data",start=sp+1,count=cnt))
  np <- as.integer(diff(ip[s0:(e0+1L)])); cof <- rep(s0:e0, times=np); m <- keep[cof]
  if (!any(m)) { rm(ci,cd); next }
  loc <- match(cof[m], gi)
  M  <- sparseMatrix(i=loc, j=as.integer(ci[m])+1L, x=cd[m], dims=c(length(gi),nG))
  Mn <- M; Mn@x <- rep(1,length(M@x))
  aa <- isA[gi]; nn <- isN[gi]; mm <- isM[gi]
  if (any(aa)) {
    f <- factor(di[gi][aa], levels=1:nD); D <- fac2sparse(f, drop.unused.levels=FALSE)
    SdA <- SdA + D %*% M[aa,,drop=FALSE]; NdA <- NdA + D %*% Mn[aa,,drop=FALSE]
    CdA <- CdA + as.numeric(rowSums(D))
    ib <- bi[gi][aa]; ok <- !is.na(ib)
    if (any(ok)) { fb <- factor(ib[ok], levels=1:NB); Db <- fac2sparse(fb, drop.unused.levels=FALSE)
      SbA <- SbA + Db %*% M[aa,,drop=FALSE][ok,,drop=FALSE]
      NbA <- NbA + Db %*% Mn[aa,,drop=FALSE][ok,,drop=FALSE]
      CbA <- CbA + as.numeric(rowSums(Db)) } }
  if (any(nn)) {
    f <- factor(di[gi][nn], levels=1:nD); D <- fac2sparse(f, drop.unused.levels=FALSE)
    SdN <- SdN + D %*% M[nn,,drop=FALSE]; CdN <- CdN + as.numeric(rowSums(D)) }
  if (any(mm)) {
    f <- factor(di[gi][mm], levels=1:nD); D <- fac2sparse(f, drop.unused.levels=FALSE)
    SdM <- SdM + D %*% M[mm,,drop=FALSE]; CdM <- CdM + as.numeric(rowSums(D))
    ib <- bi[gi][mm]; ok <- !is.na(ib)
    if (any(ok)) { fb <- factor(ib[ok], levels=1:NB); Db <- fac2sparse(fb, drop.unused.levels=FALSE)
      SbM <- SbM + Db %*% M[mm,,drop=FALSE][ok,,drop=FALSE]
      CbM <- CbM + as.numeric(rowSums(Db)) } }
  rm(M,Mn,ci,cd); cat(sprintf("    %d / %d\r", e0, nC))
}
cat(sprintf("\n    done (%.1f min)\n", as.numeric(difftime(Sys.time(),t0,units="mins")))); gc()

A  <- as.matrix(SdA)/pmax(CdA,1); colnames(A) <- genes    # astro donor mean
DA <- as.matrix(NdA)/pmax(CdA,1); colnames(DA) <- genes    # astro detection
N  <- as.matrix(SdN)/pmax(CdN,1); colnames(N) <- genes    # neuron donor mean
MG <- as.matrix(SdM)/pmax(CdM,1); colnames(MG) <- genes    # microglia donor mean
BA <- as.matrix(SbA)/pmax(CbA,1); colnames(BA) <- genes    # astro bin mean
BM <- as.matrix(SbM)/pmax(CbM,1); colnames(BM) <- genes    # microglia bin mean
dcps <- as.numeric(tapply(cps[keep], factor(dn[keep], levels=don), function(z) base::mean(z, na.rm=TRUE)))
dbrk <- as.integer(tapply(brk[keep], factor(dn[keep], levels=don), function(z) z[1]))
u <- which(CdA >= 20 & CdN >= 20 & is.finite(dcps))
cat(sprintf("  donors with >=20 astro AND >=20 exc.neurons: %d\n", length(u)))

cmp <- function(M_, gs) { g <- intersect(gs, colnames(M_)); rowMeans(M_[, g, drop=FALSE], na.rm=TRUE) }
pctv <- function(v){ e<-base::mean(v[I_E]); l<-base::mean(v[I_L])
                     if(!is.finite(e)||abs(e)<1e-9) NA_real_ else (l-e)/e*100 }

## ======================= G1: reproduce P2 ===================================
cat("\n[G1] reproducing P2's published bin-level %changes from raw\n")
G1T <- data.table(gene=c("SLC16A3","HK2","LDHA","PDK1","SLC16A1","SLC2A1"),
                  published=c(-43.2,-35.2,-20.8,-16.5,-11.1,-7.4))
G1T[, recomputed := sapply(gene, function(g) round(pctv(BA[,g]),1))]
G1T[, delta := round(recomputed - published, 1)]
G1T[, ok := abs(delta) <= 1]
print(G1T)
G1 <- isTRUE(all(G1T$ok))
cat(sprintf("  G1 >>> %s\n", if (G1) "PASS. The reader matches P2. Everything below is on the same footing."
            else "FAIL. Reader does not reproduce P2. STOP -- do not use these numbers."))
if (!G1) { sink(); stop("G1 FAILED") }

## ======================= [2] P2 composites, donor level =====================
cat("\n[2] P2's composites at the donor level\n")
d <- data.table(
  donor    = don[u],
  cps      = dcps[u],
  braak    = as.numeric(dbrk[u]),
  MCT4     = A[u,"SLC16A3"],
  ANLS     = cmp(A,ANLS_GENES)[u],                 # P2: SLC2A1 + LDHA + SLC16A1
  VATP_a10 = cmp(A,VATPASE_ASTRO)[u],              # P2: 10 subunits, astro
  VATP_n6  = cmp(N,VATPASE_NEUR)[u],               # P2: 6 subunits, NEURON  <-- the fix
  ED_energy= cmp(A,ED_GENES)[u],
  LDHB_n   = N[u,"LDHB"],
  LAMP1_n  = N[u,"LAMP1"],
  # my additional outcomes -- flagged as NOT part of P2's published composites
  NaK_a    = cmp(A,NAK_ASTRO)[u],
  GLYC_a   = cmp(A,GLYC_ASTRO)[u],
  MITO_a   = cmp(A,MITO_ASTRO)[u],
  GLUT_a   = cmp(A,GLUT_ASTRO)[u],
  ACT_n    = cmp(N,ACT_NEUR)[u])
d[, VATP_n10 := cmp(N,VATPASE_ASTRO)[u]]           # the 10-subunit version I wrongly used
fwrite(d, file.path(OUT,"49_donor_level.csv"))
cat(sprintf("  donors: %d\n", nrow(d)))

fit1 <- function(y, x) { o <- is.finite(y)&is.finite(x); f <- lm(y[o] ~ x[o])
  co <- summary(f)$coefficients; c(beta=co[2,1], p=co[2,4], n=sum(o)) }
S <- rbindlist(lapply(c("MCT4","ANLS","VATP_a10","VATP_n6","VATP_n10","ED_energy",
                        "NaK_a","GLYC_a","MITO_a","ACT_n"), function(v){
  a <- fit1(d[[v]], d$cps); b <- fit1(d[[v]], d$braak)
  data.table(variable=v, beta_CPS=round(a[1],4), p_CPS=signif(a[2],3),
             beta_Braak=round(b[1],4), p_Braak=signif(b[2],3)) }))
print(S); fwrite(S, file.path(OUT,"49_donor_slopes.csv"))

cat("\n  --- P2 reports ANLS composite beta = -0.295, p = 0.016 (bin-level, n = 8) ---\n")
anls_bin <- BA[, intersect(ANLS_GENES, colnames(BA)), drop=FALSE]
anls_bv  <- rowMeans(anls_bin); anls_nb <- anls_bv / anls_bv[1]
fb <- lm(anls_nb ~ BINS); cb <- summary(fb)$coefficients
cat(sprintf("  bin-level  (n=8, normalized to Bin 0.2): beta = %+.3f  p = %.4f\n", cb[2,1], cb[2,4]))
cat(sprintf("  donor-level(n=%d)                      : beta = %+.4f  p = %.4f\n",
            nrow(d), S[variable=="ANLS"]$beta_CPS, S[variable=="ANLS"]$p_CPS))
cat("  If these disagree, the manuscript must say WHICH ONE it is reporting.\n")

## ======================= [3] MEDIATION, P2 definitions ======================
cat("\n[3] MEDIATION -- both directions, P2's composites\n")
med <- function(y, m, x, ylab, mlab, xlab){
  o <- is.finite(y)&is.finite(m)&is.finite(x); if (sum(o)<20) return(NULL)
  f0 <- lm(y[o] ~ x[o]); f1 <- lm(y[o] ~ x[o] + m[o])
  c0 <- summary(f0)$coefficients; c1 <- summary(f1)$coefficients
  b0 <- c0[2,1]; b1 <- c1[2,1]
  data.table(axis=xlab, outcome=ylab, mediator=mlab, n=sum(o),
             beta=round(b0,4), p=signif(c0[2,4],3),
             beta_adj=round(b1,4), p_adj=signif(c1[2,4],3),
             pct_mediated=round(100*(b0-b1)/b0,1)) }
MED <- rbindlist(list(
  # FORWARD: MCT4 as the mediator.  ** VATP_n6 is P2's own neuronal composite **
  med(d$VATP_n6,  d$MCT4, d$cps, "Neuronal V-ATPase (6 subunits, P2)", "Astrocytic MCT4","CPS"),
  med(d$LAMP1_n,  d$MCT4, d$cps, "Neuronal LAMP1",                     "Astrocytic MCT4","CPS"),
  med(d$LDHB_n,   d$MCT4, d$cps, "Neuronal LDHB",                      "Astrocytic MCT4","CPS"),
  med(d$GLYC_a,   d$MCT4, d$cps, "Astrocytic glycolytic program",      "Astrocytic MCT4","CPS"),
  med(d$NaK_a,    d$MCT4, d$cps, "Astrocytic Na+/K+-ATPase",           "Astrocytic MCT4","CPS"),
  med(d$MITO_a,   d$MCT4, d$cps, "Astrocytic mitochondrial ETC",       "Astrocytic MCT4","CPS"),
  med(d$ACT_n,    d$MCT4, d$cps, "Neuronal activity genes",            "Astrocytic MCT4","CPS"),
  # REVERSE
  med(d$MCT4, d$VATP_n6, d$cps, "Astrocytic MCT4","Neuronal V-ATPase (6, P2)","CPS"),
  med(d$MCT4, d$ACT_n,   d$cps, "Astrocytic MCT4","Neuronal activity genes",  "CPS"),
  med(d$MCT4, d$NaK_a,   d$cps, "Astrocytic MCT4","Astrocytic Na+/K+-ATPase", "CPS"),
  med(d$MCT4, d$GLUT_a,  d$cps, "Astrocytic MCT4","Astrocytic glutamate uptake","CPS"),
  med(d$MCT4, d$MITO_a,  d$cps, "Astrocytic MCT4","Astrocytic mito ETC",      "CPS"),
  # and on Braak, the shared neuropathology axis
  med(d$VATP_n6, d$MCT4, d$braak, "Neuronal V-ATPase (6, P2)","Astrocytic MCT4","Braak"),
  med(d$MCT4, d$VATP_n6, d$braak, "Astrocytic MCT4","Neuronal V-ATPase (6, P2)","Braak")), fill=TRUE)
print(MED); fwrite(MED, file.path(OUT,"49_mediation_TABLE3.csv"))

cat("\n  --- the number that changes: 6-subunit vs the 10 I used in 42 ---\n")
m6  <- med(d$VATP_n6,  d$MCT4, d$cps, "V-ATPase 6 (P2)",  "MCT4","CPS")
m10 <- med(d$VATP_n10, d$MCT4, d$cps, "V-ATPase 10 (mine)","MCT4","CPS")
print(rbind(m6, m10))
cat("  42 reported 106.1% using the 10-subunit composite. The 6-subunit value is\n")
cat("  what belongs in the revision, because that is what P2's code computes.\n")

## ======================= [4] SENSITIVITY: the early donors ==================
cat("\n[4] SENSITIVITY -- coupling with and without the early donors\n")
pcor <- function(x,y,z){ o <- is.finite(x)&is.finite(y)&is.finite(z)
  rx <- residuals(lm(x[o]~z[o])); ry <- residuals(lm(y[o]~z[o]))
  ct <- cor.test(rx,ry); c(r=unname(ct$estimate), p=ct$p.value, n=sum(o)) }
gm_a <- rowMeans(A[u,,drop=FALSE]); gm_n <- rowMeans(N[u,,drop=FALSE])
SENS <- rbindlist(lapply(list(
  list(rep(TRUE,nrow(d)), "all donors"),
  list(d$cps >= 0.2,      "CPS >= 0.2 (early donors dropped)"),
  list(d$cps >= 0.3,      "CPS >= 0.3")), function(z){
  k <- z[[1]]; if (sum(k) < 20) return(NULL)
  a <- pcor(d$MCT4[k], d$VATP_n6[k], d$cps[k])
  b <- pcor(d$MCT4[k], d$VATP_n6[k], d$cps[k])   # gmean-adjusted below
  o <- is.finite(d$MCT4)&k
  rx <- residuals(lm(d$MCT4[o] ~ d$cps[o] + gm_a[o]))
  ry <- residuals(lm(d$VATP_n6[o] ~ d$cps[o] + gm_n[o]))
  ct <- cor.test(rx,ry)
  data.table(subset=z[[2]], n=a[3],
             partial_r_CPSadj=round(a[1],3), p_CPSadj=signif(a[2],3),
             partial_r_CPS_gmean=round(unname(ct$estimate),3),
             p_CPS_gmean=signif(ct$p.value,3)) }), fill=TRUE)
print(SENS); fwrite(SENS, file.path(OUT,"49_sensitivity_coupling.csv"))
cat("\n  P2 reports partial r = +0.466. If the CPS>=0.2 row loses significance, the\n")
cat("  revision must SAY SO. It is a robustness limit, not a contradiction --\n")
cat("  I earlier called it a Methods/code contradiction and that was WRONG.\n")

## ======================= [5] detection-matched null + ambient ===============
cat("\n[5] detection-matched null and the ambient control\n")
det <- colMeans(DA); b02 <- BA[1,]
ok  <- which(det >= 0.02 & b02 >= 0.01 & is.finite(b02))
pcta<- apply(BA[,ok,drop=FALSE], 2, pctv)
dv  <- det[ok]; d4 <- det[match("SLC16A3", genes)]
band<- which(dv >= d4*0.75 & dv <= d4*1.25)
m4p <- pcta[match("SLC16A3", genes[ok])]
pctl<- 100*base::mean(pcta[band] < m4p, na.rm=TRUE)
cat(sprintf("  MCT4 detection %.4f | matched band %.3f-%.3f | n = %d genes\n",
            d4, d4*0.75, d4*1.25, length(band)))
cat(sprintf("  MCT4 %%change %.1f%%  ->  percentile %.1f within the band\n", m4p, pctl))
NULLTAB <- data.table(gene="SLC16A3", detection=round(d4,4), n_matched=length(band),
                      pct_change=round(m4p,1), percentile=round(pctl,2),
                      band_median=round(stats::median(pcta[band],na.rm=TRUE),1))
print(NULLTAB); fwrite(NULLTAB, file.path(OUT,"49_detection_null.csv"))

AMB <- data.table(gene=c("SLC16A3","HK2","AQP4","PTGDS"))
AMB[, astro_pct := sapply(gene, function(g) round(pctv(BA[,g]),1))]
AMB[, micro_pct := sapply(gene, function(g) if (g %in% colnames(BM)) round(pctv(BM[,g]),1) else NA_real_)]
AMB[, ratio := round(astro_pct/pmax(abs(micro_pct),1e-6),2)]
cat("\n  ambient control (a bleed-through artefact would move BOTH cell types alike):\n")
print(AMB); fwrite(AMB, file.path(OUT,"49_ambient_control.csv"))

cat("\n\n=== 49 done. Send back:\n")
cat("    49_mediation_P2.txt, 49_mediation_TABLE3.csv, 49_donor_slopes.csv,\n")
cat("    49_sensitivity_coupling.csv, 49_detection_null.csv, 49_ambient_control.csv,\n")
cat("    49_donor_level.csv ===\n")
sink()
