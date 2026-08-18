# =============================================================================
# 47_braak_common_axis.R
#
# THE MISTAKE I MADE, AND THE CEO CAUGHT IT.
#
#   I compared SEA-AD's CPS to ROSMAP's pseudotime as if they were the same ruler.
#   They are not, and the cohorts are not the same kind of study.
#
#     SEA-AD  : an AD PROGRESSION cohort. CPS is built to span the disease
#               trajectory. It is a tool for finding the INFLECTION.
#     ROSMAP  : Religious Orders Study + Memory and Aging Project. A COMMUNITY
#               AGEING cohort followed to autopsy. Most participants are normal
#               or mildly impaired. pseudotime median 0.160, 3rd quartile 0.280.
#               It is an END-STATE SNAPSHOT of pathological burden, not a
#               progression dynamic.
#
#   On top of that: SEA-AD is MTG (intermediate vulnerability), ROSMAP is DLPFC
#   (late). Two different rulers AND two different regions.
#
#   In 45 I ran a "kill test" across those two rulers and then declared the
#   hypothesis DEAD. That verdict was not earned. I am withdrawing it until the
#   axes are matched.
#
# THE ONLY HONEST COMMON AXIS IS NEUROPATHOLOGY.
#
#   SEA-AD obs has Braak. ROSMAP clinical has braaksc. Both are the same
#   standard staging. Also CERAD (SEA-AD "CERAD score" / ROSMAP ceradsc) and
#   cognition (SEA-AD "Cognitive Status" / ROSMAP cogdx).
#
#   Align on Braak and the question finally becomes answerable:
#     same Braak, MTG MCT4 != DLPFC MCT4   ->  it is a REGION effect
#     same Braak, MTG MCT4 == DLPFC MCT4   ->  it was the RULER all along
#
# JOBS
#   A  Braak/CERAD/cognition as the shared axis. MCT4 slope on Braak, both
#      cohorts, donor-level. Same direction or not.
#   B  LEVERAGE DIAGNOSIS. Every reactive marker (GFAP, SERPINA3, CD44, OSMR, C3)
#      crashes 50-90% at BIN 0.3 and again at BIN 0.6, in lockstep. That is not
#      biology. Find out which donors own those bins. If they are leverage, the
#      "GFAP peaks at bin 0.5" observation dies with them -- and I will say so.
#   C  THE CEO'S HYPOTHESIS, TESTED DIRECTLY. Donor-level:
#         does a MORE REACTIVE donor have LESS MCT4?
#      Controlled for Braak. Both cohorts. Mediation in BOTH directions.
#      SEA-AD subtypes say reactive Astro_3 has LOW MCT4 (0.0354).
#      ROSMAP states say reactive-adjacent Ast.5 has the HIGHEST MCT4 (0.0707).
#      Those cannot both be the whole story. Donor level will arbitrate.
#
# GATES
#   G1  SANITY (P2's six numbers) + Braak parses in both cohorts
#   G2  enough donors per Braak stage to compare at all (>=5)
#   G3  LEVERAGE: are bins 0.3/0.6 donor-dominated?
#   G4  REGION vs RULER: MCT4 ~ Braak slope, sign agreement across cohorts
#   G5  REACTIVITY -> MCT4: direction, both cohorts, Braak-controlled
#
# set.seed(42).  Run with the working directory set to the repository root.
# =============================================================================

# Paths are RELATIVE to the repository root. Set the working directory there,
# or override with the environment variables SEAAD_H5AD / ROSMAP_ASTRO /
# ROSMAP_CLIN / P2_OUT_DIR. Raw data are not redistributable; see README.md.

suppressPackageStartupMessages({ library(rhdf5); library(Matrix); library(data.table) })
set.seed(42)
OUT <- Sys.getenv("P2_OUT_DIR", unset = "output/tables"); dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
SEA  <- Sys.getenv("SEAAD_H5AD", unset = "data/raw/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad")
ROS  <- Sys.getenv("ROSMAP_ASTRO", unset = "data/raw/ROSMAP/astrocytes.h5Seurat")
RCLIN<- Sys.getenv("ROSMAP_CLIN", unset = "data/raw/ROSMAP/ROSMAP_Green_clinical_merged.csv")
BS <- 5000L
MIN_DONORS_STAGE <- 5L
BINS_S <- round(seq(0.2,0.9,0.1),1); KEY_S <- sprintf("%.1f", BINS_S); NBS <- 8L
I_E <- 1:3; I_L <- 5:7

DRIVER <- "SLC16A3"
REACT  <- c("GFAP","VIM","SERPINA3","CD44","OSMR","C3","SERPING1","GBP2",
            "MT1X","MT2A","HSPB1","CSTB","S100A10","EMP1")
HOMEO  <- c("SLC1A2","SLC1A3","GLUL","AQP4","KCNJ10","GJA1","SLC6A11")
ANLS   <- c("SLC16A3","HK2","LDHA","SLC2A1","SLC16A1")
FOCUS  <- unique(c(REACT,HOMEO,ANLS,"PTGDS","GAPDH","PKM","BSG","GJB6","ATP1A2"))

sink(file.path(OUT,"47_braak.txt"), split=TRUE)
cat("=== 47: Braak as the common axis. CPS and pseudotime are NOT the same ruler. ===\n\n")

parse_braak <- function(x) {
  if (is.numeric(x)) return(as.integer(round(x)))
  s <- toupper(trimws(as.character(x)))
  s <- gsub("BRAAK", "", s); s <- trimws(s)
  r <- c("0"=0L,"I"=1L,"II"=2L,"III"=3L,"IV"=4L,"V"=5L,"VI"=6L)
  out <- r[s]
  num <- suppressWarnings(as.integer(s))
  out[is.na(out) & !is.na(num)] <- num[is.na(out) & !is.na(num)]
  as.integer(out) }

## ======================= [1] SEA-AD ========================================
cat("[1] SEA-AD MTG astrocytes\n")
rd <- function(h5,k){
  v <- tryCatch({ ca <- as.vector(h5read(h5,paste0("obs/__categories/",k)))
                  co <- as.integer(h5read(h5,paste0("obs/",k))); co[co<0] <- NA; ca[co+1L] },
                error=function(e) NULL)
  if (!is.null(v)) return(v)
  tryCatch(as.vector(h5read(h5,paste0("obs/",k))), error=function(e) NULL) }
gS  <- as.character(h5read(SEA,"var/_index")); nGS <- length(gS)
ipS <- as.numeric(h5read(SEA,"X/indptr",bit64conversion="double")); nCS <- length(ipS)-1L
cpsS<- suppressWarnings(as.numeric(rd(SEA,"Continuous Pseudo-progression Score")))
scS <- as.character(rd(SEA,"Subclass")); dnS <- as.character(rd(SEA,"Donor ID"))
brS_raw <- rd(SEA,"Braak"); ceS_raw <- rd(SEA,"CERAD score"); cgS_raw <- rd(SEA,"Cognitive Status")
cat("  raw Braak values in SEA-AD: "); print(head(sort(unique(as.character(brS_raw))), 12))
cat("  raw CERAD values          : "); print(head(sort(unique(as.character(ceS_raw))), 8))
cat("  raw Cognitive Status      : "); print(head(sort(unique(as.character(cgS_raw))), 8))
brS <- parse_braak(brS_raw)
cat("  parsed Braak: "); print(table(brS[grepl("^Astro",scS)], useNA="ifany"))

isA <- grepl("^Astro", scS, ignore.case=TRUE)
biS <- match(sprintf("%.1f", round(cpsS,1)), KEY_S)
donS<- sort(unique(dnS[is.finite(cpsS)])); nDS <- length(donS); diS <- match(dnS, donS)
kS  <- isA & is.finite(cpsS) & !is.na(diS)
cat(sprintf("  astrocytes %d | donors %d\n", sum(kS), nDS))

SdS <- Matrix(0,nDS,nGS,sparse=TRUE); CdS <- numeric(nDS)
SbS <- Matrix(0,NBS,nGS,sparse=TRUE); CbS <- numeric(NBS)
DBc <- matrix(0L, NBS, nDS)                       # bin x donor cell counts -> leverage
t0 <- Sys.time()
for (s0 in seq(1L,nCS,by=BS)) {
  e0 <- min(s0+BS-1L,nCS); gi <- (s0:e0)[kS[s0:e0]]; if (!length(gi)) next
  sp <- ipS[s0]; cnt <- ipS[e0+1L]-sp; if (cnt<=0) next
  ci <- h5read(SEA,"X/indices",start=sp+1,count=cnt,bit64conversion="double")
  cd <- as.numeric(h5read(SEA,"X/data",start=sp+1,count=cnt))
  np <- as.integer(diff(ipS[s0:(e0+1L)])); cof <- rep(s0:e0, times=np); m <- kS[cof]
  if (!any(m)) { rm(ci,cd); next }
  loc <- match(cof[m], gi)
  M <- sparseMatrix(i=loc, j=as.integer(ci[m])+1L, x=cd[m], dims=c(length(gi),nGS))
  fd <- factor(diS[gi], levels=1:nDS); Dd <- fac2sparse(fd, drop.unused.levels=FALSE)
  SdS <- SdS + Dd %*% M; CdS <- CdS + as.numeric(rowSums(Dd))
  ib <- biS[gi]; ok <- !is.na(ib)
  if (any(ok)) {
    fb <- factor(ib[ok], levels=1:NBS); Db <- fac2sparse(fb, drop.unused.levels=FALSE)
    SbS <- SbS + Db %*% M[ok,,drop=FALSE]; CbS <- CbS + as.numeric(rowSums(Db))
    tb <- table(factor(ib[ok], levels=1:NBS), factor(diS[gi][ok], levels=1:nDS))
    DBc <- DBc + matrix(as.integer(tb), NBS, nDS)
  }
  rm(M,ci,cd); cat(sprintf("    %d / %d\r", e0, nCS))
}
cat(sprintf("\n  done (%.1f min)\n", as.numeric(difftime(Sys.time(),t0,units="mins")))); gc()
DMS <- as.matrix(SdS)/pmax(CdS,1); colnames(DMS) <- gS
BMS <- as.matrix(SbS)/pmax(CbS,1); colnames(BMS) <- gS
dcpsS<- as.numeric(tapply(cpsS[kS], factor(dnS[kS], levels=donS), mean))
dbrS <- as.integer(tapply(brS[kS],  factor(dnS[kS], levels=donS), function(x) x[1]))
dceS <- as.character(tapply(as.character(ceS_raw)[kS], factor(dnS[kS], levels=donS), function(x) x[1]))
dcgS <- as.character(tapply(as.character(cgS_raw)[kS], factor(dnS[kS], levels=donS), function(x) x[1]))
uS   <- which(CdS >= 20 & is.finite(dcpsS))

pctv <- function(v){ e<-mean(v[I_E]); l<-mean(v[I_L]); if(!is.finite(e)||abs(e)<1e-9) NA_real_ else (l-e)/e*100 }
G1T <- data.table(gene=c("SLC16A3","HK2","LDHA","PDK1","SLC16A1","SLC2A1"),
                  pub=c(-43.2,-35.2,-20.8,-16.5,-11.1,-7.4))
G1T[, here := sapply(gene, function(g) round(pctv(BMS[,g]),1))]
G1T[, ok := abs(here-pub) <= 1]
G1a <- isTRUE(all(G1T$ok)); print(G1T)
cat(sprintf("  SANITY >>> %s\n", if(G1a) "PASS" else "FAIL")); if(!G1a){sink();stop("G1 FAIL")}

## ======================= [2] ROSMAP ========================================
cat("\n[2] ROSMAP DLPFC astrocytes\n")
CL <- fread(RCLIN); setnames(CL, tolower(names(CL)))
cat("  ROSMAP braaksc: "); print(table(CL$braaksc, useNA="ifany"))
cat("  ROSMAP ceradsc: "); print(table(CL$ceradsc, useNA="ifany"))
cat("  ROSMAP cogdx  : "); print(table(CL$cogdx,   useNA="ifany"))
gR  <- as.vector(h5read(ROS,"assays/RNA/features")); nGR <- length(gR)
ipR <- as.numeric(h5read(ROS,"assays/RNA/counts/indptr",bit64conversion="double"))
fl  <- c("individualID","state","nCount_RNA")
MDR <- as.data.table(lapply(fl, function(f) as.vector(h5read(ROS,paste0("meta.data/",f)))))
setnames(MDR, fl); MDR[, cell := 1:.N]; nCR <- nrow(MDR)
MDR <- merge(MDR, CL[, .(individualid, pseudotime, braaksc, ceradsc, cogdx)],
             by.x="individualID", by.y="individualid", all.x=TRUE, sort=FALSE)
setorder(MDR, cell)
kR  <- is.finite(MDR$braaksc) & is.finite(MDR$nCount_RNA) & MDR$nCount_RNA>0
dnR <- sort(unique(MDR$individualID[kR])); nDR <- length(dnR); diR <- match(MDR$individualID, dnR)
kR  <- kR & !is.na(diR)
cat(sprintf("  astrocytes %d | donors with Braak %d\n", sum(kR), nDR))
ncR <- as.numeric(MDR$nCount_RNA)
SdR <- Matrix(0,nDR,nGR,sparse=TRUE); CdR <- numeric(nDR)
t0 <- Sys.time()
for (s0 in seq(1L,nCR,by=BS)) {
  e0 <- min(s0+BS-1L,nCR); gi <- (s0:e0)[kR[s0:e0]]; if (!length(gi)) next
  sp <- ipR[s0]; cnt <- ipR[e0+1L]-sp; if (cnt<=0) next
  ci <- h5read(ROS,"assays/RNA/counts/indices",start=sp+1,count=cnt,bit64conversion="double")
  cd <- as.numeric(h5read(ROS,"assays/RNA/counts/data",start=sp+1,count=cnt))
  np <- as.integer(diff(ipR[s0:(e0+1L)])); cof <- rep(s0:e0, times=np); m <- kR[cof]
  if (!any(m)) { rm(ci,cd); next }
  loc <- match(cof[m], gi)
  M <- sparseMatrix(i=loc, j=as.integer(ci[m])+1L, x=cd[m], dims=c(length(gi),nGR))
  nc <- ncR[gi]; M@x <- log1p(M@x / nc[M@i + 1L] * 1e4)
  fd <- factor(diR[gi], levels=1:nDR); Dd <- fac2sparse(fd, drop.unused.levels=FALSE)
  SdR <- SdR + Dd %*% M; CdR <- CdR + as.numeric(rowSums(Dd))
  rm(M,ci,cd); cat(sprintf("    %d / %d\r", e0, nCR))
}
cat(sprintf("\n  done (%.1f min)\n", as.numeric(difftime(Sys.time(),t0,units="mins")))); gc()
DMR <- as.matrix(SdR)/pmax(CdR,1); colnames(DMR) <- gR
dbrR <- as.integer(tapply(MDR$braaksc[kR], factor(MDR$individualID[kR], levels=dnR), function(x) x[1]))
dptR <- as.numeric(tapply(MDR$pseudotime[kR], factor(MDR$individualID[kR], levels=dnR), function(x) x[1]))
dceR <- as.numeric(tapply(MDR$ceradsc[kR], factor(MDR$individualID[kR], levels=dnR), function(x) x[1]))
uR   <- which(CdR >= 20 & is.finite(dbrR))
cat(sprintf("  usable donors: SEA-AD %d | ROSMAP %d\n", length(uS), length(uR)))

## ======================= G2: enough donors per stage? ======================
cat("\n[G2] donors per Braak stage\n")
tS <- table(factor(dbrS[uS], levels=0:6)); tR <- table(factor(dbrR[uR], levels=0:6))
ST <- data.table(braak=0:6, SEAAD=as.integer(tS), ROSMAP=as.integer(tR))
ST[, both_usable := SEAAD >= MIN_DONORS_STAGE & ROSMAP >= MIN_DONORS_STAGE]
print(ST); fwrite(ST, file.path(OUT,"47_braak_counts.csv"))
G2 <- sum(ST$both_usable) >= 3
cat(sprintf("  G2 >>> %s (%d stages comparable)\n",
            if (G2) "PASS" else "FAIL - too few shared stages; the axes cannot be aligned",
            sum(ST$both_usable)))

## ======================= G3: LEVERAGE at bins 0.3 / 0.6 ====================
cat("\n[G3] who owns bins 0.3 and 0.6? (reactive markers crash there in lockstep)\n")
LEV <- rbindlist(lapply(1:NBS, function(b){
  cnt <- DBc[b,]; tot <- sum(cnt); if (!tot) return(NULL)
  o <- order(cnt, decreasing=TRUE)
  data.table(bin=BINS_S[b], n_cells=tot, n_donors=sum(cnt>0),
             top1_share=round(cnt[o[1]]/tot,3),
             top2_share=round(sum(cnt[o[1:2]])/tot,3),
             top3_share=round(sum(cnt[o[1:3]])/tot,3),
             top_donor=donS[o[1]]) }))
print(LEV); fwrite(LEV, file.path(OUT,"47_bin_leverage.csv"))
lev36 <- LEV[bin %in% c(0.3,0.6)]
oth   <- LEV[!bin %in% c(0.3,0.6)]
G3 <- isTRUE(mean(lev36$top3_share) > 1.3*mean(oth$top3_share))
cat(sprintf("\n  top-3 donor share: bins 0.3/0.6 = %.3f | other bins = %.3f\n",
            mean(lev36$top3_share), mean(oth$top3_share)))
cat(sprintf("  G3 >>> %s\n", if (G3)
  "LEVERAGE CONFIRMED. Bins 0.3/0.6 are donor-dominated. The reactive-marker\n         crashes there -- and the 'GFAP peaks at bin 0.5' reading -- are NOT\n         trustworthy. That observation is WITHDRAWN."
  else "not obviously donor-dominated. The bin-level swings need another explanation."))

## ======================= G4: REGION or RULER? ==============================
cat("\n[G4] MCT4 ~ Braak, both cohorts, donor-level\n")
fit <- function(y, x, lab, n_lab) {
  ok <- is.finite(y) & is.finite(x); if (sum(ok) < 15) return(NULL)
  f <- lm(y[ok] ~ x[ok]); co <- summary(f)$coefficients
  data.table(cohort=lab, axis=n_lab, n=sum(ok), beta=round(co[2,1],4), p=signif(co[2,4],3)) }
AX <- rbindlist(list(
  fit(DMS[uS,DRIVER], dcpsS[uS],           "SEA-AD MTG",   "CPS"),
  fit(DMS[uS,DRIVER], as.numeric(dbrS[uS]),"SEA-AD MTG",   "Braak"),
  fit(DMR[uR,DRIVER], dptR[uR],            "ROSMAP DLPFC", "pseudotime"),
  fit(DMR[uR,DRIVER], as.numeric(dbrR[uR]),"ROSMAP DLPFC", "Braak")), fill=TRUE)
print(AX); fwrite(AX, file.path(OUT,"47_mct4_axes.csv"))
bS <- AX[cohort=="SEA-AD MTG"   & axis=="Braak"]$beta
bR <- AX[cohort=="ROSMAP DLPFC" & axis=="Braak"]$beta
G4 <- isTRUE(is.finite(bS) && is.finite(bR) && sign(bS)==sign(bR))

# and the full chain on the Braak axis, both cohorts
CH <- rbindlist(lapply(intersect(FOCUS, intersect(gS,gR)), function(g){
  a <- fit(DMS[uS,g], as.numeric(dbrS[uS]), "SEA","Braak")
  b <- fit(DMR[uR,g], as.numeric(dbrR[uR]), "ROS","Braak")
  if (is.null(a)||is.null(b)) return(NULL)
  data.table(gene=g, MTG_beta=a$beta, MTG_p=a$p, DLPFC_beta=b$beta, DLPFC_p=b$p,
             same_sign=sign(a$beta)==sign(b$beta)) }))
setorder(CH, MTG_beta); print(CH); fwrite(CH, file.path(OUT,"47_chain_braak.csv"))

sink(file.path(OUT,"47_VERDICT.txt"))
cat("=== ON THE SHARED NEUROPATHOLOGY AXIS: REGION EFFECT, OR WRONG RULER? ===\n\n")
print(AX); cat("\n")
cat(sprintf("  MCT4 ~ Braak   MTG   beta %+.4f  p %.3g\n", bS, AX[cohort=="SEA-AD MTG" & axis=="Braak"]$p))
cat(sprintf("  MCT4 ~ Braak   DLPFC beta %+.4f  p %.3g\n", bR, AX[cohort=="ROSMAP DLPFC" & axis=="Braak"]$p))
cat(sprintf("\n  G4 >>> %s\n", if (G4)
  "SAME SIGN on Braak. The CPS-vs-pseudotime discrepancy was a RULER problem,\n         not a contradiction. The two cohorts agree once the axis is shared."
  else "OPPOSITE SIGN even on Braak. This is NOT a ruler artefact. MTG and DLPFC\n         astrocytes genuinely differ. P2's finding is REGION-SPECIFIC and must be\n         stated as such -- with MTG's intermediate vulnerability as the reason."))
cat("\n  chain agreement on Braak (same sign in both cohorts):\n")
cat(sprintf("    %d / %d genes\n", sum(CH$same_sign, na.rm=TRUE), nrow(CH)))
print(CH[gene %in% c("SLC16A3","HK2","LDHA","GFAP","SERPINA3","CD44","SLC1A2","PTGDS")])
sink(); writeLines(readLines(file.path(OUT,"47_VERDICT.txt")))

## ======================= G5: the CEO's hypothesis ==========================
cat("\n[G5] does a MORE REACTIVE donor have LESS MCT4? (Braak-controlled)\n")
rsc <- function(DM, gl){ gl <- intersect(gl, colnames(DM)); rowMeans(DM[, gl, drop=FALSE]) }
RS_S <- rsc(DMS, REACT); RS_R <- rsc(DMR, REACT)
HS_S <- rsc(DMS, HOMEO); HS_R <- rsc(DMR, HOMEO)

react_mct4 <- function(DM, react, homeo, braak, prog, u, lab) {
  m <- DM[u, DRIVER]; r <- react[u]; h <- homeo[u]; b <- as.numeric(braak[u]); p <- prog[u]
  ok <- is.finite(m)&is.finite(r)&is.finite(b)
  f0 <- lm(m[ok] ~ r[ok]);              c0 <- summary(f0)$coefficients
  f1 <- lm(m[ok] ~ r[ok] + b[ok]);      c1 <- summary(f1)$coefficients
  f2 <- lm(h[ok] ~ r[ok] + b[ok]);      c2 <- summary(f2)$coefficients
  data.table(cohort=lab, n=sum(ok),
             MCT4_on_react=round(c0[2,1],4),        p_raw=signif(c0[2,4],3),
             MCT4_on_react_adjBraak=round(c1[2,1],4), p_adj=signif(c1[2,4],3),
             HOMEO_on_react_adjBraak=round(c2[2,1],4), p_homeo=signif(c2[2,4],3)) }
RM <- rbindlist(list(
  react_mct4(DMS, RS_S, HS_S, dbrS, dcpsS, uS, "SEA-AD MTG"),
  react_mct4(DMR, RS_R, HS_R, dbrR, dptR,  uR, "ROSMAP DLPFC")), fill=TRUE)
print(RM); fwrite(RM, file.path(OUT,"47_reactivity_mct4.csv"))

# mediation, both directions, on the Braak axis
med <- function(DM, react, braak, u, lab){
  y <- DM[u,DRIVER]; r <- react[u]; b <- as.numeric(braak[u])
  ok <- is.finite(y)&is.finite(r)&is.finite(b)
  f0 <- lm(y[ok] ~ b[ok]); f1 <- lm(y[ok] ~ b[ok] + r[ok])
  g0 <- lm(r[ok] ~ b[ok]); g1 <- lm(r[ok] ~ b[ok] + y[ok])
  b0 <- coef(f0)[2]; b1 <- coef(f1)[2]; a0 <- coef(g0)[2]; a1 <- coef(g1)[2]
  data.table(cohort=lab,
             react_mediates_MCT4 = round(100*(b0-b1)/b0,1),
             MCT4_mediates_react = round(100*(a0-a1)/a0,1)) }
MD <- rbindlist(list(med(DMS,RS_S,dbrS,uS,"SEA-AD MTG"),
                     med(DMR,RS_R,dbrR,uR,"ROSMAP DLPFC")), fill=TRUE)
print(MD); fwrite(MD, file.path(OUT,"47_reactivity_mediation.csv"))

sink(file.path(OUT,"47_VERDICT.txt"), append=TRUE)
cat("\n\n=== DOES REACTIVITY EXPLAIN THE MCT4 LOSS? ===\n\n")
print(RM); cat("\n"); print(MD)
cat("\n  The CEO's model: astrocytes turn reactive -> they stop supporting the neuron.\n")
cat("  Read the sign of MCT4_on_react_adjBraak:\n")
cat("    NEGATIVE  -> reactive astrocytes carry LESS MCT4. The model holds and the\n")
cat("                 target is the reactive transition itself (NF-kB / STAT3 axis).\n")
cat("    POSITIVE  -> reactive astrocytes carry MORE MCT4 (Warburg). Then MCT4 is\n")
cat("                 not what they lose -- look at HOMEO_on_react_adjBraak instead:\n")
cat("                 if THAT is negative, they make lactate but stop DELIVERING it.\n")
cat("  If the two cohorts disagree in sign, say so. Do not average them.\n")
sink()
cat("\n=== 47 done. Send back 47_VERDICT.txt, 47_braak.txt, 47_braak_counts.csv,\n")
cat("    47_bin_leverage.csv, 47_mct4_axes.csv, 47_chain_braak.csv,\n")
cat("    47_reactivity_mct4.csv, 47_reactivity_mediation.csv ===\n")
sink()
