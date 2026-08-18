# =============================================================================
# 46_composition_vs_percell.R
#
# THE PROBLEM, STATED BY THE CEO AND HE IS RIGHT.
#
#   Astrogliosis is UNIVERSAL in AD post-mortem brain. GFAP immunostaining is a
#   standard neuropathological finding and it tracks Braak stage. MTG is involved
#   by Braak IV-V. So SEA-AD MTG astrocytes showing GFAP -11.5% (detection-matched
#   percentile 0.4 -- SIGNIFICANTLY DOWN) is not a curiosity. It contradicts
#   neuropathology. Something is wrong with how we are reading this.
#
#   ROSMAP DLPFC, by contrast, is textbook reactive:
#     SERPINA3 +265%  CD44 +142%  OSMR +125%  GBP2 +82%  GFAP +35%  C3 +43%
#   ...and its MCT4/HK2/LDHA/GLUT1 all go UP with it -- exactly what a reactive,
#   Warburg-like astrocyte does.
#
# THE MOST LIKELY EXPLANATION IS COMPOSITION.
#
#   snRNA-seq measures PER-CELL expression. Astrocytes are a mixture of subtypes.
#     mean(bin) = SUM over subtypes [ weight_s(bin) x expression_s(bin) ]
#   If the WEIGHTS shift with disease, the mean moves with NO per-cell change.
#
#   I already caught this exact trap once, on HK2:
#     bulk HK2 = (microglia count) x (per-cell HK2) -> the literature's "HK2 rise"
#     was microgliosis, not upregulation.
#   I did not apply the same decomposition INSIDE the astrocyte. That was my error.
#
#   AND IT CUTS BOTH WAYS. If GFAP -11.5% is composition, then MCT4 -43.2% may be
#   composition too:
#       NOT "astrocytes lose MCT4"
#       BUT "MCT4-high astrocytes disappear"
#   Those are different diseases and different drug targets.
#
# WHAT THIS SCRIPT DOES. Exact shift-share decomposition, both cohorts:
#     mean_early = SUM w_s(early) e_s(early)
#     mean_late  = SUM w_s(late)  e_s(late)
#   counterfactuals:
#     composition-only : SUM w_s(late)  e_s(early)   <- weights moved, cells did not
#     per-cell-only    : SUM w_s(early) e_s(late)    <- cells moved, weights did not
#   plus WITHIN-SUBTYPE donor-level regression, which is the honest test:
#   does MCT4 fall INSIDE a subtype, holding identity fixed?
#
#   P2 claims "MCT4 declined in 4 of 6 astrocyte subtypes". We verify that from raw.
#
# PRE-REGISTERED VERDICT (fixed before any number is seen)
#   G3: per-cell effect / total effect >= 0.50  -> P2's reading survives
#                                      <  0.50  -> COMPOSITION DOMINATES.
#                                                  P2's central claim must be
#                                                  RESTATED, not defended.
#
# GATES
#   G1  SANITY: P2's six published %changes reproduce (delta <= 1)
#   G2  is there a composition shift at all? (subtype weights ~ CPS)
#   G3  THE DECIDER: per-cell vs composition share of the MCT4 decline
#   G4  within-subtype: in how many subtypes does MCT4 actually fall? (P2 says 4/6)
#   G5  ROSMAP: same decomposition. Does Ast.5 (mean pseudotime 0.766) explain
#       the +9.7%?
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
MIN_CELLS_DONOR_SUBTYPE <- 10L      # a donor x subtype cell needs this many nuclei
MIN_DONORS_SUBTYPE      <- 20L      # a subtype needs this many donors to be regressed
PERCELL_SHARE_GATE      <- 0.50

NB <- 8L
BINS_S <- round(seq(0.2,0.9,0.1),1); KEY_S <- sprintf("%.1f", BINS_S)
I_E <- 1:3; I_L <- 5:7

DRIVER <- "SLC16A3"
REACT_CLASSIC <- c("GFAP","VIM","SYNM")
REACT_A1      <- c("C3","SERPING1","GBP2","PSMB8","FKBP5")
REACT_DAA     <- c("SERPINA3","CD44","OSMR","MT1X","MT2A","HSPB1","CSTB")
REACT_A2      <- c("S100A10","EMP1","CD109")
HOMEO         <- c("SLC1A2","SLC1A3","GLUL","AQP4","KCNJ10","GJA1","SLC6A11")
ANLS          <- c("SLC16A3","HK2","LDHA","SLC2A1","GAPDH","PKM","SLC16A1")
REACT_ALL <- unique(c(REACT_CLASSIC,REACT_A1,REACT_DAA,REACT_A2))
FOCUS <- unique(c(REACT_ALL, HOMEO, ANLS, "PTGDS","BSG","GJB6","ATP1A2"))

sink(file.path(OUT,"46_decomposition.txt"), split=TRUE)
cat("=== 46: is the MCT4 decline PER-CELL, or is it COMPOSITION? ===\n\n")

## ============ generic decomposition =========================================
# w_e, w_l : subtype weights early / late   (sum to 1)
# e_e, e_l : subtype mean expression early / late
shift_share <- function(w_e, w_l, e_e, e_l) {
  ok <- is.finite(w_e) & is.finite(w_l) & is.finite(e_e) & is.finite(e_l)
  w_e <- w_e[ok]; w_l <- w_l[ok]; e_e <- e_e[ok]; e_l <- e_l[ok]
  w_e <- w_e/sum(w_e); w_l <- w_l/sum(w_l)
  m_e  <- sum(w_e*e_e); m_l <- sum(w_l*e_l)
  m_c  <- sum(w_l*e_e)      # composition moved, per-cell frozen
  m_p  <- sum(w_e*e_l)      # per-cell moved, composition frozen
  if (!is.finite(m_e) || abs(m_e) < 1e-9) return(rep(NA_real_,4))
  tot  <- (m_l - m_e)/m_e*100
  comp <- (m_c - m_e)/m_e*100
  perc <- (m_p - m_e)/m_e*100
  inter<- (m_l - m_c - m_p + m_e)/m_e*100
  c(total=tot, composition=comp, per_cell=perc, interaction=inter)
}

## ============ SEA-AD reader (h5ad, cell-major, log-norm) ====================
cat("[1] SEA-AD MTG astrocytes\n")
rd_obs <- function(h5,k){
  v <- tryCatch({ ca <- as.vector(h5read(h5,paste0("obs/__categories/",k)))
                  co <- as.integer(h5read(h5,paste0("obs/",k))); co[co<0] <- NA; ca[co+1L] },
                error=function(e) NULL)
  if (!is.null(v)) return(v); as.vector(h5read(h5,paste0("obs/",k))) }
genesS <- as.character(h5read(SEA,"var/_index")); nGS <- length(genesS)
ipS <- as.numeric(h5read(SEA,"X/indptr",bit64conversion="double")); nCS <- length(ipS)-1L
cpsS <- suppressWarnings(as.numeric(rd_obs(SEA,"Continuous Pseudo-progression Score")))
scS  <- as.character(rd_obs(SEA,"Subclass"))
supS <- as.character(rd_obs(SEA,"Supertype"))
dnS  <- as.character(rd_obs(SEA,"Donor ID"))
isA  <- grepl("^Astro", scS, ignore.case=TRUE)
biS  <- match(sprintf("%.1f", round(cpsS,1)), KEY_S)
donS <- sort(unique(dnS[is.finite(cpsS)])); nDS <- length(donS); diS <- match(dnS, donS)
stS  <- sort(unique(supS[isA & !is.na(supS)])); nST_S <- length(stS); siS <- match(supS, stS)
keepS<- isA & is.finite(cpsS) & !is.na(biS) & !is.na(diS) & !is.na(siS)
cat(sprintf("  astrocytes %d | donors %d | subtypes %d: %s\n",
            sum(keepS), nDS, nST_S, paste(stS, collapse=", ")))

K_sb <- (siS-1L)*NB  + biS          # subtype x bin
K_sd <- (siS-1L)*nDS + diS          # subtype x donor
N_sb <- nST_S*NB; N_sd <- nST_S*nDS
Ssb <- Matrix(0,N_sb,nGS,sparse=TRUE); Csb <- numeric(N_sb)
Ssd <- Matrix(0,N_sd,nGS,sparse=TRUE); Csd <- numeric(N_sd)
Db  <- Matrix(0,N_sb,nGS,sparse=TRUE)
t0 <- Sys.time()
for (s0 in seq(1L,nCS,by=BS)) {
  e0 <- min(s0+BS-1L,nCS); gi <- (s0:e0)[keepS[s0:e0]]; if (!length(gi)) next
  sp <- ipS[s0]; cnt <- ipS[e0+1L]-sp; if (cnt<=0) next
  ci <- h5read(SEA,"X/indices",start=sp+1,count=cnt,bit64conversion="double")
  cd <- as.numeric(h5read(SEA,"X/data",start=sp+1,count=cnt))
  np <- as.integer(diff(ipS[s0:(e0+1L)])); cof <- rep(s0:e0, times=np); m <- keepS[cof]
  if (!any(m)) { rm(ci,cd); next }
  loc <- match(cof[m], gi)
  M <- sparseMatrix(i=loc, j=as.integer(ci[m])+1L, x=cd[m], dims=c(length(gi),nGS))
  Mn <- M; Mn@x <- rep(1,length(M@x))
  f1 <- factor(K_sb[gi], levels=1:N_sb); D1 <- fac2sparse(f1, drop.unused.levels=FALSE)
  Ssb <- Ssb + D1 %*% M; Db <- Db + D1 %*% Mn; Csb <- Csb + as.numeric(rowSums(D1))
  f2 <- factor(K_sd[gi], levels=1:N_sd); D2 <- fac2sparse(f2, drop.unused.levels=FALSE)
  Ssd <- Ssd + D2 %*% M; Csd <- Csd + as.numeric(rowSums(D2))
  rm(M,Mn,ci,cd); cat(sprintf("    %d / %d\r", e0, nCS))
}
cat(sprintf("\n  done (%.1f min)\n", as.numeric(difftime(Sys.time(),t0,units="mins")))); gc()
EB <- as.matrix(Ssb)/pmax(Csb,1); colnames(EB) <- genesS     # (subtype x bin) mean
DB <- as.matrix(Db )/pmax(Csb,1); colnames(DB) <- genesS     # detection
ED <- as.matrix(Ssd)/pmax(Csd,1); colnames(ED) <- genesS     # (subtype x donor) mean
rSB <- function(s) (s-1L)*NB + 1:NB
rSD <- function(s) (s-1L)*nDS + 1:nDS
NCB <- sapply(1:nST_S, function(s) Csb[rSB(s)])              # NB x nST  cells
NCD <- sapply(1:nST_S, function(s) Csd[rSD(s)])              # nD x nST
dcpsS <- as.numeric(tapply(cpsS[keepS], factor(dnS[keepS], levels=donS), mean))

## ---- G1 SANITY: reconstruct the pooled bin means from subtypes -------------
pooled <- function(g){ e <- sapply(1:nST_S, function(s) EB[rSB(s), g])   # NB x nST
                       rowSums(e * NCB) / rowSums(NCB) }
pctv <- function(v){ e<-mean(v[I_E]); l<-mean(v[I_L]); if(!is.finite(e)||abs(e)<1e-9) NA_real_ else (l-e)/e*100 }
G1T <- data.table(gene=c("SLC16A3","HK2","LDHA","PDK1","SLC16A1","SLC2A1"),
                  pub =c(-43.2,-35.2,-20.8,-16.5,-11.1,-7.4))
G1T[, here := sapply(gene, function(g) round(pctv(pooled(g)),1))]
G1T[, ok := is.finite(here) & abs(here-pub) <= 1]
G1 <- isTRUE(all(G1T$ok))
cat("\n[G1] rebuilding P2's numbers from the subtype decomposition\n"); print(G1T)
cat(sprintf("  G1 >>> %s\n", if(G1) "PASS - the decomposition is arithmetically consistent with P2."
            else "FAIL - subtype weights do not reconstruct P2. STOP."))
if (!G1) { sink(); stop("G1 FAILED") }

## ---- G2: does the composition shift? --------------------------------------
cat("\n[G2] do astrocyte subtype proportions move with CPS?\n")
Wb <- t(apply(NCB, 1, function(r) r/sum(r)))                 # NB x nST weights
COMP <- data.table(bin=BINS_S); for (s in 1:nST_S) COMP[[stS[s]]] <- round(Wb[,s],4)
COMP[, n_cells := rowSums(NCB)]
print(COMP); fwrite(COMP, file.path(OUT,"46_seaad_composition.csv"))
frac_d <- NCD / pmax(rowSums(NCD),1)                          # nD x nST
uD <- which(rowSums(NCD) >= 20 & is.finite(dcpsS))
CS <- rbindlist(lapply(1:nST_S, function(s){
  f <- lm(frac_d[uD,s] ~ dcpsS[uD]); co <- summary(f)$coefficients
  data.table(subtype=stS[s], n_cells=sum(NCD[,s]),
             frac_early=round(mean(Wb[I_E,s]),4), frac_late=round(mean(Wb[I_L,s]),4),
             slope=round(co[2,1],4), p=signif(co[2,4],3)) }))
setorder(CS, -slope); print(CS); fwrite(CS, file.path(OUT,"46_seaad_subtype_slopes.csv"))
G2 <- any(CS$p < 0.05, na.rm=TRUE)
cat(sprintf("\n  G2 >>> %s\n", if (G2) "composition DOES shift. The decomposition below is necessary."
            else "no significant composition shift. Then per-cell change is the only explanation."))

## ---- reactive score per subtype -------------------------------------------
cat("\n[3] which SEA-AD subtype is the reactive one?\n")
sc_sub <- function(gs){ gs <- intersect(gs, genesS)
  sapply(1:nST_S, function(s) mean(colMeans(EB[rSB(s), gs, drop=FALSE]), na.rm=TRUE)) }
RS <- data.table(subtype=stS, n_cells=colSums(NCD),
                 classic=round(sc_sub(REACT_CLASSIC),3), A1=round(sc_sub(REACT_A1),3),
                 DAA=round(sc_sub(REACT_DAA),3), homeostatic=round(sc_sub(HOMEO),3),
                 MCT4=round(sc_sub("SLC16A3"),4), HK2=round(sc_sub("HK2"),4))
RS[, reactive_index := round(DAA + A1 - homeostatic, 3)]
setorder(RS, -reactive_index); print(RS); fwrite(RS, file.path(OUT,"46_seaad_reactive.csv"))

## ---- G3: THE DECIDER ------------------------------------------------------
cat("\n[G3] SHIFT-SHARE: composition vs per-cell\n")
dec1 <- function(g) {
  w_e <- colMeans(NCB[I_E,,drop=FALSE]); w_l <- colMeans(NCB[I_L,,drop=FALSE])
  e_e <- sapply(1:nST_S, function(s) mean(EB[rSB(s)[I_E], g]))
  e_l <- sapply(1:nST_S, function(s) mean(EB[rSB(s)[I_L], g]))
  shift_share(w_e, w_l, e_e, e_l) }
DEC <- rbindlist(lapply(intersect(FOCUS, genesS), function(g){
  d <- dec1(g); data.table(gene=g, total=round(d[1],1), composition=round(d[2],1),
                           per_cell=round(d[3],1), interaction=round(d[4],1),
                           percell_share=round(abs(d[3])/pmax(abs(d[1]),1e-9),2)) }))
setorder(DEC, total); print(DEC); fwrite(DEC, file.path(OUT,"46_seaad_shiftshare.csv"))
m4 <- DEC[gene==DRIVER]
gf <- DEC[gene=="GFAP"]
G3 <- isTRUE(nrow(m4)==1 && is.finite(m4$percell_share) && m4$percell_share >= PERCELL_SHARE_GATE)

sink(file.path(OUT,"46_VERDICT.txt"))
cat("=== IS P2's MCT4 DECLINE PER-CELL, OR COMPOSITION? ===\n\n")
print(DEC[gene %in% c("SLC16A3","HK2","LDHA","SLC2A1","GFAP","SERPINA3","CD44","SLC1A2","PTGDS")])
cat("\n  MCT4 (SLC16A3):\n")
cat(sprintf("    total        %+7.1f%%\n", m4$total))
cat(sprintf("    composition  %+7.1f%%   (subtype weights moved, cells unchanged)\n", m4$composition))
cat(sprintf("    per-cell     %+7.1f%%   (cells changed, weights frozen)\n", m4$per_cell))
cat(sprintf("    per-cell share of |total| = %.0f%%   (gate: >= %.0f%%)\n",
            100*m4$percell_share, 100*PERCELL_SHARE_GATE))
cat("\n  GFAP (the neuropathology check):\n")
cat(sprintf("    total %+.1f%% | composition %+.1f%% | per-cell %+.1f%%\n",
            gf$total, gf$composition, gf$per_cell))
cat("    Astrogliosis is universal in AD post-mortem cortex. If GFAP's fall is\n")
cat("    COMPOSITION, then SEA-AD is losing GFAP-high astrocytes from the sampled\n")
cat("    pool -- and the same mechanism could produce a fake MCT4 decline.\n")
cat(sprintf("\n  G3 >>> %s\n", if (G3)
  "PER-CELL DOMINATES. P2's reading holds: astrocytes really do lose MCT4,\n         not merely 'MCT4-high astrocytes disappear'."
  else "COMPOSITION DOMINATES. P2's central sentence must be RESTATED:\n         the data show 'MCT4-high astrocytes are depleted', NOT 'astrocytes\n         downregulate MCT4'. Different disease. Different target. Do not defend\n         the old wording."))
sink(); writeLines(readLines(file.path(OUT,"46_VERDICT.txt")))

## ---- G4: within-subtype regression (P2 claims 4/6) -------------------------
cat("\n[G4] within-subtype donor regression: does MCT4 fall INSIDE a subtype?\n")
WS <- rbindlist(lapply(1:nST_S, function(s){
  idx <- which(NCD[,s] >= MIN_CELLS_DONOR_SUBTYPE & is.finite(dcpsS))
  if (length(idx) < MIN_DONORS_SUBTYPE)
    return(data.table(subtype=stS[s], n_donors=length(idx), beta=NA_real_, p=NA_real_,
                      pct=NA_real_, note="too few donors"))
  y <- ED[rSD(s)[idx], DRIVER]; x <- dcpsS[idx]
  f <- lm(y ~ x); co <- summary(f)$coefficients
  v <- EB[rSB(s), DRIVER]
  data.table(subtype=stS[s], n_donors=length(idx), beta=round(co[2,1],4),
             p=signif(co[2,4],3), pct=round(pctv(v),1), note="") }))
print(WS); fwrite(WS, file.path(OUT,"46_seaad_within_subtype.csv"))
nsig <- sum(WS$p < 0.05 & WS$beta < 0, na.rm=TRUE); ntest <- sum(is.finite(WS$p))
cat(sprintf("\n  MCT4 falls significantly in %d of %d testable subtypes.\n", nsig, ntest))
cat("  P2's Methods claims: 'MCT4 declined in 4 of 6 astrocyte subtypes'.\n")
G4 <- isTRUE(nsig >= 3)
cat(sprintf("  G4 >>> %s\n", if (G4) "the within-subtype effect is real and broad."
            else "the within-subtype effect is NARROW. P2's '4 of 6' needs checking."))

## ======================= ROSMAP =============================================
cat("\n\n########## [5] ROSMAP: same decomposition ##########\n")
CL <- fread(RCLIN); setnames(CL, tolower(names(CL)))
gvR <- as.vector(h5read(ROS,"assays/RNA/features")); nGR <- length(gvR)
ipR <- as.numeric(h5read(ROS,"assays/RNA/counts/indptr", bit64conversion="double"))
fl  <- c("individualID","state","nCount_RNA")
MDR <- as.data.table(lapply(fl, function(f) as.vector(h5read(ROS, paste0("meta.data/",f)))))
setnames(MDR, fl); MDR[, cell := 1:.N]; nCR <- nrow(MDR)
MDR <- merge(MDR, CL[, .(individualid, pseudotime)], by.x="individualID",
             by.y="individualid", all.x=TRUE, sort=FALSE); setorder(MDR, cell)
ptR <- MDR$pseudotime
kpR <- is.finite(ptR) & is.finite(MDR$nCount_RNA) & MDR$nCount_RNA > 0
qR  <- quantile(ptR[kpR], seq(0,1,length.out=NB+1), na.rm=TRUE)
biR <- cut(ptR, unique(qR), include.lowest=TRUE, labels=FALSE)
dnR <- sort(unique(MDR$individualID[kpR])); nDR <- length(dnR); diR <- match(MDR$individualID, dnR)
stR <- sort(unique(MDR$state[kpR])); nST_R <- length(stR); siR <- match(MDR$state, stR)
kpR <- kpR & !is.na(biR) & !is.na(diR) & !is.na(siR)
cat(sprintf("  cells %d | donors %d | states %d\n", sum(kpR), nDR, nST_R))
KsbR <- (siR-1L)*NB + biR; KsdR <- (siR-1L)*nDR + diR
NsbR <- nST_R*NB; NsdR <- nST_R*nDR
SsbR <- Matrix(0,NsbR,nGR,sparse=TRUE); CsbR <- numeric(NsbR)
SsdR <- Matrix(0,NsdR,nGR,sparse=TRUE); CsdR <- numeric(NsdR)
ncR <- as.numeric(MDR$nCount_RNA)
t0 <- Sys.time()
for (s0 in seq(1L,nCR,by=BS)) {
  e0 <- min(s0+BS-1L,nCR); gi <- (s0:e0)[kpR[s0:e0]]; if (!length(gi)) next
  sp <- ipR[s0]; cnt <- ipR[e0+1L]-sp; if (cnt<=0) next
  ci <- h5read(ROS,"assays/RNA/counts/indices",start=sp+1,count=cnt,bit64conversion="double")
  cd <- as.numeric(h5read(ROS,"assays/RNA/counts/data",start=sp+1,count=cnt))
  np <- as.integer(diff(ipR[s0:(e0+1L)])); cof <- rep(s0:e0, times=np); m <- kpR[cof]
  if (!any(m)) { rm(ci,cd); next }
  loc <- match(cof[m], gi)
  M <- sparseMatrix(i=loc, j=as.integer(ci[m])+1L, x=cd[m], dims=c(length(gi),nGR))
  nc <- ncR[gi]; M@x <- log1p(M@x / nc[M@i + 1L] * 1e4)
  f1 <- factor(KsbR[gi], levels=1:NsbR); D1 <- fac2sparse(f1, drop.unused.levels=FALSE)
  SsbR <- SsbR + D1 %*% M; CsbR <- CsbR + as.numeric(rowSums(D1))
  f2 <- factor(KsdR[gi], levels=1:NsdR); D2 <- fac2sparse(f2, drop.unused.levels=FALSE)
  SsdR <- SsdR + D2 %*% M; CsdR <- CsdR + as.numeric(rowSums(D2))
  rm(M,ci,cd); cat(sprintf("    %d / %d\r", e0, nCR))
}
cat(sprintf("\n  done (%.1f min)\n", as.numeric(difftime(Sys.time(),t0,units="mins")))); gc()
EBR <- as.matrix(SsbR)/pmax(CsbR,1); colnames(EBR) <- gvR
EDR <- as.matrix(SsdR)/pmax(CsdR,1); colnames(EDR) <- gvR
rSBR <- function(s) (s-1L)*NB + 1:NB
rSDR <- function(s) (s-1L)*nDR + 1:nDR
NCBR <- sapply(1:nST_R, function(s) CsbR[rSBR(s)])
NCDR <- sapply(1:nST_R, function(s) CsdR[rSDR(s)])
dptR <- as.numeric(tapply(ptR[kpR], factor(MDR$individualID[kpR], levels=dnR), mean))

WbR <- t(apply(NCBR,1,function(r) r/sum(r)))
CR <- data.table(octile=1:NB, pt_lo=round(head(unique(qR),-1),3))
for (s in 1:nST_R) CR[[stR[s]]] <- round(WbR[,s],4)
print(CR); fwrite(CR, file.path(OUT,"46_rosmap_composition.csv"))

sc_subR <- function(gs){ gs <- intersect(gs, gvR)
  sapply(1:nST_R, function(s) mean(colMeans(EBR[rSBR(s), gs, drop=FALSE]), na.rm=TRUE)) }
RSR <- data.table(state=stR, n_cells=colSums(NCDR),
                  mean_pt=round(sapply(1:nST_R, function(s) sum(NCBR[,s]*(1:NB))/sum(NCBR[,s])),2),
                  classic=round(sc_subR(REACT_CLASSIC),3), A1=round(sc_subR(REACT_A1),3),
                  DAA=round(sc_subR(REACT_DAA),3), homeostatic=round(sc_subR(HOMEO),3),
                  MCT4=round(sc_subR("SLC16A3"),4))
RSR[, reactive_index := round(DAA + A1 - homeostatic,3)]
setorder(RSR, -reactive_index); print(RSR); fwrite(RSR, file.path(OUT,"46_rosmap_reactive.csv"))

decR <- function(g){
  w_e <- colMeans(NCBR[I_E,,drop=FALSE]); w_l <- colMeans(NCBR[I_L,,drop=FALSE])
  e_e <- sapply(1:nST_R, function(s) mean(EBR[rSBR(s)[I_E], g]))
  e_l <- sapply(1:nST_R, function(s) mean(EBR[rSBR(s)[I_L], g]))
  shift_share(w_e, w_l, e_e, e_l) }
DECR <- rbindlist(lapply(intersect(FOCUS, gvR), function(g){
  d <- decR(g); data.table(gene=g, total=round(d[1],1), composition=round(d[2],1),
                           per_cell=round(d[3],1), interaction=round(d[4],1),
                           percell_share=round(abs(d[3])/pmax(abs(d[1]),1e-9),2)) }))
setorder(DECR, total); print(DECR); fwrite(DECR, file.path(OUT,"46_rosmap_shiftshare.csv"))

WSR <- rbindlist(lapply(1:nST_R, function(s){
  idx <- which(NCDR[,s] >= MIN_CELLS_DONOR_SUBTYPE & is.finite(dptR))
  if (length(idx) < MIN_DONORS_SUBTYPE)
    return(data.table(state=stR[s], n_donors=length(idx), beta=NA_real_, p=NA_real_))
  f <- lm(EDR[rSDR(s)[idx], DRIVER] ~ dptR[idx]); co <- summary(f)$coefficients
  data.table(state=stR[s], n_donors=length(idx), beta=round(co[2,1],4), p=signif(co[2,4],3)) }))
print(WSR); fwrite(WSR, file.path(OUT,"46_rosmap_within_subtype.csv"))

m4R <- DECR[gene==DRIVER]
sink(file.path(OUT,"46_VERDICT.txt"), append=TRUE)
cat("\n\n=== ROSMAP DLPFC: same decomposition ===\n")
print(DECR[gene %in% c("SLC16A3","HK2","LDHA","GFAP","SERPINA3","CD44","SLC1A2")])
cat(sprintf("\n  MCT4 total %+.1f%% = composition %+.1f%% + per-cell %+.1f%% + interaction %+.1f%%\n",
            m4R$total, m4R$composition, m4R$per_cell, m4R$interaction))
cat("\n  within-state MCT4 slopes:\n"); print(WSR)
cat("\n  READ IT LIKE THIS:\n")
cat("    If ROSMAP's +9.7% is mostly COMPOSITION (a reactive state expanding) and\n")
cat("    the within-state slopes are flat or negative, then ROSMAP never contradicted\n")
cat("    P2 -- it was measuring a different thing (which cells are present), and the\n")
cat("    two cohorts can both be right.\n")
cat("    If ROSMAP's +9.7% is PER-CELL, then astrocytes in DLPFC genuinely raise MCT4\n")
cat("    while MTG astrocytes lower it, and P2's mechanism is region-specific.\n")
sink()
writeLines(readLines(file.path(OUT,"46_VERDICT.txt")))

cat("\n=== 46 done. Send back 46_VERDICT.txt, 46_decomposition.txt,\n")
cat("    46_seaad_shiftshare.csv, 46_seaad_within_subtype.csv, 46_seaad_reactive.csv,\n")
cat("    46_seaad_composition.csv, 46_seaad_subtype_slopes.csv,\n")
cat("    46_rosmap_shiftshare.csv, 46_rosmap_within_subtype.csv, 46_rosmap_reactive.csv ===\n")
sink()
