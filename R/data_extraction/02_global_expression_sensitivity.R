# =============================================================================
# 02_global_expression_sensitivity.R
#
# Sensitivity analysis: do the donor-level cross-cellular couplings
# (astrocytic MCT4 -> neuronal V-ATPase / LAMP1 / LDHB) survive control for
# transcriptome-wide global transcriptional activity, beyond CPS?
#
# Three nested control sets:
#   [A] CPS only            (primary analysis; donor-level partial r = +0.466)
#   [B] CPS + log UMI       (library depth)
#   [C] CPS + GLOBAL MEAN   (transcriptome-wide mean log-normalized expression,
#                            per donor within each cell type — the genome-wide
#                            control reported in the manuscript)
#
# The per-cell global mean (mean log-normalized expression over ALL genes,
# including zeros) is computed in the SAME block-reader pass as the marker
# panel (no extra pass), then aggregated to donor level and used as a covariate.
#
# Input:  SEA-AD h5ad file (requires data access agreement)
#         -> https://portal.brain-map.org
# Output (written to output/tables/):
#   donor_merged_global.csv        — donor-level MCT4, neuron markers, gmean, logUMI
#   crosscell_global_results.csv   — zero-order r and partial r under [A]/[B]/[C]
#
# Run: Rscript R/data_extraction/02_global_expression_sensitivity.R
#      (requires rhdf5; large file, run where the h5ad is available locally)
# =============================================================================
suppressMessages({ library(rhdf5) })
options(stringsAsFactors = FALSE)
set.seed(42)

## ---- Paths (edit these) ----------------------------------------------------
H5      <- "path/to/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad"
OUTDIR  <- "output/tables"
CPS_MIN <- 0.1
BS      <- 100000L      # cell block size; lower to 4000-5000 if memory-limited

# Cell-type patterns (partial match against the obs cell-type label).
# Inspect the [diagnostic] label distribution printed below and adjust if needed.
ASTRO_PATTERN  <- "Astro"
NEURON_PATTERN <- "IT|ET|CT|L5/6 NP|L6b|L6 IT Car3|Exc|Glut"

# Marker panels: supply side (astrocyte) and demand side (neuron)
ASTRO_PANEL  <- c("SLC16A3")                                 # MCT4 (supply-side hub)
NEURON_VATP  <- c("ATP6V1A","ATP6V1B2","ATP6V1C1","ATP6V1E1","ATP6V1H",
                  "ATP6V0A1","ATP6V0B","ATP6V0C","ATP6V0D1","ATP6V0E1")  # V-ATPase x10
NEURON_PANEL <- c(NEURON_VATP, "LDHB", "LAMP1")              # + LDHB, LAMP1

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

## ---- Helpers ---------------------------------------------------------------
read_obs <- function(H5, key){
  base <- paste0("obs/", key)
  v <- tryCatch({
    co <- as.integer(h5read(H5, paste0(base,"/codes"))); co[co < 0] <- NA
    ca <- as.vector(h5read(H5, paste0(base,"/categories"))); ca[co + 1L]
  }, error=function(e) NULL)
  if(!is.null(v)) return(v)
  v <- tryCatch({
    ca <- as.vector(h5read(H5, paste0("obs/__categories/", key)))
    co <- as.integer(h5read(H5, base)); co[co < 0] <- NA
    ca[co + 1L]
  }, error=function(e) NULL)
  if(!is.null(v)) return(v)
  as.vector(h5read(H5, base))
}
pick_obs_key <- function(H5, cands){
  keys <- tryCatch(h5ls(H5)$name[h5ls(H5)$group=="/obs"], error=function(e) character(0))
  for(k in cands) if(k %in% keys) return(k)
  all_obs <- tryCatch(unique(sub("^/obs/?","",h5ls(H5,recursive=TRUE)$group)), error=function(e) character(0))
  for(k in cands) if(k %in% all_obs) return(k)
  NA_character_
}
get_genes <- function(H5){
  for(f in c("var/feature_name","var/gene_name","var/Gene","var/gene_symbols","var/_index","var/index")){
    g <- tryCatch(as.vector(h5read(H5, f)), error=function(e) NULL)
    if(!is.null(g) && length(g) > 1000) { attr(g,"field")<-f; return(g) }
  }
  stop("Could not find a var gene-symbol field — inspect h5ls(H5)")
}

## ---- [1] Structure + obs (once) --------------------------------------------
cat("=== [1] Structure check ===\n")
genes  <- get_genes(H5); nG <- length(genes)
indptr <- as.numeric(h5read(H5, "X/indptr", bit64conversion="double"))   # 64-bit: nnz overflows int32
stopifnot(all(diff(indptr) >= 0))
nnz <- indptr[length(indptr)]
csc <- (length(indptr) == nG + 1L); csr <- !csc
nC  <- length(indptr) - 1L
if(!csr) stop("This script assumes CSR (rows = cells). CSC layout not handled.")
samp <- tryCatch(h5read(H5,"X/data", index=list(1:2000)),
                 error=function(e) h5read(H5,"X/data",start=1,count=2000))
is_counts <- all(abs(samp - round(samp)) < 1e-8)
cat(sprintf("  genes %d | cells %d | nnz %s | X=%s\n", nG, nC,
            format(nnz, scientific=FALSE),
            if(is_counts) "raw counts (CP10k+log1p applied)" else "normalized (used as-is)"))

cat("=== [2] Reading obs ===\n")
CPS_KEY   <- pick_obs_key(H5, c("Continuous Pseudo-progression Score","Continuous_Pseudo_progression_Score","CPS","Pseudo-progression"))
CLASS_KEY <- pick_obs_key(H5, c("Subclass","subclass","Class","class","cell_type","CellType"))
DONOR_KEY <- pick_obs_key(H5, c("Donor ID","donor_id","external_donor_name","Donor","donor","individualID","sample_id","Sample ID"))
cat(sprintf("  CPS=%s | cell-type=%s | donor=%s\n", CPS_KEY, CLASS_KEY, DONOR_KEY))
if(is.na(CPS_KEY)||is.na(CLASS_KEY)||is.na(DONOR_KEY)) stop("Missing obs key — inspect h5ls(H5)")
cps_all   <- suppressWarnings(as.numeric(read_obs(H5, CPS_KEY)))
class_all <- as.character(read_obs(H5, CLASS_KEY))
donor_all <- as.character(read_obs(H5, DONOR_KEY))
stopifnot(length(cps_all)==nC, length(class_all)==nC, length(donor_all)==nC)
cat("  [diagnostic] top-12 cell-type labels:\n"); print(head(sort(table(class_all), decreasing=TRUE), 12))

## ---- [3] Per-cell-type donor aggregation (panel log-norm + gmean + logUMI) --
# gmean = per-cell mean log-normalized expression over ALL genes (zeros included)
#         -> donor mean = transcriptome-wide global transcriptional activity
extract_ct <- function(pattern, panel, tag){
  sel <- which(grepl(pattern, class_all, ignore.case=TRUE) & is.finite(cps_all) & cps_all >= CPS_MIN)
  cat(sprintf("  [%s] '%s' & CPS>=%.2f -> %d cells\n", tag, pattern, CPS_MIN, length(sel)))
  if(length(sel) < 200) stop(sprintf("[%s] too few cells — check pattern", tag))
  found <- panel[panel %in% genes]; miss <- setdiff(panel, genes)
  if(length(miss)) cat(sprintf("    [note] not found, dropped: %s\n", paste(miss, collapse=", ")))
  pidx0 <- match(found, genes) - 1L

  rmap <- rep(NA_integer_, nC); rmap[sel] <- seq_along(sel)
  M      <- matrix(0, nrow=length(sel), ncol=length(found), dimnames=list(NULL, found))
  totUMI <- numeric(length(sel)); gmean <- numeric(length(sel))
  t0 <- Sys.time()
  for(s0 in seq(1L, nC, by=BS)){
    e0 <- min(s0+BS-1L, nC)
    bk <- which(!is.na(rmap[s0:e0])); if(length(bk)==0) next
    sp <- indptr[s0]; cnt <- indptr[e0+1L] - sp; if(cnt <= 0) next
    ci <- h5read(H5,"X/indices", start=sp+1, count=cnt, bit64conversion="double")
    cd <- h5read(H5,"X/data",    start=sp+1, count=cnt)
    for(k in bk){
      g <- s0 + k - 1L
      a <- indptr[g] - sp + 1; b <- indptr[g+1L] - sp; if(b < a) next
      cig <- ci[a:b]; cdg <- cd[a:b]
      if(is_counts){
        tu <- sum(cdg); totUMI[rmap[g]] <- tu
        gmean[rmap[g]] <- sum(log1p(cdg / tu * 1e4)) / nG   # zeros: log1p(0)=0
        h <- match(pidx0, cig); ok <- !is.na(h)
        if(any(ok)) M[rmap[g], ok] <- log1p(cdg[h[ok]] / tu * 1e4)
      } else {
        totUMI[rmap[g]] <- sum(expm1(cdg))                  # approx depth if normalized
        gmean[rmap[g]]  <- sum(cdg) / nG                     # already in log space
        h <- match(pidx0, cig); ok <- !is.na(h)
        if(any(ok)) M[rmap[g], ok] <- cdg[h[ok]]
      }
    }
    cat(sprintf("    [%s] ...%d/%d\r", tag, e0, nC))
  }
  cat(sprintf("\n    [%s] done (%.1f min)\n", tag, as.numeric(difftime(Sys.time(),t0,units="mins"))))

  d <- donor_all[sel]; cpsv <- cps_all[sel]
  agg <- function(v) tapply(v, d, mean, na.rm=TRUE)
  donors <- names(agg(gmean))
  out <- data.frame(donor_id = donors)
  for(gname in found) out[[gname]] <- as.numeric(agg(M[, gname]))
  out$gmean   <- as.numeric(agg(gmean))
  out$logUMI  <- as.numeric(agg(log(pmax(totUMI,1))))
  out$cps     <- as.numeric(agg(cpsv))
  out$n_cells <- as.numeric(tapply(rep(1L,length(sel)), d, sum))
  out
}

cat("=== [3] Extract astrocytes ===\n");        astro  <- extract_ct(ASTRO_PATTERN,  ASTRO_PANEL,  "Astro")
cat("=== [4] Extract excitatory neurons ===\n"); neuron <- extract_ct(NEURON_PATTERN, NEURON_PANEL, "Neuron")

## ---- [5] Merge + neuronal V-ATPase composite -------------------------------
vcols <- intersect(NEURON_VATP, names(neuron))
neuron$VATPase_n <- rowMeans(neuron[, vcols, drop=FALSE], na.rm=TRUE)
neuron2 <- neuron[, c("donor_id","VATPase_n","LDHB","LAMP1","gmean","logUMI","cps","n_cells")]
names(neuron2) <- c("donor_id","VATPase_n","LDHB_n","LAMP1_n","gmean_n","logUMI_n","cps_n","ncell_n")
astro2 <- astro[, c("donor_id","SLC16A3","gmean","logUMI","cps","n_cells")]
names(astro2) <- c("donor_id","MCT4","gmean_a","logUMI_a","cps_a","ncell_a")
df <- merge(astro2, neuron2, by="donor_id")
cat(sprintf("=== [5] merged donors: %d ===\n", nrow(df)))
write.csv(df, file.path(OUTDIR, "donor_merged_global.csv"), row.names=FALSE)

## ---- [6] Partial correlation: [A] CPS / [B] +logUMI / [C] +gmean -----------
pcorZ <- function(x, y, Z){
  ok <- complete.cases(x, y, Z)
  if(sum(ok) < 10) return(c(r=NA, p=NA, n=sum(ok)))
  rx <- residuals(lm(x[ok] ~ Z[ok, , drop=FALSE]))
  ry <- residuals(lm(y[ok] ~ Z[ok, , drop=FALSE]))
  ct <- cor.test(rx, ry)
  c(r=unname(ct$estimate), p=ct$p.value, n=sum(ok))
}
cps <- df$cps_a
targets <- list(c("VATPase_n","Neuron V-ATPase"),
                c("LAMP1_n","Neuron LAMP1"),
                c("LDHB_n","Neuron LDHB"))
res <- data.frame()
for(tg in targets){
  y <- df[[tg[1]]]
  zA <- cbind(cps)
  zB <- cbind(cps, df$logUMI_a, df$logUMI_n)
  zC <- cbind(cps, df$gmean_a,  df$gmean_n)            # transcriptome-wide global
  A <- pcorZ(df$MCT4, y, zA)
  B <- pcorZ(df$MCT4, y, zB)
  C <- pcorZ(df$MCT4, y, zC)
  z0 <- suppressWarnings(cor.test(df$MCT4, y))
  res <- rbind(res, data.frame(
    target = tg[2],
    zero_r = round(unname(z0$estimate),3),
    A_CPS_r = round(A["r"],3),       A_CPS_p = signif(A["p"],3),
    B_logUMI_r = round(B["r"],3),    B_logUMI_p = signif(B["p"],3),
    C_GLOBAL_r = round(C["r"],3),    C_GLOBAL_p = signif(C["p"],3),
    n = A["n"]
  ))
}
cat("\n================= results =================\n")
print(res, row.names=FALSE)
write.csv(res, file.path(OUTDIR, "crosscell_global_results.csv"), row.names=FALSE)
cat(sprintf("\nsaved: %s\n", file.path(OUTDIR, "crosscell_global_results.csv")))
cat("\nInterpretation:\n")
cat("  C_GLOBAL_r staying >= +0.3 with p<0.05 -> coupling robust (not a global-expression artefact).\n")
cat("  C_GLOBAL_r collapsing toward 0 -> coupling largely explained by donor global expression.\n")
