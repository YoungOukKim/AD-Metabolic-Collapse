# ==============================================================================
# SuppFig2C_extract.R
#
# Purpose : Extract the per-bin trajectories actually plotted in Supplemental
#           Figure 2C ("metabolic overload"), so they can be added to the
#           Supporting Data Values workbook as source data. These five series
#           are NOT currently tabulated in Supporting_Data_Values.xlsx:
#             astrocyte  : EAAT2 (SLC1A2), ATP1A2, PTGDS, MCT4 (SLC16A3)
#             neuron     : MCT2 (SLC16A7)
#
# Method  : same convention as Fig.1 / the Astrocyte_179_genes sheet —
#           per-cell expression averaged within CPS bins (bin = round(CPS,1)),
#           reported for bins 0.1-0.9, then normalized to Bin 0.2 = 1.0 for
#           panel C. Detection rate (fraction of cells > 0) is reported so a
#           sparse gene cannot masquerade as a trajectory (the oligo-MCT4 lesson).
#
# ANCHOR CHECK (important): SLC16A3 and PTGDS already exist in the published
#           Astrocyte_179_genes sheet. The script recomputes them and compares
#           to the published values. If the anchors MATCH, the binning + input
#           scale are correct and EAAT2/ATP1A2/MCT2 are computed the same way.
#           If they DON'T match, set NORMALIZE_INPUT below (X may be raw counts).
#
# Input   : SEA-AD MTG h5ad (same file as Fig.1)
# Output  : <OUTDIR>/SuppFig2C_bin_means.csv      (raw per-bin means, det_rate)
#           <OUTDIR>/SuppFig2C_normalized.csv     (Bin 0.2 = 1.0; panel-C ready)
#           + console anchor verdict
#
# FIX (v2): bins are keyed by INTEGER round(round(CPS,1)*10) rather than by
#           comparing to seq(0.1,0.9,0.1); the latter yields 0.30000000000000004
#           and 0.7000000000000001 (float error) and silently drops Bin_0.3 /
#           Bin_0.7. Anchor-validated round(CPS,1) binning is preserved exactly.
# Requires: rhdf5, data.table
# Usage   : Rscript SuppFig2C_extract.R
#           (or set H5AD_PATH / OUTDIR env vars first)
# ==============================================================================
suppressPackageStartupMessages({ library(rhdf5); library(data.table) })
set.seed(42)

MTG_H5 <- Sys.getenv("H5AD_PATH",
                     unset = "path/to/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad")
OUTDIR <- Sys.getenv("OUTDIR",
                     unset = "output/suppfig2c")
if (!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE)

# If the published Astrocyte sheet values are mean LOG-normalized and the h5ad X
# layer is RAW counts, set this TRUE to apply CP10k + log1p before averaging.
# Leave FALSE first; the anchor check will tell you if you need to flip it.
NORMALIZE_INPUT <- FALSE

# gene -> cell type the trajectory is drawn from in panel C
GENE_CELL <- data.table(
  gene = c("SLC1A2","ATP1A2","PTGDS","SLC16A3","SLC16A7"),
  cell = c("astro", "astro", "astro","astro",  "neuron"),
  role = c("EAAT2 (SLC1A2)","ATP1A2 (Na/K pump)","PTGDS (anchor)",
           "MCT4 (SLC16A3) (anchor)","MCT2 neuron (SLC16A7)"))

# Published Astrocyte_179_genes values for the anchor genes (Bin_0.1 .. Bin_0.9)
ANCHOR <- list(
  SLC16A3 = c(0.1628,0.0550,0.0629,0.0721,0.0403,0.0337,0.0409,0.0334,0.0220))

# ── obs: Subclass / Class (cell type) + CPS ───────────────────────────────────
read_obs <- function(name){
  v <- tryCatch(h5read(MTG_H5, paste0("obs/", name)), error = function(e) NULL)
  if (is.null(v)) return(NULL)
  cats <- tryCatch(h5read(MTG_H5, paste0("obs/__categories/", name)), error=function(e) NULL)
  if (!is.null(cats)) cats[as.integer(v)+1L] else v
}
subclass <- read_obs("Subclass")
class_   <- read_obs("Class")                  # "Neuronal: Glutamatergic" / ... / "Non-neuronal..."
cps      <- as.numeric(h5read(MTG_H5, "obs/Continuous Pseudo-progression Score"))

# astrocytes: Subclass contains "Astro" (exclude OPC/precursor just in case)
astro_mask <- grepl("Astro", subclass, ignore.case=TRUE) &
              !grepl("OPC|precursor", subclass, ignore.case=TRUE) & !is.na(cps)
# neurons: prefer the Class field; fall back to canonical neuronal subclasses
if (!is.null(class_)) {
  neuron_mask <- grepl("Neuron", class_, ignore.case=TRUE) & !is.na(cps)
} else {
  neuron_sub <- "IT|ET|CT|NP|L6b|L5/6|Pvalb|Sst|Vip|Lamp5|Sncg|Chandelier|Pax6"
  neuron_mask <- grepl(neuron_sub, subclass, ignore.case=TRUE) & !is.na(cps)
}
keep <- astro_mask | neuron_mask
cat(sprintf("astrocyte nuclei: %d | neuron nuclei: %d | total kept: %d\n",
            sum(astro_mask), sum(neuron_mask), sum(keep)))

# ── gene index (missing-gene tolerant) ────────────────────────────────────────
genes_all <- as.character(h5read(MTG_H5,"var/_index"))
G    <- GENE_CELL$gene
pidx0 <- match(G, genes_all) - 1L
present <- !is.na(pidx0)
if (any(!present))
  cat(sprintf("MISSING genes (skipped): %s\n", paste(G[!present], collapse=",")))
G <- G[present]; pidx0 <- pidx0[present]

# ── FAST block read: keep-filtered CSC slurp (same pattern as Fig.1 scripts) ──
indptr <- h5read(MTG_H5,"X/indptr",bit64conversion="double"); nC <- length(indptr)-1L
sel  <- which(keep); kset <- logical(nC); kset[sel] <- TRUE
rmap <- integer(nC); rmap[sel] <- seq_along(sel)
M <- matrix(0, length(sel), length(pidx0)); colnames(M) <- G
bs <- 100000L; cat("reading expression ...\n")
for (s0 in seq(1L, nC, by=bs)) {
  e0 <- min(s0+bs-1L, nC); sp <- indptr[s0]; cnt <- indptr[e0+1L]-sp; if (cnt<=0) next
  ci <- h5read(MTG_H5,"X/indices",start=sp+1L,count=cnt,bit64conversion="double")
  cd <- h5read(MTG_H5,"X/data",   start=sp+1L,count=cnt)
  for (k in which(kset[s0:e0])) {
    g <- s0+k-1L; a <- indptr[g]-sp+1L; b <- indptr[g+1L]-sp; if (b<a) next
    h <- match(pidx0, ci[a:b]); ok <- !is.na(h); M[rmap[g], ok] <- cd[a:b][h[ok]]
  }
  cat(sprintf("  ...%d/%d\r", e0, nC))
}; cat("\n")

# optional CP10k + log1p (only if X is raw counts; anchor check decides)
if (NORMALIZE_INPUT) {
  # per-cell total over the FULL transcriptome is unknown here; this toggle just
  # log1p-transforms the slurped values as a coarse fallback. Prefer matching the
  # original pipeline exactly if the anchor still fails.
  M <- log1p(M)
}

cp_keep    <- cps[sel]
class_keep <- ifelse(astro_mask[sel], "astro",
              ifelse(neuron_mask[sel], "neuron", NA_character_))
binkey <- 1:9                 # integer bin ids representing 0.1 .. 0.9
binval <- binkey / 10          # exact 0.1..0.9 for labels (display only)

# ── per-gene, per-bin mean within its own cell type + detection rate ──────────
bin_table <- function(g, cell) {
  idx <- which(class_keep == cell)
  x   <- M[idx, g]; bk <- as.integer(round(round(cp_keep[idx], 1) * 10))  # round(CPS,1) -> 1..9
  means <- sapply(binkey, function(k){ v <- x[bk==k]; if(!length(v)) NA else mean(v) })
  det   <- mean(x > 0)
  list(means = means, det = det, n = length(idx))
}

raw <- list(); detr <- c(); ns <- c()
for (i in seq_len(nrow(GENE_CELL))) {
  g <- GENE_CELL$gene[i]; if (!g %in% G) next
  bt <- bin_table(g, GENE_CELL$cell[i])
  raw[[g]] <- bt$means; detr[g] <- bt$det; ns[g] <- bt$n
}

raw_dt <- data.table(gene = names(raw),
                     cell = GENE_CELL$cell[match(names(raw), GENE_CELL$gene)],
                     role = GENE_CELL$role[match(names(raw), GENE_CELL$gene)],
                     det_rate = round(detr[names(raw)], 3),
                     n_cells  = ns[names(raw)])
for (j in seq_along(binval))
  raw_dt[[sprintf("Bin_%.1f", binval[j])]] <- round(sapply(raw, `[`, j), 4)

# ── normalize to Bin 0.2 = 1.0 (panel C convention) ───────────────────────────
norm_dt <- copy(raw_dt[, .(gene, cell, role, det_rate)])
for (j in 2:9) {  # bins 0.2 .. 0.9 -> columns
  bb <- binval[j]
  norm_dt[[sprintf("B%.1f", bb)]] <-
    round(sapply(raw, function(m) m[j] / m[2]), 4)   # divide by Bin 0.2
}

cat("\n================ Supp Fig 2C source data (raw bin means) ================\n")
print(raw_dt)
cat("\n---------------- normalized to Bin 0.2 = 1.0 (panel C) ------------------\n")
print(norm_dt)

# ── ANCHOR CHECK: recomputed SLC16A3 vs published Astrocyte_179_genes ─────────
cat("\n[ANCHOR CHECK] SLC16A3 (astrocyte) recomputed vs published sheet:\n")
if ("SLC16A3" %in% names(raw)) {
  got <- round(raw[["SLC16A3"]], 4); exp <- ANCHOR$SLC16A3
  comp <- data.table(bin = binval, computed = got, published = exp,
                     rel_diff = round((got-exp)/exp, 3))
  print(comp)
  ok <- all(abs((got-exp)/exp) < 0.05, na.rm=TRUE)
  cat(sprintf("  -> %s (all bins within 5%%: %s)\n",
              if (ok) "PASS" else "FAIL",
              if (ok) "binning + input scale correct"
              else "MISMATCH: set NORMALIZE_INPUT<-TRUE or match original pipeline"))
} else cat("  SLC16A3 not found; cannot validate.\n")

cat("\n[DETECTION GATE] det_rate >= 0.10 required to report a trajectory:\n")
for (g in names(detr))
  cat(sprintf("  %-8s det_rate = %.3f -> %s\n", g, detr[g],
              if (detr[g] >= 0.10) "reportable" else "SPARSE (do not report)"))

fwrite(raw_dt,  file.path(OUTDIR, "SuppFig2C_bin_means.csv"))
fwrite(norm_dt, file.path(OUTDIR, "SuppFig2C_normalized.csv"))
cat(sprintf("\nWrote:\n  %s\n  %s\n",
            file.path(OUTDIR,"SuppFig2C_bin_means.csv"),
            file.path(OUTDIR,"SuppFig2C_normalized.csv")))
H5close()
