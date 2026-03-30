# =============================================================================
# 01_extract_seaad.R
#
# Extract gene expression from SEA-AD h5ad and export per-cell and
# bin-level summary CSVs used by all figure scripts.
#
# Input:  SEA-AD h5ad file (requires data access agreement)
#         → https://portal.brain-map.org
#
# Output (written to data/processed/):
#   astro_comprehensive.csv        — per-cell astrocyte expression
#   neuron_comprehensive.csv       — per-cell excitatory neuron expression
#   donor_level_summary.csv        — donor-level (n=84) ANLS/V-ATPase/MCT4
#   astro_subtype_trajectories.csv — per-supertype ANLS bin trajectories
#   astro_bin_means.csv            — bin-level means (used by figure scripts)
#   neuron_bin_means.csv           — bin-level means
# =============================================================================

suppressPackageStartupMessages({
  library(rhdf5)
  library(dplyr)
  library(data.table)
})

set.seed(42)

# ── Paths (edit these) ────────────────────────────────────────────────────────
h5ad_path <- "path/to/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad"
out_path  <- "data/processed/"
dir.create(out_path, recursive = TRUE, showWarnings = FALSE)

# ── Cell type classification ──────────────────────────────────────────────────
classify_cell <- function(label) {
  label <- as.character(label)
  if (grepl("^Astro", label, ignore.case = TRUE)) return("Astrocyte")
  if (grepl("^Micro", label, ignore.case = TRUE)) return("Microglia")
  if (grepl("^Endo",  label, ignore.case = TRUE)) return("Endothelial")
  if (grepl("Sst",   label, ignore.case = TRUE)) return("SST_Inhibitory")
  if (grepl("^L[0-9]|IT$|ET$|CT$|NP$|L6b", label)) return("Excitatory_Neuron")
  if (grepl("^Oligo", label, ignore.case = TRUE)) return("Oligodendrocyte")
  return("Other")
}

# ── Target gene list ──────────────────────────────────────────────────────────
target_genes <- unique(c(
  # ANLS / Energy substrate
  "SLC2A1",  # GLUT1
  "SLC2A3",  # GLUT3
  "SLC16A1", # MCT1
  "SLC16A3", # MCT4 ← primary ANLS hub
  "SLC16A7", # MCT2 (neuronal importer)
  "LDHA", "LDHB", "HK1", "HK2", "PKM",
  "PFKFB3", "GPI", "ENO1", "ENO2", "GAPDH", "PDK1",

  # V-ATPase subunits (10)
  "ATP6V1A", "ATP6V1B2", "ATP6V1C1", "ATP6V1E1", "ATP6V1H",
  "ATP6V0A1", "ATP6V0B", "ATP6V0C", "ATP6V0D1", "ATP6V0E1",

  # Lysosomal / autophagy
  "LAMP1", "LAMP2", "CTSB", "CTSD", "LIPA", "TFEB",
  "BECN1", "ATG5", "ATG7", "LC3B",

  # Iron pathway
  "TFRC", "FTH1", "FTL", "CP", "SLC40A1",

  # Mitochondria / mTOR
  "VDAC1", "NDUFS1",
  "MTOR", "LAMTOR1", "LAMTOR2", "LAMTOR3", "LAMTOR4", "LAMTOR5",

  # pH / signaling
  "SLC9A6",  # NHE6

  # Inflammation / checkpoint
  "PTGDS", "LCN2", "HMOX1", "SOX9", "STAT3",

  # Astrocyte metabolic overload markers
  "SLC1A2",  # EAAT2
  "ATP1A2",

  # Neuronal trophic (Paper 3)
  "NGFR", "BDNF", "NTRK1", "NTRK2",

  # Tau
  "MAPT"
))

# ── Extract from h5ad ─────────────────────────────────────────────────────────
message(">>> Reading gene index...")
gene_names_all <- h5read(h5ad_path, "var/_index")
gene_idx_map   <- match(target_genes, gene_names_all) - 1L
valid           <- !is.na(gene_idx_map)
target_genes    <- target_genes[valid]
gene_idx_map    <- gene_idx_map[valid]
message(sprintf("Found %d / %d target genes", sum(valid), length(valid)))

message(">>> Reading metadata...")
cps_raw      <- h5read(h5ad_path, "obs/Continuous Pseudo-progression Score")
sc_cat       <- h5read(h5ad_path, "obs/__categories/Subclass")
sc_idx       <- h5read(h5ad_path, "obs/Subclass") + 1L
sc_vec       <- sc_cat[sc_idx]
cell_type    <- sapply(sc_vec, classify_cell)
cps_bin      <- round(cps_raw, 1)

# Supertype for astrocyte subtypes
st_cat <- tryCatch(h5read(h5ad_path, "obs/__categories/Supertype"), error = function(e) NULL)
st_idx <- tryCatch(h5read(h5ad_path, "obs/Supertype") + 1L,        error = function(e) NULL)
st_vec <- if (!is.null(st_cat)) st_cat[st_idx] else rep(NA, length(cps_raw))

n_cells <- length(cps_raw)
message(sprintf("Total cells: %s", format(n_cells, big.mark = ",")))

message(">>> Extracting expression matrix (block read)...")
indptr     <- h5read(h5ad_path, "X/indptr", bit64conversion = "double")
block_size <- 200000
n_blocks   <- ceiling(n_cells / block_size)
expr_mat   <- matrix(0, nrow = n_cells, ncol = length(target_genes))
colnames(expr_mat) <- target_genes

for (b in seq_len(n_blocks)) {
  sc   <- (b - 1) * block_size + 1
  ec   <- min(b * block_size, n_cells)
  sp   <- indptr[sc] + 1
  ep   <- indptr[ec + 1]
  nr   <- ep - sp
  if (nr > 0) {
    idx <- h5read(h5ad_path, "X/indices", start = sp, count = nr)
    dat <- h5read(h5ad_path, "X/data",    start = sp, count = nr)
    nnz <- diff(indptr[sc:(ec + 1)])
    off <- 1
    for (i in seq_along(nnz)) {
      if (nnz[i] > 0) {
        ir <- off:(off + nnz[i] - 1)
        gh <- which(idx[ir] %in% gene_idx_map)
        if (length(gh) > 0)
          for (g in gh) {
            gc <- which(gene_idx_map == idx[ir[g]])
            expr_mat[sc + i - 1, gc] <- log1p(dat[ir[g]])
          }
        off <- off + nnz[i]
      }
    }
  }
  if (b %% 5 == 0) message(sprintf("  block %d / %d", b, n_blocks))
}

# ── Assemble data.table ───────────────────────────────────────────────────────
meta <- data.table(
  cps       = cps_raw,
  bin       = cps_bin,
  cell_type = cell_type,
  subclass  = sc_vec,
  supertype = st_vec
)
for (g in target_genes) meta[[g]] <- expr_mat[, g]

# ── Export per-cell CSVs ──────────────────────────────────────────────────────
message(">>> Saving per-cell CSVs...")

# Astrocytes
astro <- meta[cell_type == "Astrocyte"]
fwrite(astro, file.path(out_path, "astro_comprehensive.csv"))
message(sprintf("  astro: %s cells", format(nrow(astro), big.mark = ",")))

# Excitatory neurons
exc_labels <- c("L2/3 IT","L4 IT","L5 IT","L5 ET","L5/6 NP",
                "L6 IT","L6 IT Car3","L6 CT","L6b")
neuron <- meta[subclass %in% exc_labels]
fwrite(neuron, file.path(out_path, "neuron_comprehensive.csv"))
message(sprintf("  neuron: %s cells", format(nrow(neuron), big.mark = ",")))

# ── Bin-level means ───────────────────────────────────────────────────────────
message(">>> Computing bin-level means...")
bins_use <- seq(0.2, 0.9, 0.1)

astro_bin <- astro[bin %in% bins_use,
                   lapply(.SD, mean, na.rm = TRUE),
                   by = bin, .SDcols = target_genes][order(bin)]
fwrite(astro_bin, file.path(out_path, "astro_bin_means.csv"))

neuron_bin <- neuron[bin %in% bins_use,
                     lapply(.SD, mean, na.rm = TRUE),
                     by = bin, .SDcols = target_genes][order(bin)]
fwrite(neuron_bin, file.path(out_path, "neuron_bin_means.csv"))

# ── Donor-level summary ───────────────────────────────────────────────────────
message(">>> Donor-level summary...")
anls_genes    <- c("SLC2A1", "LDHA", "SLC16A1")
vatpase_genes <- c("ATP6V1A","ATP6V1B2","ATP6V0D1","ATP6V0A1",
                   "ATP6V1C1","ATP6V1E1","ATP6V1H","ATP6V0C","ATP6V0E1","ATP6V0B")
vatpase_n     <- c("ATP6V1A","ATP6V1B2","ATP6V0A1","ATP6V0C","ATP6V0D1","ATP6V1E1")

# Add donor metadata from h5ad if available
# (Braak, ABC score, cognitive status — stored in obs annotations)
donor_meta_cols <- c("Braak stage (categorized)", "Overall AD neuropathological Change",
                     "Cognitive Status", "Age at Death", "Sex")
obs_cols <- tryCatch(h5ls(h5ad_path, "/obs")$name, error = function(e) character(0))

astro[, ANLS    := rowMeans(.SD, na.rm = TRUE), .SDcols = intersect(anls_genes, names(astro))]
astro[, VATpase := rowMeans(.SD, na.rm = TRUE), .SDcols = intersect(vatpase_genes, names(astro))]
astro[, MCT4    := SLC16A3]
neuron[, VATPase_n := rowMeans(.SD, na.rm = TRUE), .SDcols = intersect(vatpase_n, names(neuron))]

# Round CPS to 2 decimals as donor proxy
astro[, donor_id := round(cps, 2)]
neuron[, donor_id := round(cps, 2)]

donor <- merge(
  astro[, .(ANLS = mean(ANLS), MCT4 = mean(MCT4), VATpase = mean(VATpase),
             mean_cps = mean(cps), n_astro = .N), by = donor_id],
  neuron[, .(VATPase_n = mean(VATPase_n, na.rm = TRUE)), by = donor_id],
  by = "donor_id", all.x = TRUE
)
fwrite(donor, file.path(out_path, "donor_level_summary.csv"))

# ── Astrocyte subtype trajectories ────────────────────────────────────────────
message(">>> Astrocyte subtype trajectories...")
if (all(!is.na(astro$supertype))) {
  astro[, ANLS := rowMeans(.SD, na.rm = TRUE), .SDcols = intersect(anls_genes, names(astro))]
  sub_traj <- astro[bin %in% bins_use,
                    .(ANLS = mean(ANLS), n_cells = .N),
                    by = .(bin, supertype)][order(bin, supertype)]
  setnames(sub_traj, "supertype", "subtype")
  fwrite(sub_traj, file.path(out_path, "astro_subtype_trajectories.csv"))
}

message("\n>>> Extraction complete. Output: ", out_path)
