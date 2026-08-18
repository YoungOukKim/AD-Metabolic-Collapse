# =============================================================================
# 01_extract_seaad.R  [CORRECTED: donor-level aggregation now by real Donor ID]
#
# Extract gene expression from SEA-AD h5ad and export per-cell and
# bin-level summary CSVs used by all figure scripts.
#
# Input:  SEA-AD h5ad file (requires data access agreement)
#         -> https://portal.brain-map.org
#
# Output (written to data/processed/):
#   astro_comprehensive.csv         -  per-cell astrocyte expression
#   neuron_comprehensive.csv        -  per-cell excitatory neuron expression
#   donor_level_summary.csv         -  donor-level (n=84) ANLS/V-ATPase/MCT4
#   astro_subtype_trajectories.csv  -  per-supertype ANLS bin trajectories
#   astro_bin_means.csv             -  bin-level means (used by figure scripts)
#   neuron_bin_means.csv            -  bin-level means
# =============================================================================

# Paths are RELATIVE to the repository root. Set the working directory there,
# or override with the environment variables SEAAD_H5AD / ROSMAP_ASTRO /
# ROSMAP_CLIN / P2_OUT_DIR. Raw data are not redistributable; see README.md.

suppressPackageStartupMessages({
  library(rhdf5)
  library(dplyr)
  library(data.table)
})

set.seed(42)

# -- Paths (edit these) --------------------------------------------------------
h5ad_path <- "path/to/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad"
out_path  <- "data/processed/"
dir.create(out_path, recursive = TRUE, showWarnings = FALSE)

# -- Cell type classification --------------------------------------------------
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

# -- Target gene list ----------------------------------------------------------
target_genes <- unique(c(
  # ANLS / Energy substrate
  "SLC2A1",  # GLUT1
  "SLC2A3",  # GLUT3
  "SLC16A1", # MCT1
  "SLC16A3", # MCT4 <- primary ANLS hub
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

# -- Extract from h5ad ---------------------------------------------------------
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

message(">>> Fast block read (astrocytes + excitatory neurons only) ...")
# Only astrocytes and excitatory neurons are used by any downstream output
# (per-cell CSVs, bin means, donor summary), so we read just those cells and
# fill via a single match() per cell instead of a nested per-gene which().
# Outputs are identical to the previous version; this is purely faster.
indptr     <- h5read(h5ad_path, "X/indptr", bit64conversion = "double")

# -- Storage-convention guard (added 2026-07) ----------------------------------
# SEA-AD stores X already normalised (CP10k + log1p). An earlier version of this
# script applied log1p unconditionally at the per-cell fill below, which logs an
# already-logged matrix. 02_global_expression_sensitivity.R tests before logging;
# this script now does the same, so both extraction paths use one convention.
.samp     <- h5read(h5ad_path, "X/data", start = 1,
                    count = min(2e5, indptr[length(indptr)]))
is_counts <- all(abs(.samp - round(.samp)) < 1e-8)
message(sprintf("  X storage: %s",
                if (is_counts) "raw counts -> CP10k + log1p applied"
                else           "already normalised -> used as stored"))
exc_labels <- c("L2/3 IT","L4 IT","L5 IT","L5 ET","L5/6 NP",
                "L6 IT","L6 IT Car3","L6 CT","L6b")
keep_mask  <- (cell_type == "Astrocyte") | (sc_vec %in% exc_labels)
block_size <- 100000
n_blocks   <- ceiling(n_cells / block_size)
expr_mat   <- matrix(0, nrow = n_cells, ncol = length(target_genes))
colnames(expr_mat) <- target_genes
message(sprintf("  target cells to read: %s of %s",
                format(sum(keep_mask), big.mark = ","), format(n_cells, big.mark = ",")))

for (blk in seq_len(n_blocks)) {
  s0 <- (blk - 1L) * block_size + 1L
  e0 <- min(blk * block_size, n_cells)
  sp <- indptr[s0]; cnt <- indptr[e0 + 1L] - sp
  if (cnt > 0) {
    ci <- h5read(h5ad_path, "X/indices", start = sp + 1L, count = cnt, bit64conversion = "double")
    cd <- h5read(h5ad_path, "X/data",    start = sp + 1L, count = cnt)
    for (k in which(keep_mask[s0:e0])) {
      g  <- s0 + k - 1L
      a  <- indptr[g] - sp + 1L
      bb <- indptr[g + 1L] - sp
      if (bb < a) next
      h  <- match(gene_idx_map, ci[a:bb]); ok <- !is.na(h)
      if (any(ok)) {
        vals <- cd[a:bb][h[ok]]
        expr_mat[g, ok] <- if (is_counts)
          log1p(vals / max(sum(cd[a:bb]), 1) * 1e4) else vals
      }
    }
  }
  cat(sprintf("  block %d / %d  (%d cells)\r", blk, n_blocks, e0)); flush.console()
}
cat("\n")

# -- Assemble data.table -------------------------------------------------------
meta <- data.table(
  cps       = cps_raw,
  bin       = cps_bin,
  cell_type = cell_type,
  subclass  = sc_vec,
  supertype = st_vec
)
for (g in target_genes) meta[[g]] <- expr_mat[, g]

# -- Attach donor + clinical metadata, then drop SEA-AD reference donors -------
# Reference (neurotypical) donors carry no Continuous Pseudo-progression Score;
# they appear as 3 extra donors with NA CPS/Braak. Excluding them gives the
# 84-donor AD pseudo-progression cohort used throughout (astro -> 67,419 cells).
donor_cat <- h5read(h5ad_path, "obs/__categories/Donor ID")
donor_vec <- donor_cat[h5read(h5ad_path, "obs/Donor ID") + 1L]
get_obs_cat <- function(field) {
  cats <- tryCatch(h5read(h5ad_path, paste0("obs/__categories/", field)), error = function(e) NULL)
  if (is.null(cats)) return(rep(NA_character_, length(donor_vec)))
  cats[h5read(h5ad_path, paste0("obs/", field)) + 1L]
}
read_first <- function(cands) {          # robust to obs field-name variants (esp. Braak)
  for (f in cands) { v <- get_obs_cat(f); if (!all(is.na(v))) { message(sprintf("  Braak field matched: %s", f)); return(v) } }
  message("  WARNING: no Braak field matched -> braak_num will be NA (check obs field name)")
  rep(NA_character_, length(donor_vec))
}
meta[, donor := donor_vec]
meta[, braak := read_first(c("Braak", "Braak stage", "Braak Stage", "Braak stage (categorized)"))]
meta[, cerad := get_obs_cat("CERAD score")]
meta[, abc   := get_obs_cat("Overall AD neuropathological Change")]
meta[, cog   := get_obs_cat("Cognitive Status")]
meta <- meta[!is.na(cps)]                # reference donors have no CPS -> excluded (87 -> 84 donors)
message(sprintf("  cells after dropping reference donors (no CPS): %s", format(nrow(meta), big.mark = ",")))

# -- Export per-cell CSVs ------------------------------------------------------
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

# -- Bin-level means -----------------------------------------------------------
message(">>> Computing bin-level means...")
bins_use <- round(seq(0.2, 0.9, 0.1), 1)   # round to avoid float %in% mismatch (was dropping 0.3, 0.6)

astro_bin <- astro[bin %in% bins_use,
                   lapply(.SD, mean, na.rm = TRUE),
                   by = bin, .SDcols = target_genes][order(bin)]
fwrite(astro_bin, file.path(out_path, "astro_bin_means.csv"))

neuron_bin <- neuron[bin %in% bins_use,
                     lapply(.SD, mean, na.rm = TRUE),
                     by = bin, .SDcols = target_genes][order(bin)]
fwrite(neuron_bin, file.path(out_path, "neuron_bin_means.csv"))

# -- Donor-level summary -------------------------------------------------------
message(">>> Donor-level summary...")
anls_genes    <- c("SLC2A1", "LDHA", "SLC16A1")
vatpase_genes <- c("ATP6V1A","ATP6V1B2","ATP6V0D1","ATP6V0A1",
                   "ATP6V1C1","ATP6V1E1","ATP6V1H","ATP6V0C","ATP6V0E1","ATP6V0B")
vatpase_n     <- c("ATP6V1A","ATP6V1B2","ATP6V0A1","ATP6V0C","ATP6V0D1","ATP6V1E1")

# Add donor metadata from h5ad if available
# (Braak, ABC score, cognitive status  -  stored in obs annotations)
donor_meta_cols <- c("Braak stage (categorized)", "Overall AD neuropathological Change",
                     "Cognitive Status", "Age at Death", "Sex")
obs_cols <- tryCatch(h5ls(h5ad_path, "/obs")$name, error = function(e) character(0))

astro[, ANLS    := rowMeans(.SD, na.rm = TRUE), .SDcols = intersect(anls_genes, names(astro))]
astro[, VATpase := rowMeans(.SD, na.rm = TRUE), .SDcols = intersect(vatpase_genes, names(astro))]
astro[, MCT4    := SLC16A3]
neuron[, VATPase_n := rowMeans(.SD, na.rm = TRUE), .SDcols = intersect(vatpase_n, names(neuron))]

# ============================================================================
# CORRECTED donor-level aggregation.
# Previous version used `donor_id := round(cps, 2)`, which collapsed cells into
# ~9 CPS bins and broadcast bin means onto 84 donor rows (pseudo-replication:
# the real independent unit was ~9 bins, not 84 donors). This inflated both the
# effect size and the effective sample size, and made the "CPS-adjusted partial
# correlation" vacuous (the aggregation key was itself ~CPS). We now aggregate
# by the REAL "Donor ID" from obs, with each donor contributing its own value.
# ============================================================================
# (donor + clinical metadata attached above, before reference-donor filtering)

astro  <- meta[cell_type == "Astrocyte"]
neuron <- meta[subclass %in% exc_labels]
astro[,  ANLS      := rowMeans(.SD, na.rm = TRUE), .SDcols = intersect(anls_genes,    names(astro))]
astro[,  VATpase   := rowMeans(.SD, na.rm = TRUE), .SDcols = intersect(vatpase_genes, names(astro))]
astro[,  MCT4      := SLC16A3]
neuron[, VATPase_n := rowMeans(.SD, na.rm = TRUE), .SDcols = intersect(vatpase_n, names(neuron))]
ed_genes <- c("SLC2A1","LDHA","SLC16A1","PKM","HK1")
astro[,  ED_energy := rowMeans(.SD, na.rm = TRUE), .SDcols = intersect(ed_genes, names(astro))]

mode1 <- function(x){ x <- x[!is.na(x)]; if(!length(x)) return(NA_character_); names(sort(table(x), decreasing=TRUE))[1] }
donor_astro <- astro[, .(ANLS = mean(ANLS), MCT4 = mean(MCT4), VATpase = mean(VATpase),
                         ED_energy = mean(ED_energy), TFRC_val = mean(TFRC),
                         mean_cps = mean(cps), n_astro = .N,
                         braak = mode1(braak), cerad = mode1(cerad),
                         abc = mode1(abc), cognitive = mode1(cog)), by = donor]
donor_neuron <- neuron[, .(VATPase_n = mean(VATPase_n, na.rm = TRUE),
                           LDHB_n  = mean(LDHB,  na.rm = TRUE),
                           LAMP1_n = mean(LAMP1, na.rm = TRUE),
                           n_exc   = .N), by = donor]
donor <- merge(donor_astro, donor_neuron, by = "donor", all.x = TRUE)
donor <- donor[n_astro >= 20 & !is.na(n_exc) & n_exc >= 20]   # min cells per class
setnames(donor, "donor", "donor_id")
# ordinal numeric codes (correct severity order)
braak_map <- c("Braak 0"=0,"Braak I"=1,"Braak II"=2,"Braak III"=3,"Braak IV"=4,"Braak V"=5,"Braak VI"=6)
cerad_map <- c("Absent"=0,"Sparse"=1,"Moderate"=2,"Frequent"=3)
abc_map   <- c("Not AD"=0,"Low"=1,"Intermediate"=2,"High"=3)
donor[, braak_num := braak_map[braak]]
donor[, cerad_num := cerad_map[cerad]]
donor[, abc_num   := abc_map[abc]]
donor[, ED_ratio  := ED_energy / VATPase_n]
fwrite(donor, file.path(out_path, "donor_level_summary.csv"))
donor[, abc_score := abc]; donor[, cognitive := cognitive]   # aliases for downstream scripts
fwrite(donor, file.path(out_path, "donor_level_summary.csv"))
pc <- function(x,y,z){ok<-complete.cases(x,y,z); r0<-cor(x[ok],y[ok]); rr<-cor(residuals(lm(x[ok]~z[ok])),residuals(lm(y[ok]~z[ok]))); c(r0,rr)}
for (yy in c("VATPase_n","LDHB_n","LAMP1_n")) if (yy %in% names(donor)) {
  v<-pc(donor$MCT4, donor[[yy]], donor$mean_cps)
  message(sprintf("   MCT4 ~ %s (true donor): zero r=%+.3f, CPS-partial r=%+.3f", yy, v[1], v[2]))
}
message(sprintf(">>> donor_level_summary.csv written: %d TRUE donors (by Donor ID, no bin-broadcast)", nrow(donor)))

# -- Astrocyte subtype trajectories --------------------------------------------
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
