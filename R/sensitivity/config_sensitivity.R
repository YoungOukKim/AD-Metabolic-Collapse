# ---------------------------------------------------------------------
# config.R - the only file that needs editing before a run.
# All paths are relative to the repository root.
# ---------------------------------------------------------------------

# Place (or symlink) the SEA-AD MTG h5ad here.
# Download: https://portal.brain-map.org/explore/seattle-alzheimers-disease
H5AD    <- file.path("data", "SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad")

OUT_DIR <- file.path("output", "sensitivity")

# Cell block size for the single-pass reader. 2000-5000 is a safe range;
# larger blocks read fewer times but hold several GB in memory per block.
BLOCK_SIZE <- 2000L

# Accumulate the full donor x gene matrix (needed for 03_specificity_screen.R).
# Set FALSE for a faster run that only produces the donor-level table.
FULL_GENE_MATRIX <- TRUE

SEED <- 20260720
