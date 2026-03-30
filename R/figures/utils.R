# =============================================================================
# utils.R — Shared helpers for all figure scripts
#
# Source this at the top of each figure script:
#   source("R/figures/utils.R")
#
# IMPORTANT — path configuration:
#   All scripts use relative paths from the repository root.
#   Run scripts with working directory set to the repo root, e.g.:
#     setwd("D:/path/to/AD-Metabolic-Collapse")
#     source("R/figures/Fig1_Dissociation.R")
#
# For Fig5 and Fig6 (ADNI proteomics), set paths in those scripts directly.
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(scales)
})

# ── Paths (relative to repo root) ─────────────────────────────────────────────
# Bin-level means — included in repository, sufficient for Fig1-4 and ED figs
DATA_BIN <- "data/sample/"

# Figure output — created automatically if it does not exist
FIG_OUT <- "output/figures/"
dir.create(FIG_OUT, recursive = TRUE, showWarnings = FALSE)

# NOTE: Fig5 and Fig6 require ADNI data (not included).
# Set EMORY_PATH and ADNIMERGE_PATH at the top of those scripts.

# ── Theme ─────────────────────────────────────────────────────────────────────
theme_paper <- theme_bw(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 13, hjust = 0),
    axis.title       = element_text(size = 11),
    axis.text        = element_text(size = 9),
    legend.text      = element_text(size = 8),
    legend.title     = element_text(size = 9, face = "bold"),
    legend.position  = "bottom",
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "#EBF5FB"),
    strip.text       = element_text(face = "bold")
  )

dx_colors <- c("CN" = "#4DAF4A", "MCI" = "#FF7F00", "DEM" = "#E41A1C")

# ── Gene composite definitions ────────────────────────────────────────────────
ANLS_GENES    <- c("SLC2A1", "LDHA", "SLC16A1")
VATPASE_GENES <- c("ATP6V1A","ATP6V1B2","ATP6V0D1","ATP6V0A1",
                   "ATP6V1C1","ATP6V1E1","ATP6V1H","ATP6V0C","ATP6V0E1","ATP6V0B")
VATPASE_N     <- c("ATP6V1A","ATP6V1B2","ATP6V0A1","ATP6V0C","ATP6V0D1","ATP6V1E1")

add_composites_astro <- function(df) {
  df$ANLS      <- rowMeans(df[, intersect(ANLS_GENES,    names(df)), drop=FALSE], na.rm=TRUE)
  df$VATpase   <- rowMeans(df[, intersect(VATPASE_GENES, names(df)), drop=FALSE], na.rm=TRUE)
  df$MCT4      <- df$SLC16A3
  ed_genes     <- c("SLC2A1","LDHA","SLC16A1","PKM","HK1")
  df$ED_energy <- rowMeans(df[, intersect(ed_genes, names(df)), drop=FALSE], na.rm=TRUE)
  df$ED_ratio  <- df$ED_energy / df$VATpase
  df
}

add_composites_neuron <- function(df) {
  df$VATPase <- rowMeans(df[, intersect(VATPASE_N, names(df)), drop=FALSE], na.rm=TRUE)
  df
}

# ── ADNI data loader (used by Fig5 and Fig6) ──────────────────────────────────
#
# Usage:
#   em <- load_adni_proteomics(
#     emory_path     = "D:/work/emory_results/",
#     adnimerge_path = "D:/work/ADNIMERGE2/ADNIMERGE2/data/"
#   )
#
load_adni_proteomics <- function(emory_path, adnimerge_path) {

  # 1. Emory TMT-MS proteomics ------------------------------------------------
  # Searches for CSV or RDA/RDS in the given folder.
  # Typical ADNI filename: EmoryCSFProteomics_TMT_*.csv
  csv_files <- list.files(emory_path,
                          pattern = "(?i)(emory|proteom|tmt).*\\.csv$",
                          full.names = TRUE)
  rda_files <- list.files(emory_path,
                          pattern = "(?i)(emory|proteom|tmt).*\\.(rda|rds)$",
                          full.names = TRUE)

  if (length(csv_files) > 0) {
    message("Loading Emory CSV: ", basename(csv_files[1]))
    prot <- read.csv(csv_files[1], stringsAsFactors = FALSE)
  } else if (length(rda_files) > 0) {
    message("Loading Emory RDA/RDS: ", basename(rda_files[1]))
    if (grepl("\\.rds$", rda_files[1], ignore.case = TRUE)) {
      prot <- readRDS(rda_files[1])
    } else {
      env <- new.env()
      load(rda_files[1], envir = env)
      prot <- get(ls(env)[1], envir = env)
    }
  } else {
    stop(
      "Emory TMT-MS file not found in: ", emory_path, "\n",
      "Download from https://adni.loni.usc.edu\n",
      "Search 'Emory' or 'TMT' in the ADNI study data inventory.\n",
      "Expected pattern: *emory*, *proteom*, or *tmt* in filename."
    )
  }

  # Replace below-detection sentinel (-4 → NA)
  prot_cols <- names(prot)[grepl("_[A-Z][0-9]", names(prot))]
  for (col in prot_cols)
    if (is.numeric(prot[[col]])) prot[[col]][prot[[col]] == -4] <- NA

  # 2. DXSUM — diagnosis per visit --------------------------------------------
  dxsum_file <- file.path(adnimerge_path, "DXSUM.rda")
  if (!file.exists(dxsum_file))
    stop("DXSUM.rda not found: ", dxsum_file)
  load(dxsum_file)   # → object DXSUM

  # Map to CN / MCI / DEM
  if ("DXCHANGE" %in% names(DXSUM)) {
    dx_map <- c("1"="CN","2"="MCI","3"="MCI","4"="MCI",
                "5"="CN","6"="MCI","7"="DEM","8"="DEM")
    DXSUM$DX <- dx_map[as.character(DXSUM$DXCHANGE)]
  } else if ("DIAGNOSIS" %in% names(DXSUM)) {
    dx_map <- c("CN"="CN","MCI"="MCI","Dementia"="DEM","AD"="DEM")
    DXSUM$DX <- dx_map[DXSUM$DIAGNOSIS]
  } else {
    stop("Cannot find DXCHANGE or DIAGNOSIS column in DXSUM.rda")
  }

  # Most recent visit per subject
  DXSUM <- DXSUM[order(DXSUM$RID, DXSUM$VISCODE), ]
  DXSUM_last <- DXSUM[!duplicated(DXSUM$RID, fromLast = TRUE), ]

  # 3. ADSL — baseline demographics (AGE, SEX) --------------------------------
  adsl_file <- file.path(adnimerge_path, "ADSL.rda")
  if (!file.exists(adsl_file))
    stop("ADSL.rda not found: ", adsl_file)
  load(adsl_file)    # → object ADSL

  keep_cols <- intersect(c("RID","AGE","PTGENDER","PTEDUCAT"), names(ADSL))
  ADSL_sub  <- ADSL[!duplicated(ADSL$RID), keep_cols]

  # 4. RID harmonisation ------------------------------------------------------
  if (!"RID" %in% names(prot) && "SUBJID" %in% names(prot))
    prot$RID <- as.integer(gsub("[^0-9]", "", prot$SUBJID))
  if (!"RID" %in% names(prot))
    stop("Cannot find RID column in Emory proteomics file.")

  # 5. Merge ------------------------------------------------------------------
  merged <- merge(prot,  DXSUM_last[, c("RID","DX")], by="RID", all.x=TRUE)
  merged <- merge(merged, ADSL_sub,                    by="RID", all.x=TRUE)

  if ("PTGENDER" %in% names(merged) && !"SEX" %in% names(merged))
    merged$SEX <- merged$PTGENDER

  merged <- merged[!is.na(merged$DX) & merged$DX %in% c("CN","MCI","DEM"), ]
  merged$DX <- factor(merged$DX, levels = c("CN","MCI","DEM"))

  cat(sprintf(
    "ADNI loaded: n=%d  CN=%d  MCI=%d  DEM=%d\n",
    nrow(merged),
    sum(merged$DX=="CN"), sum(merged$DX=="MCI"), sum(merged$DX=="DEM")
  ))

  merged
}

# ── Statistics helpers ────────────────────────────────────────────────────────
partial_cor <- function(x, y, z) {
  ok <- complete.cases(x, y, z)
  if (sum(ok) < 10) return(list(r=NA, p=NA, rx=NA, ry=NA))
  rx <- residuals(lm(x[ok] ~ z[ok]))
  ry <- residuals(lm(y[ok] ~ z[ok]))
  ct <- cor.test(rx, ry)
  list(r=ct$estimate, p=ct$p.value, rx=rx, ry=ry)
}

partial_cor_double <- function(x, y, z1, z2) {
  ok <- complete.cases(x, y, z1, z2)
  if (sum(ok) < 10) return(list(r=NA, p=NA))
  rx <- residuals(lm(x[ok] ~ z1[ok] + z2[ok]))
  ry <- residuals(lm(y[ok] ~ z1[ok] + z2[ok]))
  ct <- cor.test(rx, ry)
  list(r=ct$estimate, p=ct$p.value)
}

partial_cor_triple <- function(x, y, z1, z2, z3) {
  ok <- complete.cases(x, y, z1, z2, z3)
  if (sum(ok) < 10) return(list(r=NA, p=NA))
  rx <- residuals(lm(x[ok] ~ z1[ok] + z2[ok] + z3[ok]))
  ry <- residuals(lm(y[ok] ~ z1[ok] + z2[ok] + z3[ok]))
  ct <- cor.test(rx, ry)
  list(r=ct$estimate, p=ct$p.value)
}

norm_base <- function(vals, bins) vals / vals[which.min(bins)]

get_slope <- function(df, col, bins) {
  d    <- df[df$bin %in% bins, ]
  base <- d[[col]][d$bin == min(d$bin)]
  if (is.na(base) || base == 0) return(list(beta=NA, se=NA, p=NA))
  d$norm <- d[[col]] / base
  fit    <- lm(norm ~ bin, data=d)
  s      <- summary(fit)$coefficients
  list(beta=s[2,1], se=s[2,2], p=s[2,4])
}

sig_star <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else "n.s."
}

# ── Save helper ───────────────────────────────────────────────────────────────
save_fig <- function(plot, filename, width=14, height=10) {
  path <- file.path(FIG_OUT, filename)
  ggsave(path, plot, width=width, height=height, dpi=300, bg="white")
  message("Saved: ", path)
}
