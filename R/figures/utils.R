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
# For Fig5, Fig6 and Table2 (ADNI proteomics), set paths in those scripts.
#
# -----------------------------------------------------------------------------
# REVISION NOTE (distribution-robust reanalysis)
#   CSF TMT-MS protein abundances are strongly right-skewed; Pearson correlations
#   on untransformed abundance are dominated by a few high-abundance samples.
#   All CSF protein-protein associations are therefore computed PRIMARILY by
#   Spearman rank correlation, partial correlations by residualizing the
#   RANK-transformed variables on the control set (see spearman_partial), and
#   key Tau associations are validated against an immunoassay (Roche Elecsys)
#   and an independent aptamer platform (SomaScan). See CHANGES.md.
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(scales)
})

# ── Paths (relative to repo root) ─────────────────────────────────────────────
DATA_BIN <- "data/sample/"
FIG_OUT  <- "output/figures/"
dir.create(FIG_OUT, recursive = TRUE, showWarnings = FALSE)

# ── Theme ─────────────────────────────────────────────────────────────────────
theme_paper <- theme_bw(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 18, hjust = 0),
    axis.title       = element_text(size = 11),
    axis.text        = element_text(size = 9),
    legend.text      = element_text(size = 8),
    legend.title     = element_text(size = 9, face = "bold"),
    legend.position  = "bottom",
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "#EBF5FB"),
    strip.text       = element_text(face = "bold")
  )

dx_colors <- c("CN" = "#2e8b57", "MCI" = "#e08e0b", "DEM" = "#c0392b")

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

# ── ADNI Emory TMT-MS loader (Fig5, Fig6, Table2) ─────────────────────────────
#
#   em <- load_adni_proteomics(
#     emory_path     = "D:/work/emory_results/",
#     adnimerge_path = "D:/work/ADNIMERGE2/ADNIMERGE2/data/")
#
# Returns the Emory TMT-MS matrix merged with diagnosis (DXSUM) and demographics
# (ADSL). n = 1,105 subjects (CN = 379, MCI = 562, DEM = 164).
load_adni_proteomics <- function(emory_path, adnimerge_path) {

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
      env <- new.env(); load(rda_files[1], envir = env)
      prot <- get(ls(env)[1], envir = env)
    }
  } else {
    stop("Emory TMT-MS file not found in: ", emory_path, "\n",
         "Download from https://adni.loni.usc.edu (search 'Emory' or 'TMT').")
  }

  prot_cols <- names(prot)[grepl("_[A-Z][0-9]", names(prot))]
  for (col in prot_cols)
    if (is.numeric(prot[[col]])) prot[[col]][prot[[col]] == -4] <- NA

  dxsum_file <- file.path(adnimerge_path, "DXSUM.rda")
  if (!file.exists(dxsum_file)) stop("DXSUM.rda not found: ", dxsum_file)
  load(dxsum_file)   # → DXSUM

  if ("DXCHANGE" %in% names(DXSUM)) {
    dx_map <- c("1"="CN","2"="MCI","3"="MCI","4"="MCI",
                "5"="CN","6"="MCI","7"="DEM","8"="DEM")
    DXSUM$DX <- dx_map[as.character(DXSUM$DXCHANGE)]
  } else if ("DIAGNOSIS" %in% names(DXSUM)) {
    dx_map <- c("CN"="CN","MCI"="MCI","Dementia"="DEM","AD"="DEM")
    DXSUM$DX <- dx_map[DXSUM$DIAGNOSIS]
  } else stop("Cannot find DXCHANGE or DIAGNOSIS in DXSUM.rda")

  DXSUM <- DXSUM[order(DXSUM$RID, DXSUM$VISCODE), ]
  DXSUM_last <- DXSUM[!duplicated(DXSUM$RID, fromLast = TRUE), ]

  adsl_file <- file.path(adnimerge_path, "ADSL.rda")
  if (!file.exists(adsl_file)) stop("ADSL.rda not found: ", adsl_file)
  load(adsl_file)    # → ADSL
  keep_cols <- intersect(c("RID","AGE","PTGENDER","PTEDUCAT"), names(ADSL))
  ADSL_sub  <- ADSL[!duplicated(ADSL$RID), keep_cols]

  if (!"RID" %in% names(prot) && "SUBJID" %in% names(prot))
    prot$RID <- as.integer(gsub("[^0-9]", "", prot$SUBJID))
  if (!"RID" %in% names(prot)) stop("Cannot find RID column in Emory file.")

  merged <- merge(prot,  DXSUM_last[, c("RID","DX")], by="RID", all.x=TRUE)
  merged <- merge(merged, ADSL_sub,                    by="RID", all.x=TRUE)
  if ("PTGENDER" %in% names(merged) && !"SEX" %in% names(merged))
    merged$SEX <- merged$PTGENDER

  merged <- merged[!is.na(merged$DX) & merged$DX %in% c("CN","MCI","DEM"), ]
  merged$DX <- factor(merged$DX, levels = c("CN","MCI","DEM"))

  cat(sprintf("ADNI Emory loaded: n=%d  CN=%d  MCI=%d  DEM=%d\n",
              nrow(merged), sum(merged$DX=="CN"),
              sum(merged$DX=="MCI"), sum(merged$DX=="DEM")))
  merged
}

# ── Roche Elecsys immunoassay loader (orthogonal Tau validation) ──────────────
#
# Returns data.frame(RID, ElecsysTau) — baseline total-tau per subject.
# File: UPENNBIOMK_ROCHE_ELECSYS_*.csv (ADNI study data).
load_elecsys_tau <- function(adni_path, tau_col = "TAU", visit = "bl") {
  f <- list.files(adni_path,
                  pattern = "(?i)(upennbiomk).*(roche|elecsys).*\\.csv$",
                  full.names = TRUE)
  if (length(f) == 0)
    f <- list.files(adni_path, pattern = "(?i)elecsys.*\\.csv$", full.names = TRUE)
  if (length(f) == 0)
    stop("Elecsys immunoassay file not found in: ", adni_path)
  message("Loading Elecsys: ", basename(f[1]))
  d <- read.csv(f[1], stringsAsFactors = FALSE)
  if ("VISCODE2" %in% names(d)) d <- d[d$VISCODE2 == visit | is.na(d$VISCODE2), ]
  d[[tau_col]] <- suppressWarnings(as.numeric(gsub("[<>]", "", d[[tau_col]])))
  d <- d[!is.na(d$RID) & !is.na(d[[tau_col]]), ]
  d <- d[!duplicated(d$RID), ]
  data.frame(RID = d$RID, ElecsysTau = d[[tau_col]])
}

# ── SomaScan aptamer loader (independent proteomic platform) ──────────────────
#
# Returns data.frame(RID, HK1_s, MAPT_s, TFRC_s). SeqIDs are dataset-specific;
# defaults match the ADNI Cruchaga-lab SOMAscan7k post-QC matrix.
load_somascan <- function(adni_path,
                          seq = c(HK1_s = "X13131.5",
                                  MAPT_s = "X5854.60",
                                  TFRC_s = "X6895.1")) {
  f <- list.files(adni_path,
                  pattern = "(?i)soma.*matrix.*\\.csv$", full.names = TRUE)
  if (length(f) == 0)
    f <- list.files(adni_path, pattern = "(?i)soma.*\\.csv$", full.names = TRUE)
  if (length(f) == 0)
    stop("SomaScan matrix not found in: ", adni_path)
  message("Loading SomaScan: ", basename(f[1]))
  d <- read.csv(f[1], stringsAsFactors = FALSE, check.names = TRUE)
  have <- seq[seq %in% names(d)]
  if (length(have) == 0)
    stop("None of the SomaScan seqIDs found; edit `seq` for your matrix version.")
  out <- data.frame(RID = d$RID)
  for (nm in names(have)) out[[nm]] <- suppressWarnings(as.numeric(d[[have[nm]]]))
  out <- out[!is.na(out$RID), ]
  out[!duplicated(out$RID), ]
}

# ── Statistics helpers ────────────────────────────────────────────────────────
# Zero-order Spearman (primary) with raw Pearson available as sensitivity.
spearman_r <- function(x, y) {
  ok <- complete.cases(x, y)
  if (sum(ok) < 10) return(list(rho = NA, p = NA, n = sum(ok)))
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman"))
  list(rho = unname(ct$estimate), p = ct$p.value, n = sum(ok))
}
pearson_r <- function(x, y) {
  ok <- complete.cases(x, y)
  if (sum(ok) < 10) return(list(r = NA, p = NA, n = sum(ok)))
  ct <- cor.test(x[ok], y[ok])
  list(r = unname(ct$estimate), p = ct$p.value, n = sum(ok))
}

# Spearman PARTIAL correlation — rank-transform x, y and every covariate, then
# residualize the ranked x and y on the ranked control set by OLS and correlate
# the residuals (Pearson of residuals). This is the primary partial method used
# throughout the CSF analysis.
spearman_partial <- function(x, y, controls) {
  ctrl <- as.data.frame(controls)
  ok <- stats::complete.cases(x, y, ctrl)
  if (sum(ok) < 10) return(list(r = NA, p = NA, n = sum(ok)))
  rx <- rank(x[ok]); ry <- rank(y[ok])
  Zr <- as.data.frame(lapply(ctrl[ok, , drop = FALSE], rank))
  ex <- residuals(lm(rx ~ ., data = Zr))
  ey <- residuals(lm(ry ~ ., data = Zr))
  ct <- cor.test(ex, ey)
  list(r = unname(ct$estimate), p = ct$p.value, n = sum(ok))
}

# Legacy Pearson partial helpers (retained for backward compatibility only;
# NOT used for the primary CSF analysis — use spearman_partial instead).
partial_cor <- function(x, y, z) {
  ok <- complete.cases(x, y, z)
  if (sum(ok) < 10) return(list(r=NA, p=NA, rx=NA, ry=NA))
  rx <- residuals(lm(x[ok] ~ z[ok])); ry <- residuals(lm(y[ok] ~ z[ok]))
  ct <- cor.test(rx, ry); list(r=ct$estimate, p=ct$p.value, rx=rx, ry=ry)
}
partial_cor_double <- function(x, y, z1, z2) {
  ok <- complete.cases(x, y, z1, z2)
  if (sum(ok) < 10) return(list(r=NA, p=NA))
  rx <- residuals(lm(x[ok] ~ z1[ok] + z2[ok])); ry <- residuals(lm(y[ok] ~ z1[ok] + z2[ok]))
  ct <- cor.test(rx, ry); list(r=ct$estimate, p=ct$p.value)
}
partial_cor_triple <- function(x, y, z1, z2, z3) {
  ok <- complete.cases(x, y, z1, z2, z3)
  if (sum(ok) < 10) return(list(r=NA, p=NA))
  rx <- residuals(lm(x[ok] ~ z1[ok] + z2[ok] + z3[ok]))
  ry <- residuals(lm(y[ok] ~ z1[ok] + z2[ok] + z3[ok]))
  ct <- cor.test(rx, ry); list(r=ct$estimate, p=ct$p.value)
}

norm_base <- function(vals, bins) vals / vals[which.min(bins)]
get_slope <- function(df, col, bins) {
  d    <- df[df$bin %in% bins, ]
  base <- d[[col]][d$bin == min(d$bin)]
  if (is.na(base) || base == 0) return(list(beta=NA, se=NA, p=NA))
  d$norm <- d[[col]] / base
  fit    <- lm(norm ~ bin, data=d); s <- summary(fit)$coefficients
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
