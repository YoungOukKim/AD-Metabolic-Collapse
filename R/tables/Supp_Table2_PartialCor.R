# =============================================================================
# Supp_Table2_PartialCor.R
#
# Supplementary Table 2: Complete partial correlation analysis
# (transcriptomic) at bin and donor levels
#
# Reproduces the full table showing:
#   - Zero-order r between astrocyte and neuron markers
#   - CPS-adjusted partial r
#   - % change in |r| after adjustment
#   - Interpretation (Survives / Attenuated / Lost)
#
# Input:  data/sample/astro_bin_means.csv
#         data/sample/neuron_bin_means.csv
#         data/sample/donor_level_summary.csv
# Output: output/tables/Supp_Table2_PartialCor.csv
# =============================================================================
source("R/figures/utils.R")
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)

astro  <- read.csv(file.path(DATA_BIN, "astro_bin_means.csv"))
neuron <- read.csv(file.path(DATA_BIN, "neuron_bin_means.csv"))
donor  <- read.csv(file.path(DATA_BIN, "donor_level_summary.csv"))

astro  <- add_composites_astro(astro)
neuron <- add_composites_neuron(neuron)
astro  <- astro[astro$bin   >= 0.2, ]
neuron <- neuron[neuron$bin >= 0.2, ]
donor  <- donor[!is.na(donor$MCT4) & !is.na(donor$VATpase) &
                !is.na(donor$mean_cps), ]

# Merge for bin-level analysis
df <- data.frame(
  bin    = astro$bin,
  MCT4   = astro$MCT4,
  ANLS   = astro$ANLS,
  VATPase= neuron$VATPase,
  LDHB   = if ("LDHB" %in% names(neuron)) neuron$LDHB else NA,
  LAMP1  = if ("LAMP1" %in% names(neuron)) neuron$LAMP1 else NA,
  NDUFS1 = if ("NDUFS1" %in% names(neuron)) neuron$NDUFS1 else NA
)

# ── Helper: compute one row ───────────────────────────────────────────────────
cor_row <- function(level, astro_m, neuron_m, x, y, z, n_label) {
  ok <- complete.cases(x, y, z)
  if (sum(ok) < 5) return(NULL)

  r0 <- cor.test(x, y)
  pc <- partial_cor(x, y, z)

  pct_change <- if (!is.na(pc$r) && !is.na(r0$estimate))
    (abs(pc$r) - abs(r0$estimate)) / abs(r0$estimate) * 100 else NA

  interpret <- if (is.na(pc$r)) "NA"
  else if (abs(pc$r) >= 0.7 && pc$p < 0.05) "**Survives CPS**"
  else if (abs(pc$r) >= 0.3 && pc$p < 0.05) "Attenuated"
  else "Lost"

  data.frame(
    Level        = level,
    N            = n_label,
    Astro_marker = astro_m,
    Neuron_marker= neuron_m,
    r_zero       = round(r0$estimate, 3),
    p_zero       = formatC(r0$p.value, format="e", digits=1),
    r_partial    = round(pc$r, 3),
    p_partial    = if(!is.na(pc$p)) formatC(pc$p, format="e", digits=1) else "NA",
    Pct_change   = round(pct_change, 1),
    Interpretation = interpret,
    stringsAsFactors = FALSE
  )
}

# ── Bin-level pairs ───────────────────────────────────────────────────────────
bin_pairs <- list(
  list("MCT4",  "VATPase", df$MCT4,  df$VATPase, df$bin),
  list("MCT4",  "LDHB",    df$MCT4,  df$LDHB,    df$bin),
  list("MCT4",  "LAMP1",   df$MCT4,  df$LAMP1,   df$bin),
  list("MCT4",  "NDUFS1",  df$MCT4,  df$NDUFS1,  df$bin),
  list("ANLS",  "VATPase", df$ANLS,  df$VATPase, df$bin),
  list("ANLS",  "LDHB",    df$ANLS,  df$LDHB,    df$bin),
  list("ANLS",  "NDUFS1",  df$ANLS,  df$NDUFS1,  df$bin)
)

results_bin <- do.call(rbind, lapply(bin_pairs, function(pr) {
  cor_row("Bin", pr[[1]], pr[[2]], pr[[3]], pr[[4]], pr[[5]],
          sprintf("n=%d bins", nrow(df)))
}))

# ── Donor-level pairs ─────────────────────────────────────────────────────────
donor_pairs <- list(
  list("MCT4", "VATpase",   donor$MCT4, donor$VATpase,  donor$mean_cps),
  list("ANLS", "VATpase",   donor$ANLS, donor$VATpase,  donor$mean_cps)
)

results_donor <- do.call(rbind, lapply(donor_pairs, function(pr) {
  cor_row("Donor", pr[[1]], pr[[2]], pr[[3]], pr[[4]], pr[[5]],
          sprintf("n=%d donors", nrow(donor)))
}))

results <- rbind(results_bin, results_donor)

# ── Print ─────────────────────────────────────────────────────────────────────
cat("\n=== Supplementary Table 2: Complete Partial Correlation Analysis ===\n\n")
print(results, row.names = FALSE)

# Summary counts
cat(sprintf("\nSurvives CPS control: %d / %d pairs\n",
            sum(grepl("Survives", results$Interpretation)),
            nrow(results)))

# ── Save ──────────────────────────────────────────────────────────────────────
write.csv(results, "output/tables/Supp_Table2_PartialCor.csv", row.names = FALSE)
cat("Saved: output/tables/Supp_Table2_PartialCor.csv\n")
