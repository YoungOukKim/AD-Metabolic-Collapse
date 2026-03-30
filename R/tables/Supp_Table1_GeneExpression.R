# =============================================================================
# Supp_Table1_GeneExpression.R
#
# Supplementary Table 1: Astrocyte gene expression values across
# pseudo-progression bins for all analyzed genes
#
# Reproduces the table in the paper showing:
#   - Gene name and functional category
#   - Mean expression: Bins 0.2-0.4 (early) vs Bins 0.6-0.8 (late)
#   - % change early → late
#
# Input:  data/sample/astro_bin_means.csv
# Output: output/tables/Supp_Table1_GeneExpression.csv
# =============================================================================
source("R/figures/utils.R")
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)

astro <- read.csv(file.path(DATA_BIN, "astro_bin_means.csv"))

# ── Functional category map ───────────────────────────────────────────────────
cat_map <- c(
  # ANLS / glycolytic
  SLC16A3 = "ANLS",  SLC16A1 = "ANLS",  SLC2A1  = "ANLS",
  LDHA    = "ANLS",  HK2     = "ANLS",  HK1     = "ANLS",
  PKM     = "ANLS",  PDK1    = "ANLS",  PFKFB3  = "ANLS",
  LDHB    = "ANLS",  GPI     = "ANLS",  ENO1    = "ANLS",
  ENO2    = "ANLS",  GAPDH   = "ANLS",
  # V-ATPase
  ATP6V1A = "V-ATPase", ATP6V1B2 = "V-ATPase", ATP6V1C1 = "V-ATPase",
  ATP6V1E1= "V-ATPase", ATP6V1H  = "V-ATPase", ATP6V0A1 = "V-ATPase",
  ATP6V0B = "V-ATPase", ATP6V0C  = "V-ATPase", ATP6V0D1 = "V-ATPase",
  ATP6V0E1= "V-ATPase",
  # Lysosomal
  LAMP1   = "Lysosomal", LAMP2 = "Lysosomal",
  CTSB    = "Lysosomal", CTSD  = "Lysosomal", LIPA = "Lysosomal",
  # Iron
  TFRC    = "Iron", FTH1 = "Iron", FTL  = "Iron",
  CP      = "Iron", SLC40A1 = "Iron",
  # Mito / Ragulator
  VDAC1   = "Mitochondrial", NDUFS1  = "Mitochondrial",
  LAMTOR1 = "Ragulator", LAMTOR2 = "Ragulator", LAMTOR3 = "Ragulator",
  LAMTOR4 = "Ragulator", LAMTOR5 = "Ragulator",
  MTOR    = "mTOR-Ragulator",
  # pH
  SLC9A6  = "pH regulatory",
  # Signaling / checkpoint
  PTGDS   = "Checkpoint", LCN2  = "Inflammation",
  HMOX1   = "Oxidative stress", SOX9 = "Signaling",
  STAT3   = "Signaling",
  # Metabolic overload
  SLC1A2  = "Astrocyte transport", ATP1A2 = "Astrocyte transport",
  # Tau / trophic
  MAPT    = "Tau"
)

# ── Compute early/late means and % change ────────────────────────────────────
early_bins <- c(0.2, 0.3, 0.4)
late_bins  <- c(0.6, 0.7, 0.8)

# All numeric gene columns (exclude composite columns)
composite_cols <- c("ANLS","VATpase","MCT4","ED_energy","ED_ratio","VATPase","bin")
gene_cols <- setdiff(names(astro)[sapply(astro, is.numeric)], composite_cols)

results <- do.call(rbind, lapply(gene_cols, function(g) {
  early_val <- mean(astro[[g]][astro$bin %in% early_bins], na.rm = TRUE)
  late_val  <- mean(astro[[g]][astro$bin %in% late_bins],  na.rm = TRUE)
  pct_change <- (late_val - early_val) / early_val * 100

  # Gene label (add alias for known genes)
  alias_map <- c(SLC16A3="MCT4", SLC16A1="MCT1", SLC16A7="MCT2",
                 SLC2A1="GLUT1", SLC2A3="GLUT3", SLC9A6="NHE6",
                 SLC1A2="EAAT2", SLC40A1="FPN1")
  label <- if (g %in% names(alias_map)) paste0(g, " (", alias_map[g], ")") else g

  data.frame(
    Gene          = label,
    Category      = ifelse(g %in% names(cat_map), cat_map[g], "Other"),
    Bins_0.2_0.4  = round(early_val,  4),
    Bins_0.6_0.8  = round(late_val,   4),
    Pct_Change    = round(pct_change, 1),
    stringsAsFactors = FALSE
  )
}))

# Sort by category then % change
results <- results[order(results$Category, results$Pct_Change), ]

# ── Print summary ─────────────────────────────────────────────────────────────
cat("\n=== Supplementary Table 1: Gene Expression Summary ===\n")
cat(sprintf("Total genes: %d\n\n", nrow(results)))

# Category-level summary
cat_summary <- aggregate(Pct_Change ~ Category, data = results, FUN = mean)
cat_summary <- cat_summary[order(cat_summary$Pct_Change), ]
print(cat_summary, digits = 2)

cat("\n--- Top 10 decreasing genes ---\n")
print(head(results[order(results$Pct_Change), c("Gene","Category","Pct_Change")], 10))
cat("\n--- Top 10 increasing genes ---\n")
print(head(results[order(-results$Pct_Change), c("Gene","Category","Pct_Change")], 10))

# ── Save ──────────────────────────────────────────────────────────────────────
write.csv(results, "output/tables/Supp_Table1_GeneExpression.csv", row.names = FALSE)
cat("\nSaved: output/tables/Supp_Table1_GeneExpression.csv\n")

# Note on full 179-gene version
cat("\nNote: Full 179-gene table requires running 01_extract_seaad.R first\n")
cat("      to generate astro_comprehensive.csv with all analyzed genes.\n")
