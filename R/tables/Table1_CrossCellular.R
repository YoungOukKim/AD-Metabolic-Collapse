# =============================================================================
# Table1_CrossCellular.R
#
# Table 1: Cross-cellular correlations, partial correlations, and
#          donor-level clinical associations
#
# Reproduces all statistics in Table 1 of the paper.
#
# Input:  data/sample/astro_bin_means.csv
#         data/sample/neuron_bin_means.csv
#         data/sample/donor_level_summary.csv
# Output: printed to console (copy into manuscript)
#         output/tables/Table1_CrossCellular.csv
# =============================================================================
source("R/figures/utils.R")
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)

astro  <- read.csv(file.path(DATA_BIN, "astro_bin_means.csv"))
neuron <- read.csv(file.path(DATA_BIN, "neuron_bin_means.csv"))
donor  <- read.csv(file.path(DATA_BIN, "donor_level_summary.csv"))

astro  <- add_composites_astro(astro)
neuron <- add_composites_neuron(neuron)
astro  <- astro[astro$bin   >= 0.2,]
neuron <- neuron[neuron$bin >= 0.2,]
donor  <- donor[!is.na(donor$MCT4) & !is.na(donor$VATpase) & !is.na(donor$mean_cps),]

df <- data.frame(bin=astro$bin, MCT4=astro$MCT4, ANLS=astro$ANLS,
                 VATPase=neuron$VATPase, LDHB=neuron$LDHB, LAMP1=neuron$LAMP1)

results <- data.frame()

# ── Bin-level correlations ────────────────────────────────────────────────────
pairs_bin <- list(
  list(x="MCT4", y="VATPase", label="MCT4 -> Neuron V-ATPase"),
  list(x="MCT4", y="LDHB",    label="MCT4 -> Neuron LDHB"),
  list(x="MCT4", y="LAMP1",   label="MCT4 -> Neuron LAMP1"),
  list(x="ANLS", y="VATPase", label="ANLS -> Neuron V-ATPase")
)

for (pr in pairs_bin) {
  r0  <- cor.test(df[[pr$x]], df[[pr$y]])
  pc  <- partial_cor(df[[pr$x]], df[[pr$y]], df$bin)
  results <- rbind(results, data.frame(
    Level       = "Bin (n=8)",
    Measure     = pr$label,
    r_zero      = round(r0$estimate, 3),
    p_zero      = r0$p.value,
    r_partial   = round(pc$r, 3),
    p_partial   = pc$p,
    survives_CPS = !is.na(pc$r) && abs(pc$r) >= 0.7 && pc$p < 0.05
  ))
}

# ── Donor-level ───────────────────────────────────────────────────────────────
r0_d  <- cor.test(donor$MCT4, donor$VATpase)
pc_d  <- partial_cor(donor$MCT4, donor$VATpase, donor$mean_cps)
results <- rbind(results, data.frame(
  Level     = "Donor (n=84)",
  Measure   = "MCT4 -> Neuron V-ATPase",
  r_zero    = round(r0_d$estimate, 3),
  p_zero    = r0_d$p.value,
  r_partial = round(pc_d$r, 3),
  p_partial = pc_d$p,
  survives_CPS = !is.na(pc_d$r) && abs(pc_d$r) >= 0.7 && pc_d$p < 0.05
))

# ── Slope dissociation ────────────────────────────────────────────────────────
bins_all <- seq(0.2, 0.9, 0.1)
s_mct4 <- get_slope(astro, "MCT4",    bins_all)
s_anls <- get_slope(astro, "ANLS",    bins_all)
s_vatp <- get_slope(astro, "VATpase", bins_all)

delta_z_mct4 <- (s_mct4$beta - s_vatp$beta) / sqrt(s_mct4$se^2 + s_vatp$se^2)
delta_p_mct4 <- 2 * pnorm(-abs(delta_z_mct4))

cat("\n=== Table 1: Cross-cellular correlations ===\n")
print(results)

cat(sprintf("\n--- Slope dissociation ---\n"))
cat(sprintf("MCT4: beta=%.3f, SE=%.3f, p=%.4f\n", s_mct4$beta, s_mct4$se, s_mct4$p))
cat(sprintf("ANLS: beta=%.3f, SE=%.3f, p=%.4f\n", s_anls$beta, s_anls$se, s_anls$p))
cat(sprintf("V-ATPase: beta=%.3f, SE=%.3f, p=%.4f\n", s_vatp$beta, s_vatp$se, s_vatp$p))
cat(sprintf("Delta slope MCT4 vs V-ATPase: z=%.3f, p=%.4f\n", delta_z_mct4, delta_p_mct4))

# ── Donor-level clinical correlations ─────────────────────────────────────────
if ("braak_num" %in% names(donor)) {
  rho_anls_braak  <- cor.test(donor$ANLS,    donor$braak_num, method="spearman", exact=FALSE)
  rho_mct4_braak  <- cor.test(donor$MCT4,    donor$braak_num, method="spearman", exact=FALSE)
  rho_vatp_braak  <- cor.test(donor$VATpase, donor$braak_num, method="spearman", exact=FALSE)
  cat(sprintf("\n--- Braak stage correlations ---\n"))
  cat(sprintf("ANLS-Braak: rho=%.3f, p=%.4f\n", rho_anls_braak$estimate, rho_anls_braak$p.value))
  cat(sprintf("MCT4-Braak: rho=%.3f, p=%.4f\n", rho_mct4_braak$estimate, rho_mct4_braak$p.value))
  cat(sprintf("V-ATPase-Braak: rho=%.3f, p=%.4f\n", rho_vatp_braak$estimate, rho_vatp_braak$p.value))
}

write.csv(results, "output/tables/Table1_CrossCellular.csv", row.names=FALSE)
cat("\nSaved: output/tables/Table1_CrossCellular.csv\n")
