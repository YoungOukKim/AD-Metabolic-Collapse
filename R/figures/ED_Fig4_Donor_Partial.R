# =============================================================================
# ED_Fig4_Donor_Partial.R
#
# Extended Data Figure 4 (a-b): Donor-level partial correlation (CPS-adjusted)
#
#   a — Zero-order: MCT4 vs Neuronal V-ATPase (r = +0.807, n=84)
#   b — CPS-adjusted partial correlation (partial r = +0.817)
#
# Input:  data/sample/donor_level_summary.csv
# Output: output/figures/ED_Fig4_Donor_Partial.png
# =============================================================================
source("R/figures/utils.R")

donor_file <- file.path(DATA_BIN, "donor_level_summary.csv")
if (!file.exists(donor_file)) stop("donor_level_summary.csv not found")

donor_clean <- read.csv(donor_file)
donor_clean <- donor_clean[!is.na(donor_clean$MCT4) &
                           !is.na(donor_clean$VATpase) &
                           !is.na(donor_clean$mean_cps),]

r0_d <- cor.test(donor_clean$MCT4, donor_clean$VATpase)
pc_d <- partial_cor(donor_clean$MCT4, donor_clean$VATpase, donor_clean$mean_cps)

ed4a <- ggplot(donor_clean, aes(x = MCT4, y = VATpase)) +
  geom_smooth(method="lm", se=TRUE, color="#2166AC") +
  geom_point(size=2, alpha=0.6, color="#2166AC") +
  annotate("text", x=min(donor_clean$MCT4)*1.02, y=max(donor_clean$VATpase)*0.95,
           label=sprintf("r = %+.3f\np = %.1e\nn = %d",
                         r0_d$estimate, r0_d$p.value, nrow(donor_clean)),
           size=3.5, fontface="bold", hjust=0) +
  labs(x="Astrocytic MCT4", y="Neuronal V-ATPase", title="a  Zero-order") +
  theme_paper

resid_d <- data.frame(MCT4_res=pc_d$rx, VATpase_res=pc_d$ry)
ed4b <- ggplot(resid_d, aes(x=MCT4_res, y=VATpase_res)) +
  geom_smooth(method="lm", se=TRUE, color="#B2182B") +
  geom_point(size=2, alpha=0.6, color="#D6604D") +
  annotate("text", x=min(resid_d$MCT4_res)*0.95, y=max(resid_d$VATpase_res)*0.95,
           label=sprintf("partial r = %+.3f\np = %.1e", pc_d$r, pc_d$p),
           size=3.5, fontface="bold", hjust=0) +
  labs(x="MCT4 residual (CPS-adj)",
       y="V-ATPase residual (CPS-adj)", title="b  Partial (CPS-adjusted)") +
  theme_paper

ed_fig4 <- ed4a + ed4b
save_fig(ed_fig4, "ED_Fig4_Donor_Partial.png", width=12, height=5.5)
