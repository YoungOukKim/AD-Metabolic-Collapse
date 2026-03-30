# =============================================================================
# Fig6_Iron_Suppression.R
#
# Figure 6 (a): Iron pathway — TFRC-MAPT suppression effect
#
#   a — Bar chart: TFRC-MAPT r under progressive V-ATPase control
#       showing suppression effect (r strengthens after V1A control)
#
# Input:  ADNI Emory TMT-MS CSF proteomics
#         (same data object as Fig5_CSF_Validation.R)
# Output: output/figures/Fig6_Iron_Suppression.png
# =============================================================================
source("R/figures/utils.R")

emory_file <- "data/processed/merged_proteomics_dx.rds"
if (!file.exists(emory_file)) stop("ADNI proteomics file not found")

merged <- readRDS(emory_file)
prot_cols <- names(merged)[grepl("_[A-Z][0-9]", names(merged))]
for (col in prot_cols) if (is.numeric(merged[[col]]))
  merged[[col]][merged[[col]] == -4] <- NA

P <- list(
  V1A  = "ATP6V1A_P38606",
  MAPT = "MAPT_P10636",
  GFAP = "GFAP_P14136",
  TFRC = "TFRC_P02786"
)
P <- P[sapply(P, function(x) x %in% names(merged))]
em <- merged[merged$DX %in% c("CN","MCI","DEM"),]

# ── Compute r values ──────────────────────────────────────────────────────────
r_none <- cor(em[[P$TFRC]], em[[P$MAPT]], use = "complete.obs")
pc_gfap <- partial_cor(em[[P$TFRC]], em[[P$MAPT]], em[[P$GFAP]])
pc_v1a  <- partial_cor(em[[P$TFRC]], em[[P$MAPT]], em[[P$V1A]])

supp_dat <- data.frame(
  Control = factor(c("None","| GFAP","| V1A"),
                   levels = c("None","| GFAP","| V1A")),
  r = c(r_none, pc_gfap$r, pc_v1a$r)
)

# ── Panel a ───────────────────────────────────────────────────────────────────
fig6a <- ggplot(supp_dat, aes(x = Control, y = r, fill = Control)) +
  geom_col(color = "black", width = 0.6) +
  geom_text(aes(label = sprintf("%.3f", r)),
            vjust = -0.3, size = 4, fontface = "bold") +
  scale_fill_manual(values = c("None"="#e67e22","| GFAP"="#f39c12","| V1A"="#c0392b")) +
  geom_segment(aes(x = 1, xend = 3, y = max(r) + 0.05, yend = max(r) + 0.05),
               arrow = arrow(length = unit(0.15,"cm")), linewidth = 0.5) +
  annotate("text", x = 2, y = max(supp_dat$r) + 0.08,
           label = "Suppression effect", size = 3.5,
           fontface = "italic", color = "#c0392b") +
  ylim(0, max(supp_dat$r) + 0.15) +
  labs(x = "TFRC \u2194 MAPT partial correlation control",
       y = "Pearson r", title = "a") +
  theme_paper + theme(legend.position = "none")

save_fig(fig6a, "Fig6_Iron_Suppression.png", width = 6, height = 5)
