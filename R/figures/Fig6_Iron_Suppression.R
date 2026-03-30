# =============================================================================
# Fig6_Iron_Suppression.R
#
# Figure 6 (a): Iron pathway — TFRC-MAPT suppression effect
#
#   a — Bar chart showing TFRC-MAPT correlation under progressive control:
#       None → | GFAP → | V1A
#       Suppression effect: r strengthens after V1A control (+0.378 → +0.602)
#
# Input (requires ADNI data access — not included in repository):
#   Emory TMT-MS CSF proteomics + DXSUM.rda + ADSL.rda
#   See load_adni_proteomics() in R/figures/utils.R for details.
#
# Output: output/figures/Fig6_Iron_Suppression.png
# =============================================================================
source("R/figures/utils.R")

# ── Set your local paths here ─────────────────────────────────────────────────
EMORY_PATH     <- "D:/work/emory_results/"
ADNIMERGE_PATH <- "D:/work/ADNIMERGE2/ADNIMERGE2/data/"
# ─────────────────────────────────────────────────────────────────────────────

em <- load_adni_proteomics(EMORY_PATH, ADNIMERGE_PATH)

# ── Protein column mapping ────────────────────────────────────────────────────
P <- list(
  V1A  = "ATP6V1A_P38606",
  MAPT = "MAPT_P10636",
  GFAP = "GFAP_P14136",
  TFRC = "TFRC_P02786"
)
P_found   <- P[sapply(P, function(x) x %in% names(em))]
P_missing <- names(P)[!names(P) %in% names(P_found)]
if (length(P_missing) > 0)
  message("Warning — proteins not found: ", paste(P_missing, collapse=", "))
P <- P_found

if (!all(c("TFRC","MAPT","GFAP","V1A") %in% names(P)))
  stop("Required proteins missing. Check column names in Emory dataset.")

# ── Compute correlations ──────────────────────────────────────────────────────
r_none  <- cor(em[[P$TFRC]], em[[P$MAPT]], use="complete.obs")
pc_gfap <- partial_cor(em[[P$TFRC]], em[[P$MAPT]], em[[P$GFAP]])
pc_v1a  <- partial_cor(em[[P$TFRC]], em[[P$MAPT]], em[[P$V1A]])

cat(sprintf("TFRC-MAPT r (none):    %.3f\n", r_none))
cat(sprintf("TFRC-MAPT r (| GFAP): %.3f  p=%.2e\n", pc_gfap$r, pc_gfap$p))
cat(sprintf("TFRC-MAPT r (| V1A):  %.3f  p=%.2e  [suppression effect]\n",
            pc_v1a$r, pc_v1a$p))

supp_dat <- data.frame(
  Control = factor(c("None","| GFAP","| V1A"),
                   levels = c("None","| GFAP","| V1A")),
  r = c(r_none, pc_gfap$r, pc_v1a$r)
)

# ── Panel a ───────────────────────────────────────────────────────────────────
fig6a <- ggplot(supp_dat, aes(x=Control, y=r, fill=Control)) +
  geom_col(color="black", width=0.6) +
  geom_text(aes(label=sprintf("%.3f", r)),
            vjust=-0.3, size=4, fontface="bold") +
  scale_fill_manual(values=c("None"="#e67e22","| GFAP"="#f39c12","| V1A"="#c0392b")) +
  # Arrow showing suppression direction
  geom_segment(aes(x=1, xend=3, y=max(r)+0.05, yend=max(r)+0.05),
               arrow=arrow(length=unit(0.15,"cm")), linewidth=0.5) +
  annotate("text", x=2, y=max(supp_dat$r)+0.08,
           label="Suppression effect", size=3.5,
           fontface="italic", color="#c0392b") +
  ylim(0, max(supp_dat$r) + 0.15) +
  labs(x="TFRC \u2194 MAPT: partial correlation control",
       y="Pearson r", title="a") +
  theme_paper + theme(legend.position="none")

save_fig(fig6a, "Fig6_Iron_Suppression.png", width=6, height=5)
