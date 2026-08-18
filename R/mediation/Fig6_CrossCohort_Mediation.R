# =============================================================================
# Fig6_CrossCohort_Mediation.R
#
# Figure 6 (A-D): Cross-cohort replication, mediation, and robustness.
# NOTE (revision): this figure was Figure 6 in the submitted version. In-text
# citation order required Figure 6 = cross-cohort/mediation and Figure 6 = iron
# pathway, so the two were renumbered. Content is unchanged.
#   A  -  Astrocytic MCT4 donor-level slope on each progression axis, both cohorts.
#       The ROSMAP sign reverses between pseudotime and Braak IN THE SAME DONORS.
#   B  -  Mediation: % of the CPS effect carried by astrocytic MCT4, each forward bar
#       shown against its OWN matched null (222 detection- and effect-size-matched
#       astrocytic genes). Outcomes not exceeding the null are omitted.
#   C  -  Detection-matched null: MCT4 against the 2,117 genes of the same
#       detection rate (6.7%).
#   D  -  The metabolic chain on the shared Braak axis: MTG vs DLPFC slopes.
#
# Inputs (produced by 49_mediation_P2_definitions.R, 47_braak_common_axis.R,
#         46_composition_vs_percell.R, and 38_detection_matched_null.R):
#   tables/Table3_mediation.csv        (produced by 50_mediation_matched_null.R)
#   tables/Table4_ROSMAP.csv
#   data/47_chain_braak.csv
#   data/38_all_genes.csv
#
# Output: output/figures/Fig6_CrossCohort_Mediation.png   (14 x 10 in, 300 dpi)
# =============================================================================
# Paths are RELATIVE to the repository root. Set the working directory there,
# or override with the environment variables SEAAD_H5AD / ROSMAP_ASTRO /
# ROSMAP_CLIN / P2_OUT_DIR. Raw data are not redistributable; see README.md.

source("R/figures/utils.R")
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(dplyr) })

T3   <- read.csv("tables/Table3_mediation.csv")   # v2: incl. matched-null columns
T4   <- read.csv("tables/Table4_ROSMAP.csv",     check.names = FALSE)
CHN  <- read.csv("data/47_chain_braak.csv")
ALLG <- read.csv("data/38_all_genes.csv")

RED   <- "#c0392b"; GREEN <- "#2e8b57"; ORANGE <- "#e08e0b"; GREY <- "#7f8c8d"

## -- Panel A: the same donors, two rulers --------------------------------------
A <- T4 %>%
  mutate(label = paste0(Cohort, "\n", `Progression axis`),
         beta  = `β (astrocytic MCT4)`,
         p     = as.numeric(p),
         se    = abs(beta) / abs(qnorm(p / 2)),          # SE back out of beta and p
         lo    = beta - 1.96 * se, hi = beta + 1.96 * se,
         sig   = ifelse(p < .001, "***", ifelse(p < .01, "**",
                 ifelse(p < .05, "*", "n.s."))),
         dir   = ifelse(beta < 0, "decline", "increase"))
A$label <- factor(A$label, levels = rev(A$label))

pA <- ggplot(A, aes(beta, label, color = dir)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = .4) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = .16, color = GREY, linewidth = .6) +
  geom_point(size = 3.2) +
  geom_text(aes(label = sprintf("%s  p=%.2g  n=%d", sig, p, Donors)),
            hjust = ifelse(A$beta < 0, -0.08, 1.08), vjust = -1.1,
            size = 2.9, color = "#2c3e50", show.legend = FALSE) +
  scale_color_manual(values = c(decline = RED, increase = GREEN), guide = "none") +
  labs(x = "Astrocytic MCT4: donor-level slope (\u03b2)", y = NULL) +
  coord_cartesian(xlim = c(-0.105, 0.055)) +
  theme_paper +
  theme(panel.grid.major.y = element_blank(), legend.position = "none")

## -- Panel B: mediation vs the matched null -----------------------------------
## The mediated fraction is NOT interpretable on its own: adjusting for any gene
## strongly associated with CPS removes part of a weaker CPS->outcome effect.
## Each forward bar is therefore drawn against its own matched null (222 astrocytic
## genes matched to SLC16A3 on detection rate AND CPS effect size; see
## R/50_mediation_matched_null.R). Grey backdrop = null 95th percentile,
## black tick = null median. Outcomes not exceeding the null are not shown.
NULLG <- "#D0D3D4"

Bf <- T3 %>% filter(Direction == "MCT4 as mediator") %>%
  mutate(out = gsub(" \\(6 subunits, P2\\)", "", Outcome),
         out = gsub("^Astrocytic ", "Astro. ", out))
Br <- T3 %>% filter(Direction == "Reverse direction") %>%
  mutate(out = gsub(" \\(6, P2\\)", "", Mediator),
         out = gsub("^Astrocytic ", "Astro. ", out))

## forward rows are ordered by empirical p in the table; keep that order top-down
Bf$out <- factor(Bf$out, levels = rev(Bf$out))
Br$out <- factor(Br$out, levels = rev(Br$out))

pBf <- ggplot(Bf, aes(y = out)) +
  geom_vline(xintercept = 100, linetype = "dotted", color = GREY, linewidth = .5) +
  geom_col(aes(x = null_p95), width = .70, fill = NULLG,
           color = "#95A5A6", linewidth = .3) +
  geom_col(aes(x = pct_mediated), width = .40, fill = RED,
           color = "black", linewidth = .3) +
  geom_segment(aes(x = null_median, xend = null_median,
                   y = as.numeric(out) - .35, yend = as.numeric(out) + .35),
               color = "#2c3e50", linewidth = .7) +
  geom_text(aes(x = pct_mediated, label = sprintf("%.1f%%", pct_mediated)),
            hjust = -0.12, vjust = -0.25, size = 3.1, fontface = "bold") +
  geom_text(aes(x = pct_mediated,
                label = sprintf("null %.1f%%  p=%.3f", null_median, p_empirical)),
            hjust = -0.09, vjust = 1.45, size = 2.5, color = "#5D6D7E",
            fontface = "italic") +
  labs(x = NULL, y = NULL, subtitle = "MCT4 as MEDIATOR") +
  coord_cartesian(xlim = c(0, 175)) +
  theme_paper +
  theme(panel.grid.major.y = element_blank(),
        plot.subtitle = element_text(color = RED, face = "bold", hjust = 1))

pBr <- ggplot(Br, aes(pct_mediated, out)) +
  geom_vline(xintercept = 100, linetype = "dotted", color = GREY, linewidth = .5) +
  geom_col(width = .62, fill = GREEN, color = "black", linewidth = .3) +
  geom_text(aes(label = sprintf("%.1f%%", pct_mediated)), hjust = -0.15, size = 4.1) +
  labs(x = "% of the CPS effect that is mediated", y = NULL,
       subtitle = "REVERSE direction", caption = "n = 84 donors") +
  coord_cartesian(xlim = c(0, 175)) +
  theme_paper +
  theme(panel.grid.major.y = element_blank(),
        plot.subtitle = element_text(color = GREEN, face = "bold", hjust = 1))

pBf <- pBf + panel_tag("B")
pB <- (pBf / pBr) + patchwork::plot_layout(heights = c(3, 5))

## -- Panel C: detection-matched null ------------------------------------------
m4   <- ALLG[ALLG$gene == "SLC16A3", ]
band <- ALLG %>% filter(detect >= m4$detect * .75, detect <= m4$detect * 1.25,
                        is.finite(pct))
pctl <- 100 * mean(band$pct < m4$pct, na.rm = TRUE)

pC <- ggplot(band, aes(pct)) +
  geom_histogram(bins = 60, fill = "#BFD7EA", color = "white", linewidth = .2) +
  geom_vline(xintercept = median(band$pct, na.rm = TRUE),
             linetype = "dashed", color = GREY, linewidth = .5) +
  geom_vline(xintercept = m4$pct, color = RED, linewidth = 1.1) +
  annotate("text", x = m4$pct + 6, y = Inf, vjust = 2.2, hjust = 0,
           label = sprintf("MCT4  %.1f%%", m4$pct),
           color = RED, fontface = "bold", size = 3.4) +
  annotate("label", x = Inf, y = Inf, hjust = 1.02, vjust = 1.2, size = 4.1,
           label = sprintf("n = %s genes at the same detection rate\ndetection %.3f-%.3f (MCT4 = %.3f)\nMCT4 percentile: %.2f",
                           format(nrow(band), big.mark = ","),
                           m4$detect * .75, m4$detect * 1.25, m4$detect, pctl),
           fill = "#F4F6F7", label.size = .25) +
  labs(x = "Early\u2192late change (%)", y = "Number of genes") +
  coord_cartesian(xlim = c(-100, 120)) +
  theme_paper

## -- Panel D: the chain on the shared Braak axis ------------------------------
CHN <- CHN %>% mutate(sig_mtg = MTG_p < .05,
                      lab = ifelse(gene %in% c("SLC16A3","HK2","LDHA","PTGDS",
                                               "SLC1A2","GFAP","SERPINA3","GAPDH"),
                                   gene, NA))
pD <- ggplot(CHN, aes(MTG_beta, DLPFC_beta)) +
  geom_hline(yintercept = 0, linewidth = .4) +
  geom_vline(xintercept = 0, linewidth = .4) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = GREY) +
  geom_point(aes(color = sig_mtg), size = 2.4, alpha = .9) +
  geom_point(data = CHN[CHN$gene == "SLC16A3", ], shape = 23, size = 4.6,
             fill = RED, color = "black", stroke = .6) +
  ggrepel::geom_text_repel(aes(label = lab), size = 4.1, color = "#2c3e50",
                           na.rm = TRUE, min.segment.length = .2, seed = 42) +
  scale_color_manual(values = c(`TRUE` = ORANGE, `FALSE` = "#D5DBDB"),
                     labels = c(`TRUE` = "p < 0.05 in MTG", `FALSE` = "n.s. in MTG"),
                     name = NULL) +
  labs(x = "MTG (SEA-AD): slope on Braak",
       y = "DLPFC (ROSMAP): slope on Braak",
       caption = sprintf("%d / %d genes agree in sign", sum(CHN$same_sign), nrow(CHN))) +
  coord_cartesian(xlim = c(-.11, .11), ylim = c(-.11, .11)) +
  theme_paper + theme(legend.position = c(.22, .93),
                      legend.background = element_rect(fill = "white", color = "black",
                                                       linewidth = .3))

## One tag mechanism for all four panels. Panel B is a wrap_elements grob, whose
## tag patchwork insets by a fixed amount, so its tag position is nudged to sit
## on the same baseline as A, C and D rather than 46 px below them.
pA <- pA + panel_tag("A")
pC <- pC + panel_tag("C")
pD <- pD + panel_tag("D")
fig7 <- (pA | pB) / (pC | pD)
save_fig(fig7, "Fig6_CrossCohort_Mediation.png", width = 14, height = 10)
