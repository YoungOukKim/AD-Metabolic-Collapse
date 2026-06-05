# =============================================================================
# Fig4_Clinical.R  [CORRECTED — donor-level 3-panel; matches published Figure 4]
#
# Figure 4 (A-C): Donor-level clinical validation (true donor-level, n = 84).
#   A — Spearman rho of ANLS / MCT4 / astrocytic V-ATPase vs Braak, CERAD, ABC
#   B — MCT4 vs CPS (weighted-style donor regression; R^2, p)
#   C — MCT4 by dementia status (Mann-Whitney)
#
# NOTE: This replaces the previous 6-panel per-staging-boxplot version so that
# the script reproduces the published Figure 4. The clinical V-ATPase comparator
# is the ASTROCYTIC composite (consistent with the dissociation thesis).
#
# Input:  donor_level_summary.csv  (true donor-level; from 01_extract_seaad.R)
# Output: output/figures/Fig4_Clinical.png
# =============================================================================
source("R/figures/utils.R")

donor_file <- file.path(DATA_BIN, "donor_level_summary.csv")
if (!file.exists(donor_file)) stop("donor_level_summary.csv not found")
donor <- read.csv(donor_file)
donor <- donor[!is.na(donor$mean_cps), ]                 # drop reference donors -> n = 84

# ordinal staging codes (severity order)
roman <- c("Braak 0"=0,"Braak I"=1,"Braak II"=2,"Braak III"=3,"Braak IV"=4,"Braak V"=5,"Braak VI"=6)
cerad <- c("Absent"=0,"Sparse"=1,"Moderate"=2,"Frequent"=3)
abc   <- c("Not AD"=0,"Low"=1,"Intermediate"=2,"High"=3)
if (!"braak_num" %in% names(donor)) donor$braak_num <- roman[donor$braak]
if (!"cerad_num" %in% names(donor)) donor$cerad_num <- cerad[donor$cerad]
if (!"abc_num"   %in% names(donor)) donor$abc_num   <- abc[donor$abc]

# ── Panel A: staging correlations (grouped bars) ──────────────────────────────
sp <- function(g, s) suppressWarnings(cor(donor[[g]], donor[[s]],
                          method = "spearman", use = "complete.obs"))
stg <- c(Braak = "braak_num", CERAD = "cerad_num", ABC = "abc_num")
gene_cols <- c(ANLS = "ANLS", MCT4 = "MCT4", "V-ATPase (astro)" = "VATpase")
bar <- do.call(rbind, lapply(names(gene_cols), function(gn)
  data.frame(stage = names(stg),
             gene  = gn,
             rho   = sapply(stg, function(s) sp(gene_cols[[gn]], s)))))
bar$stage <- factor(bar$stage, levels = names(stg))
bar$gene  <- factor(bar$gene,  levels = names(gene_cols))
pA <- ggplot(bar, aes(stage, rho, fill = gene)) +
  geom_col(position = position_dodge(0.8), width = 0.7, color = "black", linewidth = 0.2) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  scale_fill_manual(values = c("ANLS" = "#E08214", "MCT4" = "#B2182B",
                               "V-ATPase (astro)" = "#2166AC"), name = NULL) +
  labs(x = NULL, y = "Spearman rho") +
  theme_paper + theme(legend.position = "top", axis.text.x = element_text(size = 10))

# ── Panel B: MCT4 vs CPS regression ───────────────────────────────────────────
fit <- lm(MCT4 ~ mean_cps, data = donor); r2 <- summary(fit)$r.squared
pval <- summary(fit)$coefficients["mean_cps", "Pr(>|t|)"]
pB <- ggplot(donor, aes(mean_cps, MCT4)) +
  geom_point(size = 2, color = "#B2182B", alpha = 0.7, stroke = 0.2) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 1) +
  annotate("text", x = -Inf, y = -Inf, hjust = -0.1, vjust = -0.6,
           label = sprintf("R^2 = %.3f\np = %.1e", r2, pval), size = 3.2) +
  labs(x = "CPS", y = "Astrocytic MCT4") +
  theme_paper

# ── Panel C: MCT4 by dementia status ──────────────────────────────────────────
donor$dem <- factor(donor$cognitive, levels = c("No dementia", "Dementia"))
dd <- donor[!is.na(donor$dem), ]
mw <- suppressWarnings(wilcox.test(MCT4 ~ dem, data = dd))
pC <- ggplot(dd, aes(dem, MCT4, fill = dem)) +
  geom_boxplot(width = 0.6, outlier.size = 0.8, alpha = 0.6) +
  scale_fill_manual(values = c("No dementia" = "#2166AC", "Dementia" = "#B2182B"),
                    guide = "none") +
  annotate("text", x = 1.5, y = Inf, vjust = 1.4,
           label = sprintf("Mann-Whitney  p = %.1e", mw$p.value), size = 3.2) +
  labs(x = NULL, y = "Astrocytic MCT4") +
  theme_paper

# ── Combine (uppercase bold panel tags) ───────────────────────────────────────
fig4 <- pA + pB + pC + plot_annotation(tag_levels = "A") &
        theme(plot.tag = element_text(size = 22, face = "bold"))
save_fig(fig4, "Fig4_Clinical.png", width = 15, height = 5)
