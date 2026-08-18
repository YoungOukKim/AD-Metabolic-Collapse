# =============================================================================
# Fig2_CrossCellular.R  [CORRECTED  -  donor-level; matches published Figure 2]
#
# Figure 2 (A-C): Cross-cellular metabolic coupling at the DONOR level (n = 84).
#   A  -  Astrocytic MCT4 vs neuronal V-ATPase (zero-order, colored by CPS)
#   B  -  Zero-order vs CPS-adjusted partial correlations, five cross-cellular pairs
#   C  -  CPS-adjusted partial-correlation residual scatter (MCT4 vs neuronal V-ATPase)
#
# NOTE: This replaces the previous bin-level (n = 9) version. Bin-level trajectory
# correlations were progression-confounded pseudo-replicates and are not used.
#
# Input:  donor_level_summary.csv  (true donor-level; from 01_extract_seaad.R)
# Output: output/figures/Fig2_CrossCellular.png
# =============================================================================
source("R/figures/utils.R")

donor_file <- file.path(DATA_BIN, "donor_level_summary.csv")
if (!file.exists(donor_file)) stop("donor_level_summary.csv not found")
donor <- read.csv(donor_file)
donor <- donor[!is.na(donor$mean_cps), ]        # drop reference donors (no CPS) -> n = 84

partial_cor <- function(x, y, z) {
  ok <- complete.cases(x, y, z); x <- x[ok]; y <- y[ok]; z <- z[ok]
  rx <- residuals(lm(x ~ z)); ry <- residuals(lm(y ~ z))
  ct <- cor.test(rx, ry)
  list(r = unname(ct$estimate), p = ct$p.value, rx = rx, ry = ry)
}

# -- Panel A: zero-order scatter -----------------------------------------------
r0 <- cor(donor$MCT4, donor$VATPase_n, use = "complete.obs")
pA <- ggplot(donor, aes(MCT4, VATPase_n, color = mean_cps)) +
  geom_point(size = 2, stroke = 0.2) +
  geom_smooth(method = "lm", se = FALSE, color = "#B2182B", linewidth = 1) +
  scale_color_viridis_c(name = "CPS") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.4,
           label = sprintf("n = 84 donors\nr = +%.3f", r0), size = 3.2) +
  labs(x = "Astrocytic MCT4 (per-donor mean)",
       y = "Neuronal V-ATPase (per-donor mean)") +
  theme_paper + theme(legend.position = "right")

# -- Panel B: zero-order vs CPS-adjusted partial, five pairs -------------------
pairs <- list(c("MCT4","LAMP1_n"), c("MCT4","VATPase_n"), c("ANLS","MCT4"),
              c("MCT4","LDHB_n"),  c("ANLS","VATPase_n"))
labs5 <- c("MCT4~LAMP1","MCT4~V-ATPase","ANLS~MCT4","MCT4~LDHB","ANLS~V-ATPase")
zz <- sapply(pairs, function(p) cor(donor[[p[1]]], donor[[p[2]]], use = "complete.obs"))
pp <- sapply(pairs, function(p) partial_cor(donor[[p[1]]], donor[[p[2]]], donor$mean_cps)$r)
bardf <- data.frame(pair = factor(rep(labs5, 2), levels = labs5),
                    type = rep(c("zero-order","CPS-adjusted partial"), each = 5),
                    r    = c(zz, pp))
pB <- ggplot(bardf, aes(pair, r, fill = type)) +
  geom_col(position = position_dodge(0.8), width = 0.7, color = "black", linewidth = 0.2) +
  scale_fill_manual(values = c("zero-order" = "#999999",
                               "CPS-adjusted partial" = "#B2182B"), name = NULL) +
  coord_cartesian(ylim = c(0, 0.7)) +
  labs(x = NULL, y = "Pearson r") +
  theme_paper + theme(legend.position = "top",
                      axis.text.x = element_text(size = 12))

# -- Panel C: CPS-adjusted residual scatter ------------------------------------
pcr <- partial_cor(donor$MCT4, donor$VATPase_n, donor$mean_cps)
resdf <- data.frame(rx = pcr$rx, ry = pcr$ry)
pC <- ggplot(resdf, aes(rx, ry)) +
  geom_point(size = 2, color = "#B2182B", alpha = 0.65, stroke = 0.2) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 1) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.4,
           label = sprintf("partial r = +%.3f\np = %.1e", pcr$r, pcr$p), size = 3.2) +
  labs(x = "Astrocytic MCT4 residual (CPS-adjusted)",
       y = "Neuronal V-ATPase residual (CPS-adjusted)") +
  theme_paper

# -- Combine (uppercase bold panel tags, matching the rest of the paper) -------
pA <- pA + panel_tag("A")
pB <- pB + panel_tag("B")
pC <- pC + panel_tag("C")

fig2 <- pA + pB + pC
save_fig(fig2, "Fig2_CrossCellular.png", width = 15, height = 5)
