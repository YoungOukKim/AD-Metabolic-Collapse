# =============================================================================
# Fig1_Dissociation.R
#
# Figure 1 (a-d): Energetic dissociation defines the energy-starved lysosome
#
#   a — % change bar chart: ANLS vs V-ATPase vs Lysosomal genes
#   b — Normalized trajectory: MCT4, ANLS, V-ATPase (Bin 0.2 = 1.0)
#   c — Slope comparison: MCT4 vs ANLS vs V-ATPase (with 95% CI)
#   d — Energy/Demand ratio trajectory
#
# Input:  data/sample/astro_bin_means.csv
# Output: output/figures/Fig1_Dissociation.png
# =============================================================================
source("R/figures/utils.R")

# ── Load data ─────────────────────────────────────────────────────────────────
astro <- read.csv(file.path(DATA_BIN, "astro_bin_means.csv"))
astro <- add_composites_astro(astro)
astro <- astro[astro$bin >= 0.2, ]

# ── Panel a: % change bar chart ───────────────────────────────────────────────
pct_genes <- c("SLC16A3","HK2","LDHA","PDK1","SLC16A1","SLC2A1","PKM","PFKFB3")

pct_df <- data.frame()
for (g in c(pct_genes, "VATpase", "LAMP1", "CTSB", "CTSD")) {
  label    <- ifelse(g == "VATpase", "V-ATPase\n(10 sub.)", g)
  vals     <- astro[[g]]
  early    <- mean(vals[astro$bin %in% c(0.2, 0.3, 0.4)], na.rm = TRUE)
  late     <- mean(vals[astro$bin %in% c(0.6, 0.7, 0.8)], na.rm = TRUE)
  cat_type <- if (g %in% pct_genes) "ANLS" else if (g == "VATpase") "V-ATPase" else "Lysosomal"
  pct_df   <- rbind(pct_df, data.frame(Gene = label,
                                       pct  = (late - early) / early * 100,
                                       Category = cat_type))
}
pct_df$Gene <- factor(pct_df$Gene, levels = rev(pct_df$Gene))

fig1a <- ggplot(pct_df, aes(x = Gene, y = pct, fill = Category)) +
  geom_col(color = "black", linewidth = 0.3, width = 0.7) + coord_flip() +
  scale_fill_manual(values = c("ANLS"      = "#e74c3c",
                               "V-ATPase"  = "#3498db",
                               "Lysosomal" = "#95a5a6")) +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  labs(x = NULL, y = "% Change (Bins 0.2\u20130.4 \u2192 0.6\u20130.8)", title = "a") +
  theme_paper + theme(legend.position = "none")

# ── Panel b: Normalized trajectories ─────────────────────────────────────────
traj_df <- data.frame(
  bin     = astro$bin,
  MCT4    = norm_base(astro$MCT4,    astro$bin),
  ANLS    = norm_base(astro$ANLS,    astro$bin),
  VATpase = norm_base(astro$VATpase, astro$bin)
) |> tidyr::pivot_longer(-bin, names_to = "Marker", values_to = "norm_expr")

fig1b <- ggplot(traj_df, aes(x = bin, y = norm_expr, color = Marker)) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.8,
              linetype = "dashed", alpha = 0.15) +
  geom_line(linewidth = 1) + geom_point(size = 2.5) +
  scale_color_manual(values = c("MCT4"    = "#922B21",
                                "ANLS"    = "#e74c3c",
                                "VATpase" = "#3498db")) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "grey50") +
  scale_x_continuous(breaks = seq(0.2, 0.9, 0.1)) +
  labs(x = "Pseudo-progression Bin", y = "Normalized (Bin 0.2 = 1.0)",
       title = "b", color = NULL) +
  theme_paper + theme(legend.position = "right")

# ── Panel c: Slope comparison ─────────────────────────────────────────────────
bins_all <- seq(0.2, 0.9, 0.1)
s_mct4 <- get_slope(astro, "MCT4",    bins_all)
s_anls <- get_slope(astro, "ANLS",    bins_all)
s_vatp <- get_slope(astro, "VATpase", bins_all)

slopes <- data.frame(
  Marker = factor(c("MCT4", "ANLS", "V-ATPase"),
                  levels = c("MCT4", "ANLS", "V-ATPase")),
  slope  = c(s_mct4$beta, s_anls$beta, s_vatp$beta),
  se     = c(s_mct4$se,   s_anls$se,   s_vatp$se),
  sig    = c(sig_star(s_mct4$p), sig_star(s_anls$p), sig_star(s_vatp$p))
)

fig1c <- ggplot(slopes, aes(x = Marker, y = slope, fill = Marker)) +
  geom_col(color = "black", width = 0.6) +
  geom_errorbar(aes(ymin = slope - 1.96 * se, ymax = slope + 1.96 * se), width = 0.2) +
  geom_text(aes(label = sig, y = slope - 1.96 * se - 0.05),
            size = 4, fontface = "bold") +
  scale_fill_manual(values = c("MCT4"    = "#922B21",
                               "ANLS"    = "#e74c3c",
                               "V-ATPase"= "#3498db")) +
  geom_hline(yintercept = 0) +
  labs(x = NULL, y = "Normalized slope", title = "c") +
  theme_paper + theme(legend.position = "none")

# ── Panel d: Energy/Demand ratio ──────────────────────────────────────────────
fig1d <- ggplot(astro, aes(x = bin, y = ED_ratio)) +
  geom_smooth(method = "lm", se = TRUE, color = "#8e44ad",
              linetype = "dashed", alpha = 0.15) +
  geom_line(color = "#8e44ad", linewidth = 1.2) +
  geom_point(color = "#8e44ad", size = 3) +
  scale_x_continuous(breaks = seq(0.2, 0.9, 0.1)) +
  labs(x = "Pseudo-progression Bin", y = "Energy/Demand Ratio", title = "d") +
  theme_paper

# ── Combine & save ────────────────────────────────────────────────────────────
fig1 <- (fig1a + fig1b) / (fig1c + fig1d)
save_fig(fig1, "Fig1_Dissociation.png", width = 14, height = 10)
