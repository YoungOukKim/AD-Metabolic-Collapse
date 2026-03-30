# =============================================================================
# Fig2_CrossCellular.R
#
# Figure 2 (a-f): Cross-cellular energy decline propagates to neuronal lysosomes
#
#   a — Cross-cellular heatmap (zero-order r + partial r, CPS-adjusted)
#   b — MCT4 vs Neuron V-ATPase scatter (bin-level)
#   c — Partial correlation scatter (CPS residuals)
#   d — Zero-order vs Partial r barplot (4 pairs)
#   e — Two-layer model schematic
#   f — V-ATPase % change: astrocyte vs neuron comparison
#
# Input:  data/sample/astro_bin_means.csv
#         data/sample/neuron_bin_means.csv
# Output: output/figures/Fig2_CrossCellular.png
# =============================================================================
source("R/figures/utils.R")

# ── Load & prepare ────────────────────────────────────────────────────────────
astro  <- read.csv(file.path(DATA_BIN, "astro_bin_means.csv"))
neuron <- read.csv(file.path(DATA_BIN, "neuron_bin_means.csv"))
astro  <- add_composites_astro(astro)
neuron <- add_composites_neuron(neuron)
astro  <- astro[astro$bin   >= 0.2, ]
neuron <- neuron[neuron$bin >= 0.2, ]

df <- data.frame(
  bin     = astro$bin,
  MCT4    = astro$MCT4,
  ANLS    = astro$ANLS,
  VATPase = neuron$VATPase,
  LDHB    = neuron$LDHB,
  LAMP1   = neuron$LAMP1
)

# ── Panel a: Heatmap ──────────────────────────────────────────────────────────
hm <- expand.grid(Astro  = c("MCT4", "ANLS"),
                  Neuron = c("VATPase", "LDHB", "LAMP1"),
                  stringsAsFactors = FALSE)
hm$r_zero <- hm$r_partial <- hm$p_partial <- NA

for (i in seq_len(nrow(hm))) {
  ct          <- cor.test(df[[hm$Astro[i]]], df[[hm$Neuron[i]]])
  pc          <- partial_cor(df[[hm$Astro[i]]], df[[hm$Neuron[i]]], df$bin)
  hm$r_zero[i]    <- ct$estimate
  hm$r_partial[i] <- pc$r
  hm$p_partial[i] <- pc$p
}
hm$sig   <- sapply(hm$p_partial, sig_star)
hm$label <- sprintf("r = %+.3f\npr = %+.3f %s", hm$r_zero, hm$r_partial, hm$sig)
hm$Astro  <- factor(hm$Astro,  levels = c("ANLS", "MCT4"))
hm$Neuron <- factor(hm$Neuron, levels = c("VATPase", "LDHB", "LAMP1"))

fig2a <- ggplot(hm, aes(x = Neuron, y = Astro, fill = r_zero)) +
  geom_tile(color = "white", linewidth = 2) +
  geom_text(aes(label = label), size = 3.2, fontface = "bold") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-0.2, 1), name = "Zero-order r") +
  labs(title = "a") +
  theme_paper + theme(axis.title = element_blank(),
                      axis.text  = element_text(size = 10, face = "bold"),
                      legend.position = "right")

# ── Panel b: Scatter ─────────────────────────────────────────────────────────
r_b  <- cor.test(df$MCT4, df$VATPase)
pc_b <- partial_cor(df$MCT4, df$VATPase, df$bin)

fig2b <- ggplot(df, aes(x = MCT4, y = VATPase)) +
  geom_smooth(method = "lm", se = TRUE, color = "grey60",
              linetype = "dashed", linewidth = 0.8) +
  geom_point(aes(fill = bin), shape = 21, size = 5, color = "black", stroke = 0.8) +
  geom_text(aes(label = sprintf("%.1f", bin)), vjust = -1.2, size = 3) +
  scale_fill_gradient2(low = "#27ae60", mid = "#f39c12", high = "#c0392b",
                       midpoint = 0.5, name = "Bin") +
  annotate("text", x = min(df$MCT4) * 1.02, y = max(df$VATPase) * 0.98,
           label = sprintf("r = %+.3f (p = %.3f)\npr = %+.3f (p = %.3f)",
                           r_b$estimate, r_b$p.value, pc_b$r, pc_b$p),
           size = 3, fontface = "bold", hjust = 0) +
  labs(x = "Astrocytic MCT4 (SLC16A3)",
       y = "Neuronal V-ATPase Composite", title = "b") +
  theme_paper + theme(legend.position = "right")

# ── Panel c: Partial scatter ──────────────────────────────────────────────────
pc_c <- partial_cor(df$MCT4, df$VATPase, df$bin)
resid_df <- data.frame(bin = df$bin, MCT4_res = pc_c$rx, VATPase_res = pc_c$ry)

fig2c <- ggplot(resid_df, aes(x = MCT4_res, y = VATPase_res)) +
  geom_smooth(method = "lm", se = TRUE, color = "#B2182B",
              linetype = "dashed", linewidth = 0.8, alpha = 0.15) +
  geom_point(aes(fill = bin), shape = 21, size = 5, color = "black", stroke = 0.8) +
  geom_text(aes(label = sprintf("%.1f", bin)), vjust = -1.2, size = 3) +
  scale_fill_gradient2(low = "#27ae60", mid = "#f39c12", high = "#c0392b",
                       midpoint = 0.5, name = "Bin") +
  annotate("text", x = min(resid_df$MCT4_res) * 0.9,
           y = max(resid_df$VATPase_res) * 0.9,
           label = sprintf("partial r = %+.3f\np = %.3f",
                           pc_c$r, pc_c$p),
           size = 3, fontface = "bold", hjust = 0) +
  labs(x = "MCT4 residual (CPS-adjusted)",
       y = "Neuron V-ATPase residual (CPS-adjusted)", title = "c") +
  theme_paper + theme(legend.position = "right")

# ── Panel d: Zero vs Partial barplot ──────────────────────────────────────────
pairs_info <- list(
  list(a = "MCT4", n = "VATPase", label = "MCT4\n\u2194VATPase"),
  list(a = "MCT4", n = "LDHB",    label = "MCT4\n\u2194LDHB"),
  list(a = "MCT4", n = "LAMP1",   label = "MCT4\n\u2194LAMP1"),
  list(a = "ANLS", n = "VATPase", label = "ANLS\n\u2194VATPase")
)
bar_data <- data.frame()
for (pr in pairs_info) {
  r0 <- cor.test(df[[pr$a]], df[[pr$n]])
  pc <- partial_cor(df[[pr$a]], df[[pr$n]], df$bin)
  bar_data <- rbind(bar_data, data.frame(pair      = pr$label,
                                         r_zero    = r0$estimate,
                                         r_partial = pc$r,
                                         p_partial = pc$p))
}
bar_data$pair <- factor(bar_data$pair, levels = bar_data$pair)
bar_long <- bar_data |>
  tidyr::pivot_longer(c(r_zero, r_partial), names_to = "type", values_to = "r") |>
  dplyr::mutate(type = ifelse(type == "r_zero", "Zero-order", "Partial (CPS-adj)"),
                type = factor(type, levels = c("Zero-order", "Partial (CPS-adj)")))
bar_long <- merge(bar_long, bar_data[, c("pair","p_partial")], by = "pair")
bar_long$sig <- ifelse(bar_long$type == "Partial (CPS-adj)",
                       sapply(bar_long$p_partial, sig_star), "")

fig2d <- ggplot(bar_long, aes(x = pair, y = r, fill = type)) +
  geom_col(position = position_dodge(0.8), width = 0.7,
           color = "black", linewidth = 0.3) +
  geom_text(aes(label = sig), position = position_dodge(0.8),
            vjust = -0.3, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = c("Zero-order"      = "#4393C3",
                               "Partial (CPS-adj)" = "#D6604D"),
                    name = NULL) +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  coord_cartesian(ylim = c(-0.1, 1.02)) +
  labs(x = NULL, y = "Pearson r", title = "d") +
  theme_paper + theme(legend.position = c(0.78, 0.92),
                      axis.text.x = element_text(size = 9, face = "bold"))

# ── Panel e: Two-layer schematic ──────────────────────────────────────────────
fig2e <- ggplot() + xlim(0, 10) + ylim(0, 10) + theme_void() +
  annotate("text", x = 5, y = 9.7, label = "e", size = 5, fontface = "bold") +
  annotate("rect",  xmin=0.3, xmax=9.7, ymin=5.5, ymax=9.2,
           fill="#FADBD8", color="#C0392B", linewidth=1.2, alpha=0.3) +
  annotate("text", x=5, y=8.7, label="ASTROCYTE (Layer 1)",
           size=4, fontface="bold", color="#922B21") +
  annotate("text", x=2.8, y=7.3,
           label="ANLS Disruption\nMCT4 \u221243%\nHK2 \u221235%",
           size=3, color="#C0392B") +
  annotate("text", x=7.2, y=7.3,
           label="Astro lysosome\nATP deprived\n(V-ATPase: \u22120.8%)",
           size=3, color="#C0392B") +
  annotate("segment", x=4.5, xend=5.5, y=7.3, yend=7.3,
           arrow=arrow(length=unit(0.25,"cm")), linewidth=1.2, color="#C0392B") +
  annotate("segment", x=2.8, xend=2.8, y=5.5, yend=4.3,
           arrow=arrow(length=unit(0.3,"cm")), linewidth=1.8, color="#E74C3C") +
  annotate("text", x=4.5, y=4.9,
           label="Lactate delivery \u2193\n(partial r = +0.856**)",
           size=3, fontface="bold", color="#E74C3C") +
  annotate("rect", xmin=0.3, xmax=9.7, ymin=1.3, ymax=4.0,
           fill="#EBF5FB", color="#2166AC", linewidth=1.2, alpha=0.3) +
  annotate("text", x=5, y=3.6, label="NEURON (Layer 2)",
           size=4, fontface="bold", color="#1A5276") +
  annotate("text", x=2.8, y=2.5,
           label="LDHB: Lactate \u2192 ATP\n(fuel source lost)",
           size=3, color="#2166AC") +
  annotate("text", x=7.2, y=2.5,
           label="Neuron lysosome\nV-ATPase unfueled\n(\u22125.4%, 7\u00D7 > astro)",
           size=3, color="#2166AC") +
  annotate("segment", x=4.5, xend=5.5, y=2.5, yend=2.5,
           arrow=arrow(length=unit(0.25,"cm")), linewidth=1.2, color="#2166AC") +
  annotate("rect", xmin=3, xmax=7, ymin=0.1, ymax=1.0,
           fill="#F9EBEA", color="#922B21", linewidth=1.5) +
  annotate("text", x=5, y=0.55, label="\u2192 PANTHOS",
           size=3.8, fontface="bold", color="#922B21") +
  annotate("segment", x=5, xend=5, y=1.3, yend=1.05,
           arrow=arrow(length=unit(0.25,"cm")), linewidth=1.2, color="#922B21")

# ── Panel f: V-ATPase decline comparison ─────────────────────────────────────
pct_a <- (mean(astro$VATpase[astro$bin %in% c(0.6,0.7,0.8)]) -
          mean(astro$VATpase[astro$bin %in% c(0.2,0.3,0.4)])) /
         mean(astro$VATpase[astro$bin %in% c(0.2,0.3,0.4)]) * 100
pct_n <- (mean(df$VATPase[df$bin %in% c(0.6,0.7,0.8)]) -
          mean(df$VATPase[df$bin %in% c(0.2,0.3,0.4)])) /
         mean(df$VATPase[df$bin %in% c(0.2,0.3,0.4)]) * 100

bar_df <- data.frame(Cell = factor(c("Astrocyte","Neuron"),
                                   levels = c("Astrocyte","Neuron")),
                     pct = c(pct_a, pct_n))

fig2f <- ggplot(bar_df, aes(x = Cell, y = pct, fill = Cell)) +
  geom_col(color = "black", width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%", pct)),
            vjust = ifelse(bar_df$pct < -1, 1.5, -0.5),
            size = 4, fontface = "bold") +
  scale_fill_manual(values = c("Astrocyte" = "#3498db", "Neuron" = "#e74c3c")) +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  annotate("segment", x=1, xend=2, y=0.3, yend=0.3, linewidth=0.4, color="grey30") +
  annotate("text", x=1.5, y=0.5,
           label = sprintf("%.0f\u00D7", abs(pct_n / pct_a)),
           size = 3.5, fontface = "bold") +
  labs(x = NULL, y = "V-ATPase % Change\n(Bins 0.2\u20130.4 \u2192 0.6\u20130.8)",
       title = "f") +
  theme_paper + theme(legend.position = "none")

# ── Combine & save ────────────────────────────────────────────────────────────
fig2 <- (fig2a + fig2b + fig2c + plot_layout(ncol=3, widths=c(1,1.1,1.1))) /
        (fig2d + fig2e + fig2f + plot_layout(ncol=3, widths=c(1.1,1.2,0.7)))

save_fig(fig2, "Fig2_CrossCellular.png", width = 18, height = 12)
