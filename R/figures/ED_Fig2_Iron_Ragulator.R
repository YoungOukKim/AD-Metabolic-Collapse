# =============================================================================
# ED_Fig2_Iron_Ragulator.R
#
# Extended Data Figure 2 (a-d): Iron, Ragulator, NHE6/PTGDS dynamics
#
#   a — TFRC-ANLS coupling scatter (bin-level, r = +0.666)
#   b — Iron gene trajectories (TFRC, FTH1, FTL, CP)
#   c — Ragulator subunit trajectories (LAMTOR1-5)
#   d — NHE6 and PTGDS dynamics (compensatory pH regulation)
#
# Input:  data/sample/astro_bin_means.csv
# Output: output/figures/ED_Fig2_Iron_Ragulator.png
# =============================================================================
source("R/figures/utils.R")

astro <- read.csv(file.path(DATA_BIN, "astro_bin_means.csv"))
astro <- add_composites_astro(astro)
astro <- astro[astro$bin >= 0.2,]

astro_bin <- astro |>
  dplyr::group_by(bin) |>
  dplyr::summarise(ANLS = mean(ANLS, na.rm=TRUE),
                   TFRC = mean(TFRC, na.rm=TRUE), .groups="drop")

r_tfrc <- cor.test(astro_bin$ANLS, astro_bin$TFRC)

# ── Panel a ───────────────────────────────────────────────────────────────────
ed2a <- ggplot(astro_bin, aes(x = ANLS, y = TFRC)) +
  geom_smooth(method = "lm", se = TRUE, color = "grey60", linetype = "dashed") +
  geom_point(aes(fill = bin), shape = 21, size = 4) +
  geom_text(aes(label = sprintf("%.1f", bin)), vjust = -1.2, size = 2.8) +
  scale_fill_gradient2(low="#27ae60", mid="#f39c12", high="#c0392b", midpoint=0.5) +
  annotate("text", x = min(astro_bin$ANLS) * 1.02,
           y = max(astro_bin$TFRC) * 0.95,
           label = sprintf("r = %+.3f\nn = %d bins", r_tfrc$estimate, nrow(astro_bin)),
           size = 3.5, fontface = "bold", hjust = 0) +
  labs(x = "ANLS composite (bin mean)", y = "TFRC (bin mean)", title = "a") +
  theme_paper + theme(legend.position = "none")

# ── Panel b: Iron trajectories ────────────────────────────────────────────────
iron_genes <- intersect(c("CP","FTH1","FTL","TFRC"), names(astro))
df_iron <- astro[, c("bin", iron_genes)] |>
  tidyr::pivot_longer(-bin, names_to = "Gene", values_to = "expr")

ed2b <- ggplot(df_iron, aes(x = bin, y = expr, color = Gene)) +
  geom_line(linewidth = 1) + geom_point(size = 2) +
  scale_x_continuous(breaks = seq(0.2, 0.9, 0.1)) +
  labs(x = "Bin", y = "Expression", title = "b") +
  theme_paper + theme(legend.position = "right")

# ── Panel c: Ragulator ────────────────────────────────────────────────────────
rag_genes <- intersect(paste0("LAMTOR",1:5), names(astro))
if (length(rag_genes) > 0) {
  df_rag <- astro[, c("bin", rag_genes)] |> tidyr::pivot_longer(-bin)
  ed2c <- ggplot(df_rag, aes(x = bin, y = value, color = name)) +
    geom_line(linewidth = 1) + geom_point(size = 2) +
    scale_x_continuous(breaks = seq(0.2, 0.9, 0.1)) +
    labs(x = "Bin", y = "Expression", title = "c") +
    theme_paper + theme(legend.position = "right")
} else {
  ed2c <- ggplot() + theme_void() +
    annotate("text", x=.5, y=.5, label="c: LAMTOR genes not in dataset")
}

# ── Panel d: NHE6 + PTGDS ────────────────────────────────────────────────────
ed2d <- ggplot(astro, aes(x = bin)) +
  geom_line(aes(y = SLC9A6, color = "NHE6"), linewidth = 1.2) +
  geom_point(aes(y = SLC9A6, color = "NHE6"), size = 2.5) +
  geom_line(aes(y = PTGDS,  color = "PTGDS"), linewidth = 1.2) +
  geom_point(aes(y = PTGDS, color = "PTGDS"), size = 2.5) +
  geom_vline(xintercept = 0.6, linetype = "dashed", color = "grey40") +
  annotate("text", x = 0.62, y = max(astro$PTGDS, na.rm=TRUE) * 0.95,
           label = "PTGDS peak", size = 3, fontface = "italic", color = "grey40") +
  scale_color_manual(values = c("NHE6"="#e67e22","PTGDS"="#2980b9")) +
  scale_x_continuous(breaks = seq(0.2, 0.9, 0.1)) +
  labs(x = "Bin", y = "Expression", title = "d", color = NULL) +
  theme_paper + theme(legend.position = "right")

ed_fig2 <- (ed2a + ed2b) / (ed2c + ed2d)
save_fig(ed_fig2, "ED_Fig2_Iron_Ragulator.png", width = 14, height = 10)
