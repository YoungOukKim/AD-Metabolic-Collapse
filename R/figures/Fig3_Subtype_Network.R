# =============================================================================
# Fig3_Subtype_Network.R
#
# Figure 3 (a-c): Subtype consistency and network architecture
#
#   a — ANLS trajectory by astrocyte supertype (6 subtypes)
#   b — Subtype proportion shift across pseudo-progression
#   c — Network degree centrality (top hubs)
#
# Input:  data/sample/astro_subtype_trajectories.csv
# Output: output/figures/Fig3_Subtype_Network.png
# =============================================================================
source("R/figures/utils.R")

sub_file <- file.path(DATA_BIN, "astro_subtype_trajectories.csv")

if (file.exists(sub_file)) {
  sub_df <- read.csv(sub_file)

  # ── Panel a: ANLS by subtype ──────────────────────────────────────────────
  fig3a <- ggplot(sub_df, aes(x = bin, y = ANLS, color = subtype)) +
    geom_line(linewidth = 1) + geom_point(size = 2) +
    scale_x_continuous(breaks = seq(0.1, 0.9, 0.1)) +
    labs(x = "Pseudo-progression Bin", y = "ANLS composite",
         title = "a", color = "Subtype") +
    theme_paper + theme(legend.position = "right")

  # ── Panel b: Proportion shift ─────────────────────────────────────────────
  prop_df <- sub_df |>
    dplyr::group_by(bin) |>
    dplyr::mutate(total = sum(n_cells), proportion = n_cells / total) |>
    dplyr::ungroup()

  fig3b <- ggplot(prop_df, aes(x = bin, y = proportion, fill = subtype)) +
    geom_area(alpha = 0.8, color = "white", linewidth = 0.3) +
    scale_x_continuous(breaks = seq(0.1, 0.9, 0.1)) +
    labs(x = "Pseudo-progression Bin", y = "Proportion",
         title = "b", fill = "Subtype") +
    theme_paper + theme(legend.position = "right")
} else {
  message("astro_subtype_trajectories.csv not found — using placeholder panels")
  fig3a <- fig3b <- ggplot() + theme_void() +
    annotate("text", x=.5, y=.5,
             label = "Run 01_extract_seaad.R to generate\nastro_subtype_trajectories.csv")
}

# ── Panel c: Network degree ───────────────────────────────────────────────────
# Pre-computed from correlation network (|r| > 0.7, 28 genes, 128 edges)
net_df <- data.frame(
  Gene     = c("TFRC","VDAC1","GLUT1","MCT1","LDHA","FTH1","PKM",
               "MCT4","LAMP1","LAMP2","HK1","V-ATPase","NHE6"),
  degree   = c(16, 15, 15, 13, 12, 11, 10, 9, 8, 7, 7, 5, 4),
  Category = c("Iron","Mito","ANLS","ANLS","ANLS","Iron","ANLS",
               "ANLS","Lyso","Lyso","ANLS","V-ATPase","pH")
)
net_df$Gene <- factor(net_df$Gene, levels = net_df$Gene[order(net_df$degree)])

fig3c <- ggplot(net_df, aes(x = Gene, y = degree, fill = Category)) +
  geom_col(color = "black", linewidth = 0.3, width = 0.7) + coord_flip() +
  geom_hline(yintercept = 14, linetype = "dashed", color = "#e74c3c") +
  annotate("text", x = 1.5, y = 14.3, label = "80th percentile",
           size = 3, color = "#e74c3c", fontface = "italic") +
  scale_fill_manual(values = c("ANLS"     = "#e74c3c",
                               "Iron"     = "#e67e22",
                               "Mito"     = "#27ae60",
                               "Lyso"     = "#95a5a6",
                               "V-ATPase" = "#3498db",
                               "pH"       = "#95a5a6")) +
  labs(x = NULL, y = "Network degree (|r| > 0.7)", title = "c") +
  theme_paper + theme(legend.position = "right")

# ── Combine & save ────────────────────────────────────────────────────────────
fig3 <- (fig3a + fig3b) / fig3c
save_fig(fig3, "Fig3_Subtype_Network.png", width = 14, height = 10)
