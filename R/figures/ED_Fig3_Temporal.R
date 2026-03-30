# =============================================================================
# ED_Fig3_Temporal.R
#
# Extended Data Figure 3 (a-c): Temporal ordering of the metabolic cascade
#
#   a — Event timeline: sequential ordering of MCT4, LMR, LCN2, PTGDS events
#   b — Normalized trajectory: MCT4 / ANLS / V-ATPase with slope annotation
#   c — Astrocyte metabolic overload: EAAT2, ATP1A2, PTGDS, MCT4, MCT2 (neuron)
#
# Input:  data/sample/astro_bin_means.csv
#         data/sample/neuron_bin_means.csv  (for panel c MCT2)
# Output: output/figures/ED_Fig3_Temporal.png
# =============================================================================
source("R/figures/utils.R")

astro  <- read.csv(file.path(DATA_BIN, "astro_bin_means.csv"))
neuron <- read.csv(file.path(DATA_BIN, "neuron_bin_means.csv"))
astro  <- add_composites_astro(astro)
astro  <- astro[astro$bin >= 0.2,]
neuron <- neuron[neuron$bin >= 0.2,]

# ── Panel a: Event timeline ───────────────────────────────────────────────────
events <- data.frame(
  event    = c("MCT4 decline onset (10%)","LMR decline onset",
               "LCN2 first detected","PTGDS / EAAT2 / ATP1A2 peak",
               "PTGDS collapse","HMOX1 surge (+67%)"),
  bin      = c(0.3, 0.4, 0.5, 0.6, 0.7, 0.7),
  category = c("Energy","Lysosomal","Inflammation","Checkpoint","Checkpoint","Stress"),
  ypos     = c(6, 5, 4, 3, 2, 1)
)
cat_colors <- c(Energy="922B21",Lysosomal="#2980b9",
                Inflammation="#e74c3c",Checkpoint="#8e44ad",Stress="#e67e22")

ed3a <- ggplot(events, aes(x = bin, y = ypos)) +
  annotate("rect", xmin=0.45, xmax=0.65, ymin=0, ymax=7,
           fill="#FEF9E7", alpha=0.6) +
  annotate("text", x=0.55, y=6.7,
           label="Metabolic\ntransition zone", size=3, fontface="italic", color="grey40") +
  geom_segment(aes(x=0.2, xend=bin-0.01, y=ypos, yend=ypos, color=category),
               linewidth=1.5, alpha=0.4) +
  geom_point(aes(color=category), size=5) +
  geom_text(aes(label=sprintf("Bin %.1f", bin)), hjust=-0.3, size=3.5, fontface="bold") +
  scale_color_manual(values=cat_colors) +
  scale_x_continuous(limits=c(0.2, 0.85), breaks=seq(0.2, 0.8, 0.1)) +
  scale_y_continuous(breaks=events$ypos, labels=events$event, limits=c(0.5, 7.2)) +
  labs(x="Pseudo-progression Bin", y=NULL, title="a", color=NULL) +
  theme_paper + theme(axis.text.y=element_text(size=10),
                      panel.grid.major.y=element_line(color="grey90"))

# ── Panel b: Trajectory with slopes ──────────────────────────────────────────
traj_df <- data.frame(
  bin     = astro$bin,
  MCT4    = norm_base(astro$MCT4,    astro$bin),
  ANLS    = norm_base(astro$ANLS,    astro$bin),
  VATpase = norm_base(astro$VATpase, astro$bin)
) |> tidyr::pivot_longer(-bin, names_to="Marker", values_to="norm_expr")

ed3b <- ggplot(traj_df, aes(x=bin, y=norm_expr, color=Marker)) +
  geom_smooth(method="lm", se=TRUE, linewidth=0.8, linetype="dashed", alpha=0.15) +
  geom_line(linewidth=1) + geom_point(size=2.5) +
  scale_color_manual(values=c("MCT4"="#922B21","ANLS"="#e74c3c","VATpase"="#3498db")) +
  geom_hline(yintercept=1, linetype="dotted", color="grey50") +
  scale_x_continuous(breaks=seq(0.2,0.9,0.1)) +
  annotate("text", x=0.7, y=0.45, label="\u0394slope p = 0.0005",
           size=3.5, fontface="bold", color="#922B21") +
  labs(x="Pseudo-progression Bin", y="Normalized (Bin 0.2 = 1.0)",
       title="b", color=NULL) +
  theme_paper + theme(legend.position="right")

# ── Panel c: Metabolic overload ───────────────────────────────────────────────
overload_genes_astro  <- intersect(c("SLC1A2","ATP1A2","PTGDS","SLC16A3"), names(astro))
overload_genes_neuron <- intersect(c("SLC16A7"), names(neuron))

if (length(overload_genes_astro) > 0) {
  ol_astro <- astro[, c("bin", overload_genes_astro)]
  ol_neuron <- neuron[, c("bin", overload_genes_neuron)]
  ol_merged <- merge(ol_astro, ol_neuron, by="bin", all=TRUE)

  # Normalize to Bin 0.2
  for (col in setdiff(names(ol_merged), "bin")) {
    base <- ol_merged[[col]][ol_merged$bin == 0.2]
    if (!is.na(base) && base > 0)
      ol_merged[[col]] <- ol_merged[[col]] / base
  }

  rename_map <- c(SLC1A2="EAAT2 (SLC1A2)", ATP1A2="ATP1A2 (Na/K pump)",
                  PTGDS="PTGDS", SLC16A3="MCT4 (SLC16A3)",
                  SLC16A7="MCT2 neuron (SLC16A7)")
  names(ol_merged) <- ifelse(names(ol_merged) %in% names(rename_map),
                             rename_map[names(ol_merged)], names(ol_merged))

  ol_long <- ol_merged |> tidyr::pivot_longer(-bin, names_to="Gene", values_to="norm_expr")
  ol_colors <- c("EAAT2 (SLC1A2)"="#E67E22","ATP1A2 (Na/K pump)"="#D35400",
                 "PTGDS"="#8E44AD","MCT4 (SLC16A3)"="#922B21",
                 "MCT2 neuron (SLC16A7)"="#2980B9")

  ed3c <- ggplot(ol_long, aes(x=bin, y=norm_expr, color=Gene)) +
    annotate("rect", xmin=0.45, xmax=0.65, ymin=-Inf, ymax=Inf, fill="#FEF9E7", alpha=0.4) +
    geom_line(linewidth=1) + geom_point(size=2.5) +
    geom_hline(yintercept=1, linetype="dotted", color="grey50") +
    scale_color_manual(values=ol_colors) +
    scale_x_continuous(breaks=seq(0.2,0.9,0.1)) +
    labs(x="Pseudo-progression Bin", y="Normalized (Bin 0.2 = 1.0)",
         title="c", color=NULL) +
    theme_paper + theme(legend.position="right", legend.text=element_text(size=8))

  ed_fig3 <- ed3a / ed3b / ed3c + patchwork::plot_layout(heights=c(1,1,1))
} else {
  ed_fig3 <- ed3a / ed3b + patchwork::plot_layout(heights=c(1,1))
}

save_fig(ed_fig3, "ED_Fig3_Temporal.png", width=12, height=14)
