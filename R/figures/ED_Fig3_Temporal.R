# =============================================================================
# ED_Fig3_Temporal.R  ->  Supplemental Figure 2
#
# Temporal ordering of the metabolic cascade (matches the published figure):
#   A  Event timeline: MCT4 -> LMR -> iron -> PTGDS (HMOX1 event removed,
#      consistent with the companion reappraisal that treats HMOX1 as a
#      detection artifact)
#   B  Lysosomal Metabolic Reserve (LMR) composite trajectory (single line;
#      decline precedes the iron-gene change-point by one bin)
#   C  Astrocyte metabolic overload: EAAT2, ATP1A2, PTGDS, MCT4, MCT2 (neuron)
#
# Layout: A on top (full width), B and C on the bottom row.
# Paths are RELATIVE to repo root. Run with working dir = repo root:
#   setwd("path/to/AD-Metabolic-Collapse"); source("R/figures/ED_Fig3_Temporal.R")
#
# Input:  data/sample/astro_bin_means.csv , data/sample/neuron_bin_means.csv
# Output: output/figures/Supplemental_Figure2.png  (300 dpi)
#         output/figures/Supplemental_Figure2.tif  (600 dpi, LZW)
# =============================================================================
source("R/figures/utils.R")

# -----------------------------------------------------------------------------
# LMR (Lysosomal Metabolic Reserve) — canonical 34-gene module is defined in
# utils.R (LMR_GENES / add_lmr) as the five functional sets tabulated in
# Supporting_Data_Values.xlsx. No editing required; the composite (rowMeans,
# normalized to Bin 0.2) reproduces the published panel B trough at Bin 0.5.
# -----------------------------------------------------------------------------

astro  <- add_composites_astro(read.csv(file.path(DATA_BIN, "astro_bin_means.csv")))
neuron <- read.csv(file.path(DATA_BIN, "neuron_bin_means.csv"))
astro  <- astro[astro$bin >= 0.2, ]
neuron <- neuron[neuron$bin >= 0.2, ]

# Note: astrocytic LCN2 mRNA is near-undetected in this snRNA-seq; the Bin-0.5
# event is the iron-gene co-decline (FTH1/FTL 10% onset), not LCN2 induction.
# CSF protein LCN2 (Table 2 ANCOVA) is a separate, valid result.

# ── A: Event timeline (HMOX1 removed) ─────────────────────────────────────────
events <- data.frame(
  event    = c("MCT4 onset","LMR onset","Iron co-decline",
               "Compensatory peak","PTGDS collapse"),
  bin      = c(0.3, 0.4, 0.5, 0.6, 0.7),
  category = c("Energy","Lysosomal","Iron","Checkpoint","Checkpoint"),
  blab     = c("Bin 0.3 (-10%)","Bin 0.4","Bin 0.5","Bin 0.6","Bin 0.7"),
  ypos     = c(5, 4, 3, 2, 1)
)
cat_colors <- c(Energy="#922B21",Lysosomal="#2980b9",Iron="#16a085",Checkpoint="#8e44ad")

p_a <- ggplot(events, aes(x = bin, y = ypos)) +
  annotate("rect", xmin=0.45, xmax=0.65, ymin=0, ymax=6, fill="#FEF9E7", alpha=0.6) +
  annotate("text", x=0.55, y=5.7, label="Metabolic\ntransition zone",
           size=3, fontface="italic", color="grey40") +
  geom_segment(aes(x=0.2, xend=bin-0.01, y=ypos, yend=ypos, color=category),
               linewidth=1.5, alpha=0.4) +
  geom_point(aes(color=category), size=5) +
  geom_text(aes(label=blab), hjust=-0.15, size=3.5, fontface="bold") +
  scale_color_manual(values=cat_colors) +
  scale_x_continuous(limits=c(0.2, 0.86), breaks=seq(0.2, 0.8, 0.1)) +
  scale_y_continuous(breaks=events$ypos, labels=events$event, limits=c(0.5, 6.2)) +
  labs(x="Pseudo-progression Bin", y=NULL, title="A", color=NULL) +
  theme_paper + theme(axis.text.y=element_text(size=10),
                      panel.grid.major.y=element_line(color="grey90"))

# ── B: LMR composite trajectory ───────────────────────────────────────────────
astro  <- add_lmr(astro)
lmr_df <- data.frame(bin = astro$bin, LMR = norm_base(astro$LMR, astro$bin))

# Built-in reproduction check: the published panel B dips at Bin 0.5 (LMR decline
# leads the iron co-decline by one bin). Prints the trough bin and PASS/FAIL; if
# it FAILs, verify the input CSV rather than editing this check.
.lmr_o <- lmr_df[order(lmr_df$bin), ]
.trough <- .lmr_o$bin[which.min(.lmr_o$LMR)]
cat(sprintf("[LMR check] trough (dip) bin = %s  (published expectation: 0.5)  -> %s\n",
            format(.trough), if (isTRUE(all.equal(.trough, 0.5))) "PASS" else "FAIL"))
if (!isTRUE(all.equal(.trough, 0.5)))
  warning("LMR trajectory does not dip at Bin 0.5; check data/sample/astro_bin_means.csv.",
          call. = FALSE)

p_b <- ggplot(lmr_df, aes(x=bin, y=LMR)) +
  geom_hline(yintercept=1, linetype="dotted", color="grey50") +
  geom_smooth(method="lm", se=FALSE, linetype="dashed", color="#2980b9", linewidth=0.9) +
  geom_line(color="#2980b9", linewidth=1) + geom_point(color="#2980b9", size=2.5) +
  geom_vline(xintercept=0.4, linetype="dashed", color="#2980b9", alpha=0.8) +
  geom_vline(xintercept=0.5, linetype="dashed", color="#16a085", alpha=0.8) +
  annotate("text", x=0.4, y=max(lmr_df$LMR), label="LMR onset\n(Bin 0.4)",
           color="#2980b9", size=2.8, vjust=0) +
  annotate("text", x=0.5, y=max(lmr_df$LMR), label="Iron co-decline\n(Bin 0.5)",
           color="#16a085", size=2.8, vjust=0) +
  scale_x_continuous(breaks=seq(0.2,0.9,0.1)) +
  labs(x="Pseudo-progression Bin", y="Normalized (Bin 0.2 = 1.0)", title="B") +
  theme_paper

# ── C: Metabolic overload ─────────────────────────────────────────────────────
ol_a <- intersect(c("SLC1A2","ATP1A2","PTGDS","SLC16A3"), names(astro))
ol_n <- intersect(c("SLC16A7"), names(neuron))
ol   <- merge(astro[, c("bin", ol_a)], neuron[, c("bin", ol_n)], by="bin", all=TRUE)
for (col in setdiff(names(ol), "bin")) {
  base <- ol[[col]][ol$bin == 0.2]
  if (!is.na(base) && base > 0) ol[[col]] <- ol[[col]] / base
}
rename_map <- c(SLC1A2="EAAT2 (SLC1A2)", ATP1A2="ATP1A2 (Na/K pump)", PTGDS="PTGDS",
                SLC16A3="MCT4 (SLC16A3)", SLC16A7="MCT2 neuron (SLC16A7)")
names(ol) <- ifelse(names(ol) %in% names(rename_map), rename_map[names(ol)], names(ol))
ol_long <- ol |> tidyr::pivot_longer(-bin, names_to="Gene", values_to="norm_expr")
ol_colors <- c("EAAT2 (SLC1A2)"="#E67E22","ATP1A2 (Na/K pump)"="#D35400","PTGDS"="#8E44AD",
               "MCT4 (SLC16A3)"="#922B21","MCT2 neuron (SLC16A7)"="#2980B9")

p_c <- ggplot(ol_long, aes(x=bin, y=norm_expr, color=Gene)) +
  annotate("rect", xmin=0.45, xmax=0.65, ymin=-Inf, ymax=Inf, fill="#FEF9E7", alpha=0.4) +
  geom_hline(yintercept=1, linetype="dotted", color="grey50") +
  geom_line(linewidth=1) + geom_point(size=2.5) +
  scale_color_manual(values=ol_colors) +
  scale_x_continuous(breaks=seq(0.2,0.9,0.1)) +
  labs(x="Pseudo-progression Bin", y="Normalized (Bin 0.2 = 1.0)", title="C", color=NULL) +
  theme_paper + theme(legend.position="right", legend.text=element_text(size=8))

# ── Assemble: A on top, B | C on the bottom ───────────────────────────────────
supp_fig2 <- p_a / (p_b | p_c) + patchwork::plot_layout(heights=c(1, 1))

ggsave(file.path(FIG_OUT, "Supplemental_Figure2.png"), supp_fig2,
       width=13, height=10, dpi=300, bg="white")
ggsave(file.path(FIG_OUT, "Supplemental_Figure2.tif"), supp_fig2,
       width=13, height=10, dpi=600, bg="white", device="tiff", compression="lzw")
message("Saved: ", file.path(FIG_OUT, "Supplemental_Figure2.{png,tif}"))
