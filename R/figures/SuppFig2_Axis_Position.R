# =============================================================================
# SuppFig2_Axis_Position.R  ->  Supplemental Figure 2
#
# Position along the pseudo-progression axis. NOT a temporal ordering: the axis
# is cross-sectional, and the revision withdrew every within-individual timing
# claim. The panels describe where each trajectory sits on the axis.
#   A  Position of each trajectory on the axis: the bin at which 10% of the
#      total change from Bin 0.2 is reached. HMOX1 is not shown, consistent
#      with the companion reappraisal that treats it as a detection artifact.
#   B  Lysosomal Metabolic Reserve (LMR) composite trajectory (single line;
#      its minimum is at Bin 0.4)
#   C  Astrocyte metabolic overload: EAAT2, ATP1A2, PTGDS, MCT4, MCT2 (neuron)
#
# Layout: A on top (full width), B and C on the bottom row.
# Paths are RELATIVE to repo root. Run with working dir = repo root:
#   setwd("path/to/AD-Metabolic-Collapse"); source("R/figures/SuppFig2_Axis_Position.R")
#
# Input:  data/sample/astro_bin_means.csv , data/sample/neuron_bin_means.csv
# Output: output/figures/Supplemental_Figure2.png  (300 dpi)
#         output/figures/Supplemental_Figure2.tif  (600 dpi, LZW)
# =============================================================================
source("R/figures/utils.R")

# -----------------------------------------------------------------------------
# LMR (Lysosomal Metabolic Reserve)  -  canonical 34-gene module is defined in
# utils.R (LMR_GENES / add_lmr) as the five functional sets tabulated in
# Supporting_Data_Values.xlsx. No editing required; the composite (rowMeans,
# normalized to Bin 0.2) reproduces the published panel B, whose minimum is at
# Bin 0.4. Sixteen of the 34 genes are below the detection floor in this panel
# and drop out of the mean; that is the published behaviour, not a data fault.
# -----------------------------------------------------------------------------

astro  <- add_composites_astro(read.csv(file.path(DATA_BIN, "astro_bin_means.csv")))
neuron <- read.csv(file.path(DATA_BIN, "neuron_bin_means.csv"))
astro  <- astro[astro$bin >= 0.2, ]
neuron <- neuron[neuron$bin >= 0.2, ]

# Note: astrocytic LCN2 mRNA is near-undetected in this snRNA-seq; the Bin-0.5
# event is the iron-gene co-decline (FTH1/FTL 10% onset), not LCN2 induction.
# CSF protein LCN2 (Table 2 ANCOVA) is a separate, valid result.

# -- A: Event timeline (HMOX1 removed) -----------------------------------------
# Positions along the axis, not events in time. Every label that asserted an
# ordering or a compensatory phase was withdrawn in the revision: "onset",
# "compensatory peak" and "collapse" are gone, MCT4 is Bin 0.5 rather than the
# Bin 0.3 reported earlier, and Bin 0.3 belongs to TFRC.
events <- data.frame(
  event    = c("TFRC","LMR","MCT4","FTH1 / FTL",
               "PTGDS / EAAT2 peak","PTGDS"),
  bin      = c(0.3, 0.4, 0.5, 0.5, 0.6, 0.7),
  category = c("Iron","Lysosomal","Energy","Iron","Checkpoint","Checkpoint"),
  blab     = c("Bin 0.3","Bin 0.4","Bin 0.5","Bin 0.5","Bin 0.6","Bin 0.7"),
  ypos     = c(6, 5, 4, 3, 2, 1)
)
cat_colors <- c(Energy="#922B21",Lysosomal="#2980b9",Iron="#16a085",Checkpoint="#8e44ad")

# The shaded band and the "metabolic transition zone" label that used to sit
# here are removed. "Transition window" was withdrawn in the revision under
# reviewer 1 point 5: the axis is cross-sectional and marking a zone on it
# asserts a stage that the data do not establish.
p_a <- ggplot(events, aes(x = bin, y = ypos)) +
  geom_segment(aes(x=0.2, xend=bin-0.01, y=ypos, yend=ypos, color=category),
               linewidth=1.5, alpha=0.4) +
  geom_point(aes(color=category), size=5) +
  geom_text(aes(label=blab), hjust=-0.15, size=4.9, fontface="bold") +
  scale_color_manual(values=cat_colors) +
  scale_x_continuous(limits=c(0.2, 0.86), breaks=seq(0.2, 0.8, 0.1)) +
  scale_y_continuous(breaks=events$ypos, labels=events$event, limits=c(0.5, 6.2)) +
  labs(x="Pseudo-progression Bin", y=NULL, title="A", color=NULL) +
  theme_paper + theme(axis.text.y=element_text(size=14),
                      panel.grid.major.y=element_line(color="grey90"))

# -- B: LMR composite trajectory -----------------------------------------------
astro  <- add_lmr(astro)
lmr_df <- data.frame(bin = astro$bin, LMR = norm_base(astro$LMR, astro$bin))

# Built-in reproduction check: the published panel B has its minimum at Bin 0.4.
# An earlier version of this check expected Bin 0.5 and therefore failed against
# the deposited data; the expectation, not the data, was wrong. Prints the bin
# and PASS/FAIL; if it FAILs, verify the input CSV rather than editing this check.
.lmr_o <- lmr_df[order(lmr_df$bin), ]
.trough <- .lmr_o$bin[which.min(.lmr_o$LMR)]
cat(sprintf("[LMR check] minimum bin = %s  (published expectation: 0.4)  -> %s\n",
            format(.trough), if (isTRUE(all.equal(.trough, 0.4))) "PASS" else "FAIL"))
if (!isTRUE(all.equal(.trough, 0.4)))
  warning("LMR trajectory minimum is not at Bin 0.4; check data/sample/astro_bin_means.csv.",
          call. = FALSE)

p_b <- ggplot(lmr_df, aes(x=bin, y=LMR)) +
  geom_hline(yintercept=1, linetype="dotted", color="grey50") +
  geom_smooth(method="lm", se=FALSE, linetype="dashed", color="#2980b9", linewidth=0.9) +
  geom_line(color="#2980b9", linewidth=1) + geom_point(color="#2980b9", size=2.5) +
  geom_vline(xintercept=0.4, linetype="dashed", color="#2980b9", alpha=0.8) +
  geom_vline(xintercept=0.5, linetype="dashed", color="#16a085", alpha=0.8) +
  annotate("text", x=0.4, y=max(lmr_df$LMR), label="LMR\n(Bin 0.4)",
           color="#2980b9", size=2.8, vjust=0) +
  annotate("text", x=0.5, y=max(lmr_df$LMR), label="Iron\n(Bin 0.5)",
           color="#16a085", size=2.8, vjust=0) +
  scale_x_continuous(breaks=seq(0.2,0.9,0.1)) +
  labs(x="Pseudo-progression Bin", y="Normalized (Bin 0.2 = 1.0)", title="B") +
  theme_paper

# -- C: Metabolic overload -----------------------------------------------------
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
  geom_hline(yintercept=1, linetype="dotted", color="grey50") +
  geom_line(linewidth=1) + geom_point(size=2.5) +
  scale_color_manual(values=ol_colors) +
  scale_x_continuous(breaks=seq(0.2,0.9,0.1)) +
  labs(x="Pseudo-progression Bin", y="Normalized (Bin 0.2 = 1.0)", title="C", color=NULL) +
  theme_paper + theme(legend.position="right", legend.text=element_text(size=12))

# -- Assemble: A on top, B | C on the bottom -----------------------------------
supp_fig2 <- p_a / (p_b | p_c) + patchwork::plot_layout(heights=c(1, 1))

ggsave(file.path(FIG_OUT, "Supplemental_Figure2.png"), supp_fig2,
       width=13, height=10, dpi=300, bg="white")
ggsave(file.path(FIG_OUT, "Supplemental_Figure2.tif"), supp_fig2,
       width=13, height=10, dpi=600, bg="white", device="tiff", compression="lzw")
message("Saved: ", file.path(FIG_OUT, "Supplemental_Figure2.{png,tif}"))
