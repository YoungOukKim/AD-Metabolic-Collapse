# =============================================================================
# Fig2D_Specificity_Control.R  ->  Figure 2, panel D
#
# The cross-cellular partial correlations in panels A-C are real associations,
# but they are not specific to MCT4: a 14-gene housekeeping panel run through the
# identical procedure, in both directions, produces correlations of comparable
# magnitude. This panel shows that directly, so the negative control is visible
# rather than merely asserted.
#
# NOTE (v15): the facet titles previously carried the internal model codes
# "[A]" and "[C1]". Those collided with the figure's own panel letters A and C,
# so a reader could take "[A]" inside panel D for panel A. The codes are gone;
# the conditions are now named in full. The underlying column d$control is
# unchanged, so hk_control_by_gene.csv does not change.
#
# Signal pairs are the three MCT4 couplings. The housekeeping set is 28
# comparisons: astrocytic MCT4 against each neuronal housekeeping gene, and each
# astrocytic housekeeping gene against the neuronal V-ATPase composite (the six
# subunits used throughout this study).
#
# Paths are RELATIVE to the repository root.
#   Input : output/tables/hk_control_by_gene.csv   (written by R/audit/P2_run_all.R)
#   Output: figures/Figure2D_specificity_control.{png,tif}
# =============================================================================
suppressMessages({ library(ggplot2); library(patchwork) })
source("R/figures/utils.R")

IN  <- Sys.getenv("HK_TABLE",  unset = "hk_control_by_gene.csv")
OUT <- Sys.getenv("FIG_OUT",   unset = "figures")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

d <- read.csv(IN, stringsAsFactors = FALSE)
d$control <- factor(d$control, levels = c("A", "C1"),
                    labels = c("CPS-adjusted", "CPS + genome-wide control"))
d$absr <- abs(d$r)

sig <- d[d$set == "signal", ]
hk  <- d[d$set == "housekeeping", ]
gap <- do.call(rbind, lapply(split(d, d$control), function(s)
  data.frame(control = s$control[1],
             gap  = mean(abs(s$r[s$set == "signal"])) - mean(abs(s$r[s$set == "housekeeping"])),
             msig = mean(abs(s$r[s$set == "signal"])),
             mhk  = mean(abs(s$r[s$set == "housekeeping"])))))
gap$lab <- sprintf("specificity gap  %.3f", gap$gap)
top <- do.call(rbind, lapply(split(hk, hk$control), function(s) s[which.max(s$absr), ]))

set.seed(42)
p <- ggplot() +
  geom_hline(data = gap, aes(yintercept = mhk),  linetype = "dashed", colour = "grey55", linewidth = 0.4) +
  geom_hline(data = gap, aes(yintercept = msig), linetype = "dashed", colour = "#922B21", linewidth = 0.4) +
  geom_jitter(data = hk, aes(x = "Housekeeping\ncontrol (n = 28)", y = absr),
              width = 0.16, height = 0, size = 2.1, alpha = 0.65, colour = "grey35") +
  geom_point(data = sig, aes(x = "MCT4 signal\npairs (n = 3)", y = absr),
             size = 3.4, colour = "#922B21") +
  geom_text(data = sig, aes(x = "MCT4 signal\npairs (n = 3)", y = absr,
                            label = sub("MCT4 - Neuron ", "", label)),
            hjust = -0.22, size = 4.2, colour = "#922B21") +
  geom_point(data = top, aes(x = "Housekeeping\ncontrol (n = 28)", y = absr),
             shape = 21, size = 3.4, fill = "white", colour = "grey15", stroke = 0.9) +
  geom_text(data = top, aes(x = "Housekeeping\ncontrol (n = 28)", y = absr,
                            label = sprintf("%s  %.3f", gene, absr)),
            hjust = -0.2, size = 4.2, colour = "grey15") +
  geom_text(data = gap, aes(x = 1.5, y = 0.72, label = lab), size = 4.5,
            fontface = "italic", colour = "grey25") +
  facet_wrap(~ control) +
  scale_y_continuous(limits = c(0, 0.78), expand = c(0, 0)) +
  labs(x = NULL, y = expression("|partial " * italic(r) * "|")) +
  theme_bw(base_size = 15) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.background = element_rect(fill = "grey95", colour = NA),
        strip.text = element_text(face = "bold", size = 14),
        axis.text.x = element_text(size = 13))

p <- p + panel_tag("D")
ggsave(file.path(OUT, "Figure2D_specificity_control.png"), p,
       width = 15, height = 5.2, dpi = 300, bg = "white")
ggsave(file.path(OUT, "Figure2D_specificity_control.tif"), p,
       width = 15, height = 5.2, dpi = 600, bg = "white", device = "tiff", compression = "lzw")
cat(sprintf("Saved: %s/Figure2D_specificity_control.{png,tif}\n", OUT))
print(gap[, c("control", "msig", "mhk", "gap")], row.names = FALSE)
