# =============================================================================
# Fig2_Assemble.R -- Figure 2 as it appears in the manuscript
#
# Figure 2 is four panels: A-C from Fig2_CrossCellular.R and D from
# Fig2D_Specificity_Control.R, stacked. They are produced by separate scripts
# because D uses a different input (the housekeeping control table), so this
# script joins them and aligns the D tag to the A tag column.
#
# Paths are RELATIVE to the repository root.
#   Input : output/figures/Fig2_ABC.png, output/figures/Figure2D_specificity_control.png
#   Output: figures/Figure2.{png,tif}
# =============================================================================
suppressMessages({ library(png); library(grid) })

abc_f <- "output/figures/Fig2_ABC.png"
d_f   <- "output/figures/Figure2D_specificity_control.png"
for (f in c(abc_f, d_f))
  if (!file.exists(f)) stop("missing ", f, " -- run Fig2_CrossCellular.R and Fig2D_Specificity_Control.R first")

abc <- readPNG(abc_f); d <- readPNG(d_f)
stopifnot(dim(abc)[2] == dim(d)[2])          # same width, or the tags cannot align

## x position of the leftmost dark pixel in the top-left corner = the panel tag
tag_x <- function(a, band = 0.12, win = 0.12) {
  g <- if (length(dim(a)) == 3) apply(a[, , 1:3, drop = FALSE], c(1, 2), mean) else a
  r <- g[seq_len(round(nrow(g) * band)), seq_len(round(ncol(g) * win)), drop = FALSE]
  cols <- which(apply(r < 0.47, 2, any))
  if (length(cols)) min(cols) else NA_integer_
}
dx  <- tag_x(abc) - tag_x(d)                  # shift D so its tag sits under A's
gap <- 40
H <- nrow(abc) + gap + nrow(d); W <- ncol(abc)
out <- array(1, c(H, W, 3))
out[seq_len(nrow(abc)), , ] <- abc[, , 1:3, drop = FALSE]
xs  <- seq_len(ncol(d)) + dx
keep <- xs >= 1 & xs <= W
out[nrow(abc) + gap + seq_len(nrow(d)), xs[keep], ] <- d[, keep, 1:3, drop = FALSE]

dir.create("figures", showWarnings = FALSE, recursive = TRUE)
png::writePNG(out, "figures/Figure2.png")
message(sprintf("Saved figures/Figure2.png (%d x %d); panel D shifted %+d px", W, H, dx))
