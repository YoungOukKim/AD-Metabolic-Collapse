# figures/

This directory holds the three figures that were rebuilt in Python during the
revision, plus one carried over from the R pipeline. It is not the full figure
set, and the naming difference is worth stating plainly:

| file | produced by | writes to |
| --- | --- | --- |
| `Figure2.{png,tif}` | `R/figures/Fig2_Redraw.py` | this directory |
| `Figure4.{png,tif}` | `R/figures/Fig4_Redraw.py` | this directory |
| `Figure6.{png,tif}` | `R/figures/Fig6_Redraw.py` | this directory |
| `Supplemental_Figure2.{png,tif}` | `R/figures/SuppFig2_Axis_Position.R` | `output/figures/`, copied here |

The other figures in the paper -- 1, 3, 5, 7 and Supplemental Figures 1 and 3 --
are produced by the R scripts in `R/figures/`, which write to `output/figures/`.
That directory is empty in a fresh clone because the scripts have not been run.

Why the split. Three figures were redrawn in Python because rendering defects
were found in them during the revision: bracketed model codes leaking into a
panel title, a clipped axis label, a label struck through by a mean line, and in
Figure 4 a panel that printed values from the superseded extraction. Each redraw
script asserts every published value before it draws anything, so a failed
assertion stops the script rather than producing a plausible wrong figure. The R
scripts were left in place and are still the record of how the panels were
originally specified.

Which figures can be regenerated from what is deposited here:

  Figures 2, 4, 6                deposited CSVs only; run the three Python scripts
  Figures 1, 3, Supp. Figures 1, 2   deposited CSVs only; run the R scripts
  Figure 5, Figure 7             need ADNI proteomics, which requires ADNI access
  Supplemental Figure 3          a schematic; it has no data of its own. Every
                                 number printed on it is listed in the
                                 SuppFig3_data sheet of the Supporting Data
                                 Values workbook, with its source sheet.
