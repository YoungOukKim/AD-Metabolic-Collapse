# AD-Metabolic-Collapse

**Revision round, npj Dementia submission 9b600979.** This repository is at the
state of the revised manuscript, not the original submission. Two script names
were changed in this round because their previous names carried claims the
revision withdrew:

| was | is now | why |
| --- | --- | --- |
| `Fig5_CSF_Validation.R` | `Fig5_CSF_Orthogonal_Test.R` | the CSF analysis is an orthogonal test, not validation of the mechanism |
| `SuppFig2_Temporal.R` | `SuppFig2_Axis_Position.R` | the axis is cross-sectional; every within-individual timing claim was withdrawn |

`CHANGELOG.md` records what changed in the deposited files and why, newest first.


**Coordinated astrocytic metabolic decline and the energy-starved lysosome hypothesis in Alzheimer's disease**

Code and derived data for the manuscript. This is the revision build: it supersedes
the initial deposit, and the changes are listed under *What changed* below.

---

## What changed in this revision

Three corrections and three additions, all author-initiated.

**Corrections**

1. **Extraction convention.** `R/data_extraction/01_extract_seaad.R` applied `log1p`
   to a matrix that SEA-AD already stores normalised, i.e. it logged an already
   logged matrix. `R/data_extraction/02_global_expression_sensitivity.R` did not.
   The manuscript therefore quoted CPS-adjusted partial correlations from one
   convention and genome-wide-controlled ones from the other. The extraction now
   tests the matrix before logging, and every reported value uses the single-log
   convention.

2. **Neuronal V-ATPase composite.** This study uses ten subunits for astrocytes and
   six for neurons. An earlier audit script used ten for both. All couplings are now
   computed on the six-subunit neuronal composite (`ATP6V1A`, `ATP6V1B2`, `ATP6V0A1`,
   `ATP6V0C`, `ATP6V0D1`, `ATP6V1E1`); the ten-subunit version is reported as a
   sensitivity analysis and changes the estimate immaterially (+0.226 vs +0.232).

3. **MCT4 temporal onset.** The manuscript stated that astrocytic MCT4 declines from
   Bin 0.3 and precedes the iron-gene co-decline. The deposited bin means show the
   opposite: SLC16A3 is 14.3% *above* its Bin 0.2 value at Bin 0.3, peaks at Bin 0.4
   (+31.1%), and first falls at Bin 0.5, concurrent with FTH1 and FTL, with TFRC
   earlier at Bin 0.3. The claim of temporal precedence is withdrawn. The event table
   in `R/figures/SuppFig2_Axis_Position.R` was hard-coded; it is now computed from the
   deposited CSV, with a `stopifnot` guard against silent drift.

**Additions**

4. **Housekeeping negative control.** Fourteen housekeeping genes run through the
   identical procedure in both directions produce partial correlations of comparable
   magnitude to the signal pairs (mean |r| 0.331 vs 0.464; the largest, astrocytic
   PPIA against neuronal V-ATPase, reaches +0.651). The specificity gap is 0.133,
   below a 0.15 threshold fixed before the analysis. Donor-level cross-cellular
   partial correlation is therefore **not** established as specific to MCT4, and the
   manuscript says so. Figure 2 gains a panel D showing the signal pairs against the
   housekeeping distribution, so the control is visible rather than asserted; the
   28 individual comparisons are in `output/tables/hk_control_by_gene.csv`.

5. **Detection-matched null and ambient control.** Against 2,042 genes matched on
   detection rate, the astrocytic MCT4 decline sits at the **0.39th percentile**.
   Microglia, the highest-expressing population for SLC16A3, decline three-fold less
   than astrocytes (-14.4% vs -43.2%), which excludes ambient RNA.

6. **Mediation, matched null and cross-cohort replication.** See `R/mediation/`.

7. **External replication in a second cohort.** `R/replication/70_ITG_replication.R` tests the
   astrocytic MCT4 decline in the inferior temporal gyrus of an independent astrocyte and
   microglia atlas (Serrano-Pozo et al., Nat Neurosci 2024; 31 donors, four pathology stages).
   The script applies that study's own quality control and gates on reproducing its published
   cohort size before any result is read. Outputs are in `output/tables/ITG_replication/`.

   Result: the decline is directionally consistent but the cohort is underpowered for an effect
   of this size (rho = -0.26, 95% CI [-0.57, +0.11], n = 31; 80% power requires |rho| >= 0.48).
   The detection-matched specificity does not reproduce. The ambient control does: microglia
   express SLC16A3 seven-fold more highly yet do not change. Neurons are absent from that
   dataset, so the cross-cellular coupling was not tested there.

   The raw matrices are not redistributable and are not included; obtain them from
   <https://ad-progression-atlas.partners.org> and set `ASTMET_DIR`.

**Renumbering.** A full pass over the manuscript found that the cross-cohort figure was
cited before the iron-pathway figure, and that the mediation table was cited before the
CSF table. Figures 6 and 7 and Tables 2 and 3 were therefore swapped so that every
figure, table, supplementary item and reference is now first cited in ascending order.
In this repository the cross-cohort figure is `figures/Figure6.png`, produced by
`R/mediation/Fig6_CrossCohort_Mediation.R`. Content is unchanged.

---

## Repository structure

```
AD-Metabolic-Collapse/
├── README.md
├── R/
│   ├── data_extraction/
│   │   └── 01_extract_seaad.R          per-cell extraction (storage-convention guard)
│   ├── audit/
│   │   └── P2_run_all.R                single pass; every control in the manuscript
│   ├── figures/
│   │   ├── SuppFig2_Axis_Position.R         Supplemental Figure 2 (onset computed, not hard-coded)
│   ├── Fig2D_Specificity_Control.R Figure 2 panel D: housekeeping specificity control
│   ├── Fig2_Redraw.py             Figure 2 as printed (matplotlib; asserts published values)
│   ├── Fig6_Redraw.py             Figure 6 as printed (matplotlib; asserts published values)
│   └── (see R/replication/ for the external cohort)
│   └── mediation/
│       ├── 46_composition_vs_percell.R shift-share: composition vs per-cell
│       ├── 47_braak_common_axis.R      cross-cohort on a common Braak axis
│       ├── 49_mediation_P2_definitions.R  mediation on this study's own composites
│       ├── 50_mediation_matched_null.R    detection- and effect-size-matched null
│       └── Fig6_CrossCohort_Mediation.R   Figure 6
├── R/replication/
│   └── 70_ITG_replication.R            external replication, inferior temporal gyrus
├── data/sample/                        bin- and donor-level summaries (single-log)
├── output/tables/                      every table the audit and mediation produce
└── figures/                            Figure 2 panel D, Figure 6, Supplemental Figure 2
```

---

## Panel labels

Every figure labels its panels with an uppercase letter set through `panel_tag()`
in `R/figures/utils.R`. Three different mechanisms were in use before this
revision -- lowercase titles, uppercase titles and patchwork `tag_levels` -- which
placed the letter at different offsets, so panels tagged one way did not line up
with panels tagged another, and the lowercase scripts did not reproduce the
uppercase letters printed in the manuscript. All scripts now share one mechanism
and one size (`PANEL_TAG_SIZE`), so running them reproduces the published figures.

Figure 2 is assembled from two scripts: `Fig2_CrossCellular.R` (panels A-C) and
`Fig2D_Specificity_Control.R` (panel D), joined by `Fig2_Assemble.R`.

### Figure 2, v15

Three rendering defects were found in the submitted Figure 2 and are fixed here.

1. The panel D facet titles carried the internal model codes `[A]` and `[C1]`,
   which collided with the panel letters A and C of the figure itself.
2. The panel C x-axis label was clipped at the canvas edge.
3. In panel D the `V-ATPase` label was struck through by the signal-mean line.
4. The panel B tag overlapped the legend swatch, and the D tag did not begin at
   the same x as the A tag.

The published figure is now rendered by `R/figures/Fig2_Redraw.py`, which reads
the same two deposited inputs the R scripts use -- `data/sample/donor_level_summary.csv`
and `output/tables/hk_control_by_gene.csv` -- and asserts every published value
before it draws anything: the zero-order r of +0.517, the partial r of +0.474 at
p = 5.86e-06 on n - k - 2 degrees of freedom, the five panel B partials, and the
specificity gaps of 0.133 and 0.052. If any assertion fails the script stops and
no figure is written.

The R scripts are retained and their labels corrected, so they remain the record
of how the panels were originally specified. They are not byte-identical to the
published rendering: the horizontal jitter of the housekeeping points in panel D
is drawn from a different generator, and font metrics differ between ggplot2 and
matplotlib. Every plotted value is the same. `Fig2_Redraw.py` is the script of
record for the figure as printed.

## Reproducing

R >= 4.3.0.

```r
install.packages(c("data.table","ggplot2","dplyr","tidyr","patchwork","scales"))
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install("rhdf5")
```

Paths are **relative to the repository root**. Set the working directory there:

```r
setwd("path/to/AD-Metabolic-Collapse")
```

or override with environment variables `SEAAD_H5AD`, `ROSMAP_ASTRO`, `ROSMAP_CLIN`,
`P2_OUT_DIR`.

```r
source("R/audit/P2_run_all.R")        # gates, couplings, controls, sample CSVs
source("R/figures/SuppFig2_Axis_Position.R")  # needs only data/sample/
```

`P2_run_all.R` accepts `selftest` as an argument and then runs its logic checks
without touching any data. Blocks are sized by non-zero entry count
(`P2_TARGET_NNZ`, default `2e7`); lower it if memory-constrained. The full pass over
1,378,211 nuclei takes roughly 80 minutes.

---

## Audit

`R/audit/P2_run_all.R` re-derives, in one pass over the raw h5ad, every control the
manuscript reports: the reproduction gate, the anchor gate, the donor-level partial
correlations under three covariate specifications and at four CPS thresholds, the
housekeeping negative control, leave-one-donor-out and bootstrap stress tests, the
cell-type comparison, the ambient control, the detection-matched null and the binned
trajectory. Every threshold is fixed in the script before any number is computed, and
the script checks its own output against the values printed in the manuscript.

Two gates must pass before anything else is interpretable, and both are printed first.

* **Reproduction gate.** Six previously published bin-level percentage changes
  (SLC16A3 -43.2, HK2 -35.2, LDHA -20.8, PDK1 -16.5, SLC16A1 -11.1, SLC2A1 -7.4) are
  recomputed from the raw object and must agree to within one percentage point. They
  agree exactly.
* **Anchor gate.** Astrocytic SLC16A3 = -0.5020 and PTGDS = -0.2485, global-corrected
  partial Spearman at the donor level.

Outputs are written to `output/tables/` and are included in this repository.

---

## Data requirements

**SEA-AD** — not included; requires a data access agreement.
`SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad`, 1,378,211 nuclei, 84 donors.
[portal.brain-map.org](https://portal.brain-map.org)

**ROSMAP DLPFC** — not included; used only by `R/mediation/46` and `47`.

**ADNI Emory TMT-MS CSF proteomics** — not included; n = 1,105.
[adni.loni.usc.edu](https://adni.loni.usc.edu)

`data/sample/` holds bin- and donor-level summaries derived from SEA-AD and is
sufficient to reproduce the figures and tables that do not require raw access.

---

## License

Code: MIT. Data: subject to the SEA-AD, ROSMAP and ADNI data use agreements and not
redistributable.

## Severely-affected donor sensitivity (revision addition)

Eleven donors are flagged `Y` in `obs/Severely Affected Donor` in the released
h5ad. They are retained in the primary analysis; `R/sensitivity/` reproduces the
donor-level results directly from the h5ad and reports what changes when they
are removed.

```bash
Rscript R/sensitivity/01_extract_h5ad.R selftest   # logic only, no file access
Rscript R/sensitivity/01_extract_h5ad.R            # ~45 min, single pass
Rscript R/sensitivity/02_sa_sensitivity.R          # seconds, reads obs only
Rscript R/sensitivity/03_specificity_screen.R BSG  # seconds
```

Step 1 writes `output/sensitivity/donor_level.csv` and `donor_by_gene.rds`;
steps 2 and 3 read those and never re-read `X`.

Two gates run before any result is used. **G0** checks the nucleus and
cell-type counts against the Methods (1,378,211 / 67,419 / 671,689 / 84).
**G1** checks the cached table against the published donor-level values
(slope -0.0737, partial r +0.474 / +0.501 / +0.418) and refuses to print a
sensitivity result if it fails.

Every script self-tests before it opens a file and exits non-zero on failure.
Stochastic quantities are never judged on a single draw: the partial-correlation
implementation is checked against theory over 2,000 replicates, and the
mediation estimator against a median over 200 replicates.

### Two things that are easy to get wrong

**AnnData categorical NA codes.** A missing category is stored as `-1`. In R,
`categories[codes + 1]` silently drops index 0, which shortens the vector and
misaligns it against every other `obs` column. In this cohort that affects
`Severely Affected Donor`, where neurotypical reference nuclei carry no value.
`obs_cat()` preserves the NA and asserts the returned length.

**Cohort filter before counting.** Neurotypical reference donors carry no CPS.
Counting cell types without excluding them overstates astrocytes by 2,590 and
excitatory neurons by 84,687; the asymmetry reflects NeuN sorting in the
reference data. G0 catches this.

## The detection-matched null (settled)

`R/audit/04_detection_matched_null.R` was run against the raw h5ad with every
filter stated. It computes two astrocytic detection rates for SLC16A3, because
the manuscript uses both:

| quantity | value | used for |
| --- | --- | --- |
| detection, all astrocytes in the cohort | 0.0612 (6.1%) | the microglia comparison in the Limitations |
| detection, astrocytes in Bins 0.2-0.9 | 0.0602 (6.0%) | the null, matched to the early-to-late change |
| early-to-late change | -43.2% | |
| matched band, +/-25% | 0.0452 - 0.0753 | |
| genes in band | 2,026 | target excluded from its own null |
| SLC16A3 percentile | 0.39 | |

Restricting to Bins 0.2-0.9 moves the detection rate by 0.001, not by the 0.006
that would have been needed to explain the earlier 6.1% versus 6.7% split. The
0.0669 and 0.066 values in the superseded files came from a different pipeline
whose filters are not recorded; across 16,160 shared genes that pipeline agrees
with this one to a median ratio of 1.003, but SLC16A3 itself sits in the top 5%
of the disagreement. The authoritative values are the ones above.

Propagated to: the Figure 6 caption, the Results text, `data/38_all_genes.csv`
(regenerated from `output/tables/detection_matched_null_all_genes.csv`),
`output/tables/SuppTable_detection_null.csv` and Figure 6C, which now renders
2,026 genes at the 0.39th percentile.

## Paths

Figure and mediation scripts read `tables/` and `data/` relative to the
repository root. Those directories now hold the files those scripts expect;
`output/tables/` remains the canonical location written by the analysis scripts.


## Figure 6

`R/figures/Fig6_CrossCohort_Mediation.py` renders Figure 6 with matplotlib and
has no R dependency:

```bash
python3 R/figures/Fig6_CrossCohort_Mediation.py
```

It reads `tables/Table4_ROSMAP.csv`, `tables/Table3_mediation.csv`,
`data/38_all_genes.csv` and `data/47_chain_braak.csv`, and writes
`figures/Figure6.png` and `.tif` at 300 dpi. The earlier renderer left its own
output under `output/figures/`; those files were removed because they were an
earlier version of the figure and did not match the printed one.

Panel B is driven entirely by `tables/Table3_mediation.csv`; adding or removing
a row there changes the bars without touching the figure code. The reverse panel
now carries six bars, the sixth being neuronal LAMP1 at 21.4%.

Panel C recomputes its band at render time from `data/38_all_genes.csv`
(detection rate +/- 25%) and now agrees with the caption: 2,026 genes at the
0.39th percentile.


## Supplemental Table 1(b): normalisation

`R/audit/78_supp_table_1b.R` recomputes the Bin 0.1 sensitivity analysis. Every
published value in the primary column reproduces exactly from the h5ad-derived
cache: MCT4 slope -1.036 (p = 0.011), astrocytic V-ATPase slope -0.044
(p = 0.244), delta-slope p = 0.0005, MCT4 -43.2%, V-ATPase -0.8%. The astrocytic
V-ATPase composite is the ten subunits defined in
`R/data_extraction/01_extract_seaad.R`, not the six used for neurons.

The two columns of the published table did not share a normalisation. Methods
state Bin 0.2 = 1.0; the sensitivity column normalised the Bin 0.1-included
trajectory to Bin 0.1 instead, which reproduces its published -0.727 exactly.
On the common scale the value is -2.152: including the high-leverage early bin
steepens the MCT4 slope rather than attenuating it, and the delta-slope contrast
strengthens from p = 0.0217 to p = 0.0048. Both columns now use Bin 0.2 = 1.0.

`output/tables/SuppTable_1b_recomputed.csv` carries all three conventions side by
side so the change is auditable.


## Supplemental Table 1(a): the ANLS row

Three of the four rows were recomputed directly from the h5ad and match the
manuscript exactly, including the attenuation percentages (-8.2%, -8.8%, -2.2%).
The fourth, ANLS against neuronal V-ATPase, needs the ANLS composite, which the
donor-level export does not carry. `R/audit/79_anls_row.R` rebuilds it from the
donor x gene matrix and verifies the row, gated on the three MCT4 rows
reproducing first.

The manuscript's embedded supplementary table already carried the current values
for all four rows. The separate `.xlsx` was a stale copy of an earlier version;
all four rows in it have now been brought into line.

ANLS composite = mean(SLC2A1, LDHA, SLC16A1). MCT4/SLC16A3 is analysed separately
and is deliberately not part of the composite, because its trajectory (-43.2%)
would dominate a mean whose other members move -7.4% to -20.8%.


## Degrees of freedom for partial correlations

Every partial correlation p-value in the manuscript was computed on n - 2
degrees of freedom. For a partial correlation controlling k variables the
correct value is n - k - 2; n - 2 is right only for the zero-order column.
The published p-values were therefore slightly anti-conservative.

`R/audit/80_partial_df_audit.R` enumerates all seventeen specifications that
appear in the paper, recomputes each from the donor-level data under both
conventions, and gates on the correlation coefficients reproducing first. This
matters because the specifications are described in prose and are not
interchangeable: three rows residualise each variable on its own cell type's
global mean, which is a different procedure from entering both means as shared
covariates and gives a different coefficient (+0.281 against +0.232).

All seventeen reproduce on r (`output/tables/partial_df_audit.csv`), including
the two ANLS rows, which need the donor x gene matrix.

No p-value crosses 0.05. Ratios range from 1.05 to 1.18. The only wording that
changed as a consequence is "all p < 0.025" at the CPS thresholds, which becomes
"all p < 0.03". The Methods now state the convention explicitly.


## Repository layout

```
config.R                       not used; see R/sensitivity/config_sensitivity.R
data/                          figure inputs; place the h5ad here (git-ignored)
output/tables/                 canonical tables written by the analysis scripts
output/figures/                Figure 6
output/sensitivity/            written by R/sensitivity and R/audit
R/data_extraction/             original SEA-AD extraction
R/mediation/                   mediation and cross-cohort analyses
R/replication/                 ITG replication
R/figures/                     figure scripts (Fig6 is Python/matplotlib)
R/sensitivity/                 severely-affected donor sensitivity, specificity screen
R/audit/                       post-hoc audits: detection null, Supp Table 1(a) and 1(b),
                               partial-correlation degrees of freedom
```

Every script is run from the repository root. Paths are relative throughout; no
absolute path appears anywhere in the code. All comments are in English.

Each analysis script self-tests before it opens a file and exits non-zero on
failure:

```bash
for f in R/sensitivity/0[123]_*.R R/audit/*.R; do Rscript "$f" selftest; done
python3 R/figures/Fig6_CrossCohort_Mediation.py
```

## Reproducing the revision

```bash
Rscript R/sensitivity/01_extract_h5ad.R        # one pass over the h5ad, ~45 min
Rscript R/sensitivity/02_sa_sensitivity.R      # severely-affected donor sensitivity
Rscript R/sensitivity/03_specificity_screen.R  # expression-matched screen
Rscript R/audit/04_detection_matched_null.R    # detection-matched null, authoritative
Rscript R/audit/78_supp_table_1b.R             # Bin 0.1 sensitivity, one normalisation
Rscript R/audit/79_anls_row.R                  # ANLS row of Supplemental Table 1(a)
Rscript R/audit/80_partial_df_audit.R          # partial-correlation degrees of freedom
python3  R/figures/Fig6_CrossCohort_Mediation.py
```

Step 1 is the only expensive step; everything after it reads its cached output.


## APOE as a candidate confound of the iron axis

`R/audit/81_apoe_confound.R` tests whether the reported findings are an APOE
composition effect. The question is prompted by Ayton et al. (Nat Commun 6:6760,
2015), who show that APOE genotype sets CSF ferritin: epsilon-4 carriers run 22%
higher.

The composition path is open in this cohort. The epsilon-4 carrier fraction
rises along the progression axis, from 11% in the lowest CPS tertile to 55% in
the highest (logistic beta = +5.63, p = 6.2e-4, 25 of 84 donors are carriers).

Thresholds were fixed before the data were read: a marker is called confounded
if adjustment for epsilon-4 moves its CPS slope by more than 20%, or a partial
correlation by more than 0.05.

| marker | unadjusted | epsilon-4 adjusted | shift | verdict |
| --- | --- | --- | --- | --- |
| MCT4 | -0.0737 | -0.0823 | 11.6% | robust |
| MCT4-V-ATPase (partial) | +0.474 | +0.458 | 0.017 | robust |
| MCT4-LAMP1 (partial) | +0.501 | +0.488 | 0.013 | robust |
| MCT4-LDHB (partial) | +0.418 | +0.400 | 0.018 | robust |
| FTH1 | -1.015 | -1.190 | 17.2% | robust |
| FTL | -0.405 | -0.482 | 18.8% | robust |
| SLC11A2 | -0.096 | -0.110 | 14.7% | robust |
| TFRC | -0.136 | -0.164 | 20.8% | confounded |
| HMOX1 | -0.0152 | -0.0115 | 24.4% | confounded |
| SLC40A1 | -0.0046 | -0.0078 | 69.2% | confounded |

The MCT4 axis survives, and adjustment steepens rather than attenuates it: APOE
was partially masking the decline. Three iron markers cross the threshold. TFRC
clears it by 0.8 percentage points and also steepens, so the decline itself
stands. HMOX1 is the only marker that weakens under adjustment, and it is the
same gene the P3 audit rejected as a cross-region anchor. SLC40A1 has an
unadjusted slope of -0.0046, so its 69.2% shift is an unstable ratio on a
near-zero denominator; the verdict stands under the pre-set rule but the row is
not used as evidence.

The self-test checks the detector itself: a confound is planted in one synthetic
marker and absent from another, and the rule must fire on the first and stay
silent on the second. Without that check, an output of "no confound" could not
be distinguished from a broken detector.

Nothing here changes a published value. It is recorded because the composition
shift is a real cohort property that a reviewer can find.

---

## Revision round (npj Dementia, submission 9b600979)

Everything added for the revision lives under `R/revision2/` and `output/revision2/`.
Both scripts run from the repository root and take every path relative to it; set the
environment variables listed in each header if your data sit elsewhere.

```
Rscript R/revision2/90_revision2_all.R      # S1-S9: reproduction gates, acid-base panel,
                                            # mitochondrial panel, pathway screen, staging axes,
                                            # leverage composition, curve features
Rscript R/revision2/91_revision2_S12_all.R  # S12: four-compartment decomposition, mediation
                                            # competition, panel-definition reconciliation,
                                            # sparse-cell guard. One pass over the h5ad.
```

### What the revision added

| Block | Question it answers | Output |
|---|---|---|
| S12-P | Does the detection-matched null separate lactate export from iron handling? | `S12_panels_from_S3.csv`, `S12_astro_ranked_by_percentile.csv` |
| S12-A / S12-D | Is the SLC16A3 decline restricted to astrocytes? | `S12_A_SLC16A3_by_compartment.csv`, `S12_ambient_symmetry.csv`, `S12_D_microglial_glycolysis.csv` |
| S12-N | Does the lactate composite decline under the definition given in Methods? | `S12_N_pathway_screen_manuscript_definition.csv` |
| S12-M | Does the mediation analysis separate MCT4 from competing astrocytic mediators? | `S12_M_mediation_competition.csv`, `S12_M_null_membership.csv`, `S12_M_reproduction_gate.csv` |
| S12-H | Do the leverage-composition chi-square tests survive a sparse-cell guard? | `S7_leverage_composition_guarded.csv` |

Verdict rules were written to `output/revision2/S12_RULES.txt` before any data file was
opened. Three of them fired against the original claims; the verdicts as recorded at run
time are in `output/revision2/S12_VERDICTS.txt`.

### Inputs

Third-party data are not redistributed. See `data/README.md` for the four inputs the
revision scripts expect and the environment variables that override their locations.
`output/sensitivity/README.md` explains the one intermediate that is regenerated rather
than committed (`donor_by_gene.rds`). Note that `91_revision2_S12_all.R` does not need it:
it rebuilds the donor x gene matrices in its own h5ad pass and cross-checks them against
that file when it is present.

| Block | Needs the h5ad? | Needs `donor_by_gene.rds`? |
|---|---|---|
| S12-P (panels from S3) | no | no |
| S12-H (sparse-cell guard) | no | no |
| S12-A / S12-D (four compartments) | yes | no |
| S12-N (pathway screen) | yes | no |
| S12-M (mediation competition) | yes | no |
| `90_revision2_all.R` S6 pathway screen and C1 gate | yes | yes |

S12-P and S12-H therefore run in seconds from the deposited
`output/revision2/S3_detection_by_celltype.csv` and
`output/revision2/S7_leverage_composition.csv` alone, with no raw data at all.

### What is deliberately not in this repository

Nothing derived from a restricted dataset at individual or per-donor resolution is
committed here, and nothing that a user could not regenerate from the source with the code
provided.

| Not included | Why |
|---|---|
| The SEA-AD h5ad, the Emory CSF matrix, the Elecsys file, ADNIMERGE2, ROSMAP | Third-party datasets under their own terms; see `data/README.md` |
| `output/sensitivity/donor_by_gene.rds` | Donor x gene matrices computed from the h5ad. Committing it would redistribute SEA-AD in derived form, and anyone running the code needs the h5ad anyway |
| `output/tables/rosmap_donor_mct4.csv` | Per-donor astrocytic MCT4 for 430 ROSMAP donors, keyed by ROSMAP identifier. ROSMAP is released under a RADC data use agreement, so per-donor derived values are not redistributable. `R/replication/82_rosmap_braak_matched_donors.R` regenerates it from the source in about a minute, and the three-row result it feeds is committed as `output/tables/82_rosmap_braak_matched.csv` |
| Per-subject rows of the `Fig5_data` sheet | ADNI individual-level data keyed by RID. The repository copy of the source-data workbook is `Supporting_Data_Values_v7.xlsx`, with that sheet replaced by a pointer; every aggregate sheet is intact |

`output/revision2/S1_donor_stage_full.csv` carries the SEA-AD donor identifiers together
with their neuropathological and demographic annotations. Every field in it is read
directly from the `obs` table of the SEA-AD h5ad, which is an open download requiring no
application, so this is a redistribution of already-public data rather than a new
disclosure. It is kept because the staging-axis analyses (S8) can then be checked without
downloading the atlas, and because two later blocks read the integer staging columns back
from it.

One field is treated differently. `Age.at.Death` is released by SEA-AD as an exact value,
including 46 donors above 89 and a maximum of 102. Ages over 89 are conventionally
aggregated, so in this copy that column is capped at `90+`. No script reads it — the
scripts take age from `obs` directly when they need it — so the cap changes nothing that
is computed here.

### Conventions fixed in this round

Empirical p values from any matched-null test are `(k + 1) / (n + 1)`, where `k` is the number
of matched control genes reaching or exceeding the observed value and `n` is the size of that
gene's matched null. The `+1` correction is used because a permutation p of exactly zero is not
attainable, and `output/revision2/S12_M_mediation_competition.csv` carries the uncorrected
`k / n` alongside it in `p_empirical_uncorrected` so the difference is visible. Under the
corrected form transferrin receptor is 0.053 rather than 0.049 and ferritin light chain is
0.030 rather than 0.027. The correction is not used to reverse a verdict. It was adopted after
the values were known, and of the eighteen mediator-outcome tests it moves exactly one across
0.05, in the direction that favours this study's own claim. The pre-registered criterion is the
position of the observed value in its own matched null: transferrin receptor's mediated fraction
of 36.2% lies above that null's 95th percentile of 35.7%, so it is recorded as exceeding, and
rule M-1 fires on two competitors rather than one.

Degrees of freedom for partial correlations are n - k - 2 throughout the manuscript, where
k is the number of covariates. Two deposited files predate that convention and are kept as the
audit trail rather than silently rewritten:

  `output/tables/couplings.csv` was written by `P2_run_all.R`, whose `pcorZ` uses n - 2. The
  n - 2 columns (`A_p`, `C1_p`, `C2_p`) are unchanged and the n - k - 2 values are added
  alongside as `A_p_nk2`, `C1_p_nk2`, `C2_p_nk2`. The manuscript quotes the n - k - 2 form.
  The small residual difference against the manuscript (6.0e-06 against 5.86e-06 for
  MCT4-V-ATPase) is display rounding: this file stores r to three decimals.

  `output/tables/partial_df_audit.csv` had a column named `p_in_manuscript`. That name was
  true when it was written and is not now, because the manuscript moved to the corrected
  values. The column holds the value as submitted and is renamed `p_in_submitted_version`.

Analyses of covariance use Type II sums of squares, as stated in Methods. Two values in an
earlier draft were Type I and have been corrected; `output/revision2/Table3_verification.csv`
records the recomputation.

### Supplementary files cited by the manuscript

| Cited as | File |
|---|---|
| Supplemental Table 1, parts a to g | `output/revision2/Supplemental_Table_1_v10.xlsx` |
| Supporting Data Values | `output/revision2/Supporting_Data_Values_v7.xlsx` (per-subject ADNI sheet removed; see above) |

Part c of Supplemental Table 1 carries the mediation competition together with the
5,000-sample bootstrap intervals on the proportion mediated and the one outcome that was
examined and dropped, the astrocytic mitochondrial electron transport chain. Those values
previously sat in a second file that the manuscript never cited; keeping the same numbers
in two places is how they drift apart, so the second copy has been removed and part c is
now the single source. `R/mediation/50_mediation_matched_null.R` still writes that file on
each run, and its contents are reproduced there.

### Two locations for the same three tables

`tables/` and `output/tables/` both hold `Table3_mediation.csv`, `Table3_reverse.csv` and
`Table4_ROSMAP.csv`, byte for byte. That is deliberate, not stray duplication:
`R/mediation/50_mediation_matched_null.R` writes to `tables/` (its `OUT_T`) and
`R/mediation/Fig6_CrossCohort_Mediation.R` reads from `tables/`, while `output/tables/`
holds the deposited copy of the same run alongside the other audit outputs. Delete either
and something breaks. The same applies to `data/sample/donor_level_summary.csv`, which is
the documented sample input and is identical to `output/tables/donor_merged_base.csv`.

### Seeds

The project seed is 42. There is one documented exception:
`R/mediation/50_mediation_matched_null.R` and the S12-M block of
`R/revision2/91_revision2_S12_all.R` use `set.seed(20260713)`, because the bootstrap
confidence intervals in Table 3 were generated with it. Changing it would change the
`boot_lo` and `boot_hi` columns only; mediated fractions, matched-null membership and
empirical p values are deterministic. All other seeds, including those inside self-test
blocks, are 42.

### Corrections carried in this round

- An earlier extraction applied `log1p` to a matrix SEA-AD already stores in normalised
  form. All values now use the single-log convention; `S2_verify_july_outputs.csv`
  reproduces the seventeen previously published summary values from the raw object.
- Partial-correlation p values now use `n - k - 2` degrees of freedom throughout;
  zero-order correlations remain `n - 2`.
- The astrocytic ANLS composite is `SLC2A1 + LDHA + SLC16A1` everywhere, the definition
  given in Methods, with MCT4 analysed separately. An intermediate revision script used a
  different set under the same name; `R/revision2/90_revision2_all.R` retains that set only
  as an explicitly labelled contrast row.
- The two columns of the Bin 0.1 sensitivity table are now normalised to a common baseline.
