# CHANGES

## Distribution-robust CSF reanalysis (Figure 5, Figure 6, Table 2)

The CSF protein-level analyses were revised after we determined that several
individual-level results in the previous version were artifacts of the strongly
right-skewed TMT-MS abundance distribution rather than genuine biological
coupling. CSF TMT-MS abundances are dominated by a small number of
high-abundance samples, and untransformed Pearson correlations on these values
are driven by that leverage.

**What changed**

- **Primary statistic.** All CSF protein–protein associations are now computed
  by **Spearman rank correlation**. Untransformed Pearson is reported only as a
  labelled sensitivity comparison.
- **Partial correlations.** Computed by **rank-transforming** each variable and
  residualising the ranked response and predictor on the ranked control set by
  OLS, then correlating the residuals (`spearman_partial()` in `utils.R`).
- **Orthogonal validation.** Key Tau associations are validated against the
  **Roche Elecsys immunoassay** (total Tau) and an **independent aptamer
  platform (SomaScan)**.

**Specific claims withdrawn**

- The individual-level **CSF V1A–Tau correlation** (untransformed Pearson
  r ≈ +0.86) does **not** survive trimming, rank-based analysis
  (Spearman ρ = +0.25), confound control (ρ | GFAP = +0.05;
  ρ | GFAP+TREM2+NfL = −0.06) or validation against immunoassay Tau
  (ρ = +0.03). CSF V1A is therefore **not** interpreted as an individual-level
  Tau marker.
- The **MCI-stage "peak"** in V1A–Tau coupling (previously r = +0.92) reflected
  the same distribution leverage and is **no longer claimed**.
- The **subunit-specific axes** (V1A "neuronal-Tau" vs V1E1 "microglial-iron"),
  which were based on untransformed Pearson correlations, are **not reproduced**
  under rank-based analysis; the subunit comparison is now reported as
  exploratory with **no robust subunit-specific axis**.
- The **TFRC↔MAPT "suppression effect"** (correlation strengthening after V1A
  control, +0.38 → +0.60) is **withdrawn**; the iron–Tau claim is limited to the
  transcriptomic and group-level proteomic evidence. TFRC–Tau concordance is
  not consistent across platforms (Elecsys vs TMT-TFRC ρ = +0.08; vs
  SomaScan-TFRC ρ = −0.44).

**What is unchanged / promoted**

- **Mechanistic thesis (transcriptomic).** MCT4 (−43%) declines far faster than
  V-ATPase (−0.8%; Δslope p = 0.0005); astrocytic MCT4 couples to neuronal
  V-ATPase after pseudo-progression adjustment (donor-level partial
  r = +0.466). Unchanged.
- **Group-level CSF V-ATPase preservation.** V1A abundance does not differ
  across CN/MCI/DEM (KW p = 0.999; ANCOVA p = 0.647), consistent with structural
  pump integrity. Unchanged.
- **Glycolysis–Tau coupling (promoted).** HK1 tracks Tau under rank-based
  analysis (Spearman ρ = +0.52), after confound control (ρ | GFAP = +0.42;
  ρ | GFAP+V1A = +0.42), and reproduces against immunoassay Tau (TMT ρ = +0.18)
  and on an independent platform (SomaScan ρ = +0.31).
- **Iron pathway (group level).** TFRC declines at the group level after age/sex
  adjustment (ANCOVA p = 0.036, −4.7% DEM vs CN); LCN2 shows diagnosis-related
  elevation (ANCOVA p = 0.020). Unchanged.

**File-level changes**

- `R/figures/utils.R` — added `spearman_r()`, `pearson_r()`, `spearman_partial()`
  (rank-based partial), and `load_elecsys_tau()` / `load_somascan()` loaders.
  Legacy Pearson `partial_cor*()` helpers retained for backward compatibility
  only; they are no longer used for the primary CSF analysis.
- `R/figures/Fig5_CSF_Validation.R` — rebuilt as a distribution-robust 8-panel
  figure (robustness ladders for V1A–Tau and HK1–Tau; cross-platform validation
  vs Elecsys; subunit comparison). Panels reporting the untransformed V1A–Tau
  correlation and the MCI-stage peak were removed.
- `R/figures/Fig6_Iron_Suppression.R` → **renamed** to
  `R/figures/Fig6_Iron_Pathway.R` — rebuilt as a 2-panel figure (group-level
  TFRC decline; TFRC–Tau cross-platform inconsistency). The "suppression effect"
  panel was removed.
- `R/tables/Table2_CSF_Proteomics.R` — Part A switched to Spearman primary (raw
  Pearson in parentheses), Elecsys-anchored rows added, "no robust
  subunit-specific axis" interpretation. Part B (ANCOVA, Type I SS) retained.
- `docs/data_dictionary.md` — added Elecsys and SomaScan dataset entries.

## Manuscript context (no code change)

- **Supplemental Figure 3 — "Parallel routes to lysosomal acidification failure"**
  was added to the manuscript: a conceptual schematic placing the energetic ESL
  route (this study) alongside the structural ApoE4-cholesterol route
  (Lee et al.) converging on PANTHOS. It is a BioRender-style diagram with **no
  analysis script** and no data dependency.
- Supplemental Figures 1–2 are unchanged in content; they are assembled from the
  existing `ED_Fig1–4` panel scripts (see README mapping table).
- All transcriptomic values in the current manuscript (including the newly stated
  neuronal MCT2 / SLC16A7 +2.3% preservation) are reproducible from the existing
  `data/sample/` CSVs; no sample-data update is required.
