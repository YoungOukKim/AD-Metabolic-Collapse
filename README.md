# AD-Metabolic-Collapse

**Astrocytic lactate shuttle disruption and the energy-starved lysosome in Alzheimer's disease**

---

## Overview

This repository contains all R code and sample data used to reproduce the figures and tables in the manuscript. The analysis integrates single-nucleus transcriptomics (SEA-AD atlas) with CSF proteomics (ADNI Emory TMT-MS) to define an **'energy-starved lysosome' (ESL)** state in Alzheimer's disease, and positions it as a **complementary metabolic route** alongside the established structural (ApoE4-cholesterol) route to lysosomal acidification failure and PANTHOS (Lee et al., *Nat Neurosci* 2022).


---

## Repository Structure

```
AD-Metabolic-Collapse/
├── README.md
├── data/sample/                      # Bin-level summary CSVs (included)
├── R/
│   ├── data_extraction/              # SEA-AD extraction + sample-CSV generation
│   ├── figures/
│   │   ├── utils.R                   # Theme, helpers, ADNI/Elecsys/SomaScan loaders
│   │   ├── Fig1_Dissociation.R       Fig2_CrossCellular.R
│   │   ├── Fig3_Subtype_Network.R    Fig4_Clinical.R
│   │   ├── Fig5_CSF_Validation.R     # Fig 5 (a-h): CSF validation (distribution-robust) *
│   │   ├── Fig6_Iron_Pathway.R       # Fig 6 (a-b): iron — group decline + cross-platform *
│   │   ├── SuppFig1_Comprehensive.R  SuppFig1_Iron_Ragulator.R   # → Supplemental Figure 1
│   │   └── SuppFig2_Temporal.R       # → Supplemental Figure 2
│   └── tables/
│       ├── Table1_CrossCellular.R
│       ├── Table2_CSF_Proteomics.R          # Table 2: subunit comparison (Spearman primary) *
│       ├── Supp_Table1_GeneExpression.R     # → Supporting Data Values
│       ├── Supp_Table2_PartialCor.R         # → Supplemental Table 1 (Part A)
│       └── Supp_Table3_Sensitivity.R        # → Supplemental Table 1 (Part B)
├── analysis/
│   └── FigS/SuppFig2C_extract.R      # source-data extraction for Supplemental Figure 2C
└── docs/data_dictionary.md
```

`*` requires ADNI data access (see Data section).

### Supplemental figures (manuscript numbering)

The manuscript presents three supplemental figures. The R scripts produce the
component panels, which are assembled into the figures below:

| Supplemental Figure (manuscript) | Content | Source |
|----------------------------------|---------|--------|
| **Supplemental Figure 1** — Astrocyte transcriptomic dissociation context | comprehensive gene changes (A), TFRC–ANLS coupling (B), iron trajectories (C), Ragulator (D), pH/EAAT2/ATP1A2 (E) | `SuppFig1_Comprehensive.R`, `SuppFig1_Iron_Ragulator.R` |
| **Supplemental Figure 2** — Temporal ordering of the metabolic cascade | event timeline (A), LMR composite trajectory (B), metabolic overload / transition zone (C) | `SuppFig2_Temporal.R` (canonical; reproduces the published figure, A on top, B\|C on bottom) |
| **Supplemental Figure 3** — Parallel routes to lysosomal acidification failure | conceptual schematic: energetic (ESL, this study) **vs** structural (ApoE4-cholesterol, Lee et al.) routes converging on PANTHOS | **Author-drawn schematic; not R-generated** |

> Supplemental Figure 3 is a conceptual diagram, drawn by the authors, and has no
> analysis script. All figure scripts follow the manuscript's "Supplemental Figure"
> numbering. The cross-cellular donor-level partial correlation (CPS-adjusted) is
> shown in main-text **Figure 2 (B–C)** and is produced by `Fig2_CrossCellular.R`.

### Supplemental tables (manuscript numbering)

| Manuscript item | Source script |
|-----------------|---------------|
| **Supporting Data Values** (179-gene per-bin expression) | `Supp_Table1_GeneExpression.R` |
| **Supplemental Table 1 — Part A** (donor-level partial correlations, n = 84) | `Supp_Table2_PartialCor.R` |
| **Supplemental Table 1 — Part B** (Bin 0.1 sensitivity analysis) | `Supp_Table3_Sensitivity.R` |

---

## Data Requirements

### SEA-AD — not included, requires data access agreement
- File: `SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad`
- Access: [portal.brain-map.org](https://portal.brain-map.org); n = 1,378,211 nuclei, 84 donors

### ADNI CSF data — not included, requires ADNI access
Used by `Fig5_CSF_Validation.R`, `Fig6_Iron_Pathway.R`, `Table2_CSF_Proteomics.R`.
Access: [adni.loni.usc.edu](https://adni.loni.usc.edu)

| Dataset | Role | Used by |
|---------|------|---------|
| Emory TMT-MS CSF proteomics (n = 1,105; 3,907 proteins) | Primary CSF proteomics | Fig 5, 6, Table 2 |
| `DXSUM.rda` + `ADSL.rda` (ADNIMERGE2) | Diagnosis + demographics merge | Fig 5, 6, Table 2 |
| Roche Elecsys immunoassay (`UPENNBIOMK_ROCHE_ELECSYS_*.csv`) | Orthogonal total-Tau validation | Fig 5 (d, g), Fig 6 (b), Table 2 (Part A) |
| SomaScan 7k post-QC matrix | Independent aptamer platform | Fig 5 (g), Fig 6 (b) |

> SomaScan seqIDs are dataset-specific; `load_somascan()` defaults to the
> Cruchaga-lab ADNI SOMAscan7k matrix (HK1 = `X13131.5`, MAPT = `X5854.60`,
> TFRC = `X6895.1`). V-ATPase V1A is **not** on the SomaScan panel.

### Sample Data — included
`data/sample/` contains **bin-level summary statistics** derived from SEA-AD,
sufficient to reproduce all figures except Fig 5–6 and all tables except Table 2
without raw-data access. See `docs/data_dictionary.md`.

---

## Reproducing the Analysis

```r
# 0) dependencies
install.packages(c("ggplot2","dplyr","tidyr","patchwork","scales","data.table"))
if (!require("BiocManager")) install.packages("BiocManager"); BiocManager::install("rhdf5")
# R >= 4.3.0

# 1) working directory = repo root
setwd("path/to/AD-Metabolic-Collapse")

# 2) figures 1–4 + supplemental panels (data/sample/ only)
source("R/figures/Fig1_Dissociation.R"); source("R/figures/Fig2_CrossCellular.R")
source("R/figures/Fig3_Subtype_Network.R"); source("R/figures/Fig4_Clinical.R")
source("R/figures/SuppFig1_Comprehensive.R"); source("R/figures/SuppFig1_Iron_Ragulator.R")
source("R/figures/SuppFig2_Temporal.R")

# 3) figures 5–6 (set EMORY_PATH / ADNIMERGE_PATH / ADNI_AUX_PATH at top of each)
source("R/figures/Fig5_CSF_Validation.R"); source("R/figures/Fig6_Iron_Pathway.R")

# 4) tables
source("R/tables/Table1_CrossCellular.R")
source("R/tables/Supp_Table1_GeneExpression.R"); source("R/tables/Supp_Table2_PartialCor.R")
source("R/tables/Supp_Table3_Sensitivity.R"); source("R/tables/Table2_CSF_Proteomics.R")
```

## ADNI Path Configuration
For `Fig5`, `Fig6`, `Table2`, set at the top of each script:
```r
EMORY_PATH     <- "path/to/emory_results/"      # Emory TMT-MS file
ADNIMERGE_PATH <- "path/to/ADNIMERGE2/data/"    # DXSUM.rda, ADSL.rda
ADNI_AUX_PATH  <- "path/to/adni_aux/"           # Elecsys + SomaScan files
```

## Expected key values (verification)
- Transcriptomic: MCT4 −43.2% (β = −1.036); HK2 −35.2%; LDHA −20.8%; GLUT1 −7.4%;
  neuronal MCT2 (SLC16A7) preserved (+2.3%); ANLS slope β = −0.295 (p = 0.016);
  MCT4→neuronal V-ATPase donor partial r = +0.466 (p = 8.0×10⁻⁶);
  ANLS composite→neuronal V-ATPase partial r = +0.305; TFRC–ANLS r = +0.666;
  neuronal V-ATPase −5.4% vs astrocytic −0.8% (Δslope p = 0.0005).
- CSF (distribution-robust): V1A KW p = 0.999; Energy/Demand KW p = 0.285;
  V1A–Tau ladder +0.86→+0.23→+0.25→+0.05→−0.06→+0.03;
  HK1–Tau ladder +0.52→+0.42→+0.42; vs Elecsys Tau HK1 TMT +0.18 / SomaScan +0.31,
  V1A +0.03, SomaScan MAPT −0.17, SomaScan TFRC −0.44;
  TFRC ANCOVA p = 0.068 (−4.7% DEM vs CN; non-significant trend); LCN2 KW p = 0.18, ANCOVA p = 0.11 (n.s.).

## License
Code: MIT License — Data: subject to SEA-AD and ADNI data use agreements (not redistributable)
