# AD-Metabolic-Collapse

**Astrocytic lactate shuttle disruption is linked to cross-cellular lysosomal energetic decoupling in Alzheimer's disease**

---

## Overview

This repository contains all R code and sample data used to reproduce the figures and tables in the manuscript. The analysis integrates single-nucleus transcriptomics (SEA-AD atlas) with CSF proteomics (ADNI Emory TMT-MS) to define an **'energy-starved lysosome' (ESL)** state in Alzheimer's disease.

---

## Repository Structure

```
AD-Metabolic-Collapse/
├── README.md
├── data/
│   └── sample/                                   # Bin-level summary CSVs (included)
│       ├── astro_bin_means.csv                   # Astrocyte bin-level expression means
│       ├── neuron_bin_means.csv                  # Excitatory neuron bin-level means
│       ├── donor_level_summary.csv               # Donor-level (n=84) summary
│       └── astro_subtype_trajectories.csv        # Astrocyte subtype ANLS trajectories
├── R/
│   ├── data_extraction/
│   │   ├── 01_extract_seaad.R                    # Extract from SEA-AD h5ad (raw data)
│   │   └── 02_generate_sample_csvs.R             # Generate data/sample/ CSVs from per-cell data
│   ├── figures/
│   │   ├── utils.R                               # Shared theme, helpers, ADNI loader
│   │   ├── Fig1_Dissociation.R                   # Fig 1 (a-d): Energetic dissociation
│   │   ├── Fig2_CrossCellular.R                  # Fig 2 (a-f): Cross-cellular propagation
│   │   ├── Fig3_Subtype_Network.R                # Fig 3 (a-c): Subtype consistency + network
│   │   ├── Fig4_Clinical.R                       # Fig 4 (a-f): Donor-level clinical validation
│   │   ├── Fig5_CSF_Validation.R                 # Fig 5 (a-h): CSF proteomic validation *
│   │   ├── Fig6_Iron_Suppression.R               # Fig 6 (a): Iron suppression effect *
│   │   ├── ED_Fig1_Comprehensive.R               # Extended Data Fig 1: All gene changes
│   │   ├── ED_Fig2_Iron_Ragulator.R              # Extended Data Fig 2: Iron + Ragulator
│   │   ├── ED_Fig3_Temporal.R                    # Extended Data Fig 3: Temporal ordering
│   │   └── ED_Fig4_Donor_Partial.R               # Extended Data Fig 4: Donor partial cor.
│   └── tables/
│       ├── Table1_CrossCellular.R                # Table 1: Cross-cellular correlations
│       ├── Table2_CSF_Proteomics.R               # Table 2: V-ATPase subunit dissociation *
│       ├── Supp_Table1_GeneExpression.R          # Supp Table 1: Gene expression (179 genes)
│       ├── Supp_Table2_PartialCor.R              # Supp Table 2: Complete partial correlation
│       └── Supp_Table3_Sensitivity.R             # Supp Table 3: Bin 0.1 sensitivity analysis
└── docs/
    └── data_dictionary.md                        # Variable definitions and CSV schemas
```

`*` requires ADNI data access (see Data section below)

---

## Data Requirements

### SEA-AD — not included, requires data access agreement

- File: `SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad`
- Access: [portal.brain-map.org](https://portal.brain-map.org)
- n = 1,378,211 nuclei, 84 donors

### ADNI Emory TMT-MS CSF Proteomics — not included, requires ADNI access

Used by: `Fig5_CSF_Validation.R`, `Fig6_Iron_Suppression.R`, `Table2_CSF_Proteomics.R`

- Access: [adni.loni.usc.edu](https://adni.loni.usc.edu)
- Search "Emory" or "TMT" in the ADNI study data inventory
- n = 1,105 subjects, 3,907 proteins
- Also requires `DXSUM.rda` and `ADSL.rda` from the ADNIMERGE2 data package

### Sample Data — included in this repository

`data/sample/` contains **bin-level summary statistics** (mean expression per pseudo-progression bin) derived from SEA-AD. These are sufficient to reproduce **all figures except Fig 5–6** and **all tables except Table 2** without requiring access to the raw datasets.

See `docs/data_dictionary.md` for full variable definitions.

---

## Reproducing the Analysis

### Step 0 — Install dependencies

```r
install.packages(c("ggplot2", "dplyr", "tidyr", "patchwork", "scales", "data.table"))

# rhdf5 required only for 01_extract_seaad.R
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install("rhdf5")
```

R >= 4.3.0

### Step 1 — Set working directory to repository root

```r
setwd("path/to/AD-Metabolic-Collapse")
```

All scripts use **relative paths** from the repository root.

### Step 2 — (Optional) Re-extract from SEA-AD raw data

Skip this step if you use the provided `data/sample/` CSVs.

```r
# Requires SEA-AD h5ad file — edit path inside the script first
source("R/data_extraction/01_extract_seaad.R")   # generates per-cell CSVs
source("R/data_extraction/02_generate_sample_csvs.R")  # generates data/sample/
```

### Step 3 — Generate figures

```r
# Figures 1–4 and Extended Data: data/sample/ CSVs only
source("R/figures/Fig1_Dissociation.R")
source("R/figures/Fig2_CrossCellular.R")
source("R/figures/Fig3_Subtype_Network.R")
source("R/figures/Fig4_Clinical.R")
source("R/figures/ED_Fig1_Comprehensive.R")
source("R/figures/ED_Fig2_Iron_Ragulator.R")
source("R/figures/ED_Fig3_Temporal.R")
source("R/figures/ED_Fig4_Donor_Partial.R")

# Figures 5–6: requires ADNI data — set paths at top of each script
source("R/figures/Fig5_CSF_Validation.R")
source("R/figures/Fig6_Iron_Suppression.R")
```

Output figures are saved to `output/figures/`.

### Step 4 — Generate tables

```r
# Tables 1, Supp 1–3: data/sample/ CSVs only
source("R/tables/Table1_CrossCellular.R")
source("R/tables/Supp_Table1_GeneExpression.R")
source("R/tables/Supp_Table2_PartialCor.R")
source("R/tables/Supp_Table3_Sensitivity.R")

# Table 2: requires ADNI data — set paths at top of the script
source("R/tables/Table2_CSF_Proteomics.R")
```

Output tables are saved to `output/tables/`.

---

## ADNI Path Configuration

For scripts requiring ADNI data (`Fig5`, `Fig6`, `Table2`), set these two lines at the top of each script:

```r
EMORY_PATH     <- "path/to/emory_results/"         # folder with Emory TMT-MS file
ADNIMERGE_PATH <- "path/to/ADNIMERGE2/data/"       # folder with DXSUM.rda, ADSL.rda
```

The `load_adni_proteomics()` function in `utils.R` handles merging of proteomics data with diagnosis (DXSUM) and demographics (ADSL) automatically.

---

## What Each Script Produces

| Script | Output | Data needed |
|--------|--------|-------------|
| `Fig1_Dissociation.R` | Fig 1 (a–d) | `data/sample/` |
| `Fig2_CrossCellular.R` | Fig 2 (a–f) | `data/sample/` |
| `Fig3_Subtype_Network.R` | Fig 3 (a–c) | `data/sample/` |
| `Fig4_Clinical.R` | Fig 4 (a–f) | `data/sample/` |
| `Fig5_CSF_Validation.R` | Fig 5 (a–h) | ADNI |
| `Fig6_Iron_Suppression.R` | Fig 6 (a) | ADNI |
| `ED_Fig1–4.R` | Extended Data Figs 1–4 | `data/sample/` |
| `Table1_CrossCellular.R` | Table 1 | `data/sample/` |
| `Table2_CSF_Proteomics.R` | Table 2 (Part A+B) | ADNI |
| `Supp_Table1_GeneExpression.R` | Supp Table 1 | `data/sample/` |
| `Supp_Table2_PartialCor.R` | Supp Table 2 | `data/sample/` |
| `Supp_Table3_Sensitivity.R` | Supp Table 3 | `data/sample/` |

---

## License

Code: MIT License  
Data: Subject to SEA-AD and ADNI data use agreements (not redistributable)
