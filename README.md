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
│   └── sample/                          # Bin-level summary CSVs (included)
│       ├── astro_bin_means.csv          # Astrocyte bin-level expression means
│       ├── neuron_bin_means.csv         # Excitatory neuron bin-level means
│       ├── donor_level_summary.csv      # Donor-level (n=84) summary
│       └── astro_subtype_trajectories.csv
├── R/
│   ├── data_extraction/
│   │   ├── 01_extract_seaad.R           # Extract from SEA-AD h5ad (raw data)
│   │   └── 02_generate_sample_csvs.R    # Generate data/sample/ CSVs
│   ├── figures/
│   │   ├── utils.R                      # Shared theme & helper functions
│   │   ├── Fig1_Dissociation.R          # Fig 1 (a-d)
│   │   ├── Fig2_CrossCellular.R         # Fig 2 (a-f)
│   │   ├── Fig3_Subtype_Network.R       # Fig 3 (a-c)
│   │   ├── Fig4_Clinical.R              # Fig 4 (a-f)
│   │   ├── Fig5_CSF_Validation.R        # Fig 5 (a-h)
│   │   ├── Fig6_Iron_Suppression.R      # Fig 6 (a)
│   │   ├── ED_Fig1_Comprehensive.R
│   │   ├── ED_Fig2_Iron_Ragulator.R
│   │   ├── ED_Fig3_Temporal.R
│   │   └── ED_Fig4_Donor_Partial.R
│   └── tables/
│       ├── Table1_CrossCellular.R
│       └── Supp_Table3_Sensitivity.R
└── docs/
    └── data_dictionary.md               # Variable definitions and CSV schemas
```

---

## Data

### SEA-AD (not included — requires data access agreement)

- File: `SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad`
- Access: [portal.brain-map.org](https://portal.brain-map.org)
- n = 1,378,211 nuclei, 84 donors

### ADNI Emory TMT-MS CSF Proteomics (not included — requires ADNI access)

- Access: [adni.loni.usc.edu](https://adni.loni.usc.edu)
- n = 1,105 subjects, 3,907 proteins

### Sample Data (included)

`data/sample/` contains **bin-level summary statistics** derived from SEA-AD.
These are sufficient to reproduce all figures without the raw h5ad file.
See `docs/data_dictionary.md` for variable definitions.

---

## Reproducing Figures

### Step 1: (Optional) Re-extract from SEA-AD raw data

```r
# Only needed if you have SEA-AD h5ad access and want to regenerate CSVs
source("R/data_extraction/01_extract_seaad.R")  # generates per-cell CSVs
source("R/data_extraction/02_generate_sample_csvs.R")  # generates data/sample/
```

### Step 2: Generate figures

```r
# All scripts use relative paths from the repository root
source("R/figures/Fig1_Dissociation.R")
source("R/figures/Fig2_CrossCellular.R")
source("R/figures/Fig3_Subtype_Network.R")
source("R/figures/Fig4_Clinical.R")
# Fig 5 & 6 require ADNI proteomics access
source("R/figures/Fig5_CSF_Validation.R")
source("R/figures/Fig6_Iron_Suppression.R")
```

Output figures are saved to `output/figures/`.

---

## Dependencies

```r
install.packages(c("ggplot2", "dplyr", "tidyr", "patchwork", "scales", "data.table"))

# rhdf5 required only for 01_extract_seaad.R
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install("rhdf5")
```

R >= 4.3.0

---

## License

Code: MIT License  
Data: Subject to SEA-AD and ADNI data use agreements (not redistributable)
