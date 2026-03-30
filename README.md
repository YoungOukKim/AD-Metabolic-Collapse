# AD-Metabolic-Collapse

**Astrocytic lactate shuttle disruption is linked to cross-cellular lysosomal energetic decoupling in Alzheimer's disease**

YoungOuk Kim, WooMyung Heo, Se Jin Park, YoungChul Kim  
BioXP Research Institute / Kangwon National University

---

## Paper

> Kim et al. (2025). *Astrocytic lactate shuttle disruption is linked to cross-cellular lysosomal energetic decoupling in Alzheimer's disease.* Nature Aging (submitted).

---

## Repository Structure

```
AD-Metabolic-Collapse/
├── README.md
├── data/
│   └── sample/                  # Sample CSVs derived from SEA-AD (bin-level means)
│       ├── astro_bin_means.csv          # Astrocyte bin-level expression means
│       ├── neuron_bin_means.csv         # Excitatory neuron bin-level expression means
│       ├── donor_level_summary.csv      # Donor-level (n=84) summary statistics
│       └── astro_subtype_trajectories.csv  # Astrocyte subtype ANLS trajectories
├── R/
│   ├── data_extraction/
│   │   └── 01_extract_seaad.R       # Extract gene expression from SEA-AD h5ad
│   ├── figures/
│   │   ├── Fig1_Dissociation.R          # Fig 1 (a-d): Energetic dissociation
│   │   ├── Fig2_CrossCellular.R         # Fig 2 (a-f): Cross-cellular propagation
│   │   ├── Fig3_Subtype_Network.R       # Fig 3 (a-c): Subtype + Network
│   │   ├── Fig4_Clinical.R              # Fig 4 (a-f): Donor-level clinical
│   │   ├── Fig5_CSF_Validation.R        # Fig 5 (a-h): CSF proteomic validation
│   │   ├── Fig6_Iron_Suppression.R      # Fig 6 (a): Iron suppression effect
│   │   ├── ED_Fig1_Comprehensive.R      # Extended Data Fig 1
│   │   ├── ED_Fig2_Iron_Ragulator.R     # Extended Data Fig 2
│   │   ├── ED_Fig3_Temporal.R           # Extended Data Fig 3
│   │   └── ED_Fig4_Donor_Partial.R      # Extended Data Fig 4
│   └── tables/
│       ├── Table1_CrossCellular.R       # Table 1: Cross-cellular correlations
│       └── Supp_Table3_Sensitivity.R    # Supplementary Table 3: Bin 0.1 sensitivity
└── docs/
    └── data_dictionary.md               # Variable definitions and CSV schemas
```

---

## Data

### SEA-AD (not included — requires data access agreement)

The primary analysis uses the **Seattle Alzheimer's Disease Brain Cell Atlas** (SEA-AD) single-nucleus RNA-seq dataset:

- File: `SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad`  
- Access: [portal.brain-map.org](https://portal.brain-map.org)  
- n = 1,378,211 nuclei, 84 donors

### ADNI Emory TMT-MS CSF Proteomics (not included — requires ADNI access)

- Access: [adni.loni.usc.edu](https://adni.loni.usc.edu)  
- n = 1,105 subjects, 3,907 proteins

### Sample Data (included)

`data/sample/` contains **bin-level summary statistics** (mean expression per pseudo-progression bin) derived from SEA-AD. These are sufficient to reproduce all figures. See `docs/data_dictionary.md` for variable definitions.

---

## Reproducing Figures

### Step 1: Data extraction from SEA-AD (if you have access)

```r
source("R/data_extraction/01_extract_seaad.R")
```

This generates the full per-cell CSVs. **If you do not have SEA-AD access**, skip this step and use the pre-computed `data/sample/` CSVs directly.

### Step 2: Generate figures individually

```r
# Set your output path at the top of each script
source("R/figures/Fig1_Dissociation.R")
source("R/figures/Fig2_CrossCellular.R")
# ... etc.
```

All scripts use **relative paths** from the repository root. Set `data_path` and `save_path` at the top of each script.

---

## Dependencies

```r
install.packages(c("ggplot2", "dplyr", "tidyr", "patchwork", "scales", "rhdf5"))
```

- R >= 4.3.0  
- `rhdf5` required only for `01_extract_seaad.R` (Bioconductor)

```r
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install("rhdf5")
```

---

## Citation

```
Kim YO, Heo WM, Park SJ, Kim YC. Astrocytic lactate shuttle disruption is linked 
to cross-cellular lysosomal energetic decoupling in Alzheimer's disease. 
Nature Aging (2025, submitted).
```

---

## License

Code: MIT License  
Data: Subject to SEA-AD and ADNI data use agreements (not redistributable)
