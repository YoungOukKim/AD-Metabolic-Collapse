# Data Dictionary

## data/sample/ — Bin-level summary CSVs (derived from SEA-AD)

These files contain **bin-level mean expression values** (per pseudo-progression bin)
derived from the full SEA-AD single-nucleus RNA-seq dataset. They are sufficient to
reproduce all figures in the paper without requiring access to the raw h5ad file.

---

### astro_bin_means.csv

Per-bin mean log-normalized expression across all astrocytes (n = 67,419).

| Column | Description |
|--------|-------------|
| `bin` | Pseudo-progression bin (0.2 to 0.9, step 0.1) |
| `SLC16A3` | MCT4 — primary astrocytic lactate exporter (**ANLS hub**) |
| `SLC16A1` | MCT1 — additional lactate transporter |
| `SLC2A1`  | GLUT1 — glucose importer |
| `LDHA`    | Lactate dehydrogenase A (pyruvate → lactate) |
| `LDHB`    | Lactate dehydrogenase B (lactate → pyruvate) |
| `HK1`, `HK2` | Hexokinases |
| `PDK1`    | Pyruvate dehydrogenase kinase 1 |
| `PKM`     | Pyruvate kinase M |
| `ATP6V1A`, `ATP6V1B2`, `ATP6V1C1`, `ATP6V1E1`, `ATP6V1H` | V-ATPase V1 sector |
| `ATP6V0A1`, `ATP6V0B`, `ATP6V0C`, `ATP6V0D1`, `ATP6V0E1` | V-ATPase V0 sector |
| `LAMP1`, `LAMP2` | Lysosomal membrane proteins |
| `CTSB`, `CTSD`   | Cathepsins B and D |
| `TFRC`    | Transferrin receptor (iron uptake) |
| `FTH1`, `FTL` | Ferritin heavy/light chain |
| `CP`      | Ceruloplasmin (ferroxidase) |
| `VDAC1`   | Mitochondria-lysosome contact marker |
| `LAMTOR1`–`LAMTOR5` | Ragulator complex (mTOR scaffold) |
| `SLC9A6`  | NHE6 — endosomal pH regulator |
| `PTGDS`   | Prostaglandin D2 synthase (metabolic checkpoint) |
| `LCN2`    | Lipocalin-2 (inflammatory effector) |
| `SLC1A2`  | EAAT2 (glutamate transporter) |
| `ATP1A2`  | Na⁺/K⁺-ATPase α2 (glutamate clearance energy) |
| `ANLS`    | Composite: mean(SLC2A1, LDHA, SLC16A1) |
| `VATpase` | Composite: mean of 10 V-ATPase subunits |
| `MCT4`    | Alias for SLC16A3 |
| `ED_energy` | Composite: mean(SLC2A1, LDHA, SLC16A1, PKM, HK1) |
| `ED_ratio`  | ED_energy / VATpase |

---

### neuron_bin_means.csv

Per-bin mean expression across excitatory neurons (n = 671,689).
Subclasses: L2/3 IT, L4 IT, L5 IT, L5 ET, L5/6 NP, L6 IT, L6 IT Car3, L6 CT, L6b.

| Column | Description |
|--------|-------------|
| `bin` | Pseudo-progression bin |
| `SLC16A7` | MCT2 — high-affinity neuronal lactate importer |
| `LDHB`    | Oxidizes lactate → pyruvate (neuronal ATP) |
| `LAMP1`   | Lysosomal membrane protein |
| `NDUFS1`  | Mitochondrial Complex I subunit |
| `ATP6V1A`, ... | V-ATPase subunits (neuronal) |
| `VATPase` | Composite: mean(ATP6V1A, ATP6V1B2, ATP6V0A1, ATP6V0C, ATP6V0D1, ATP6V1E1) |

---

### donor_level_summary.csv

Per-donor (n = 84) aggregated statistics from SEA-AD astrocytes.

| Column | Description |
|--------|-------------|
| `donor_id` | Donor identifier (rounded CPS) |
| `ANLS` | Mean ANLS composite |
| `MCT4` | Mean MCT4 expression |
| `VATpase` | Mean V-ATPase composite (astrocyte) |
| `VATPase_n` | Mean V-ATPase composite (neuron) |
| `mean_cps` | Mean continuous pseudo-progression score |
| `n_astro` | Number of astrocytes for this donor |
| `braak_num` | Braak stage (numeric: 0–6) |
| `abc_score` | NIA-AA ABC neuropathological score |
| `cognitive` | Cognitive status (Normal/MCI/Dementia) |

---

### astro_subtype_trajectories.csv

ANLS trajectory by astrocyte supertype (Astro_1 – Astro_6).

| Column | Description |
|--------|-------------|
| `bin` | Pseudo-progression bin |
| `subtype` | Astrocyte supertype label |
| `ANLS` | Mean ANLS composite for this subtype × bin |
| `n_cells` | Cell count |

---

## ADNI CSF datasets (not included; require ADNI access)

Access via: https://adni.loni.usc.edu

### Emory TMT-MS CSF proteomics — primary CSF proteomics

File: merged proteomics matrix (e.g. `merged_proteomics_dx.rds` after the
`utils.R` DXSUM/ADSL merge). n = 1,105 subjects (CN = 379, MCI = 562, DEM = 164),
3,907 proteins. Below-detection values are encoded as −4 and converted to NA on
load. Protein columns follow the `GeneName_UniProtID` convention:

- `ATP6V1A_P38606` — V-ATPase V1A subunit
- `ATP6V1E1_P36543` — V-ATPase V1E1 subunit
- `MAPT_P10636` — Tau protein
- `GFAP_P14136` — Glial fibrillary acidic protein
- `TREM2_Q9NZC2` — Triggering receptor on myeloid cells 2
- `TFRC_P02786` — Transferrin receptor
- `NEFL_P07196` — Neurofilament light chain (NfL)
- `HK1_P19367`, `LDHA_P00338`, `PKM_P14618` — Glycolytic enzymes
- `APP_P05067`, `LCN2_P80188`, `LAMP2_P13473`, `CTSB_P07858` — additional markers

### Roche Elecsys immunoassay — orthogonal total-Tau validation

File: `UPENNBIOMK_ROCHE_ELECSYS_*.csv`. Loaded by `load_elecsys_tau()`, which
returns one baseline total-Tau value per subject:

| Column | Description |
|--------|-------------|
| `RID` | ADNI subject ID (merge key) |
| `ElecsysTau` | Roche Elecsys total tau (immunoassay), baseline visit |

Used to anchor CSF protein–Tau associations against a non-mass-spectrometry
detection chemistry (Fig 5 d/g, Fig 6 b, Table 2 Part A).

### SomaScan 7k — independent aptamer platform

File: SOMAscan7k post-QC protein matrix (ADNI Cruchaga-lab). Loaded by
`load_somascan()`, which returns the seqID columns:

| Returned column | Default seqID | Protein |
|-----------------|---------------|---------|
| `HK1_s`  | `X13131.5` | HK1 (hexokinase-1) |
| `MAPT_s` | `X5854.60` | MAPT (Tau) |
| `TFRC_s` | `X6895.1`  | TFRC (transferrin receptor) |

SeqIDs are dataset-version specific — edit the `seq` argument of
`load_somascan()` if your matrix uses different identifiers. **V-ATPase V1A is
not present on the SomaScan panel**, so V1A cannot be cross-validated on this
platform.
