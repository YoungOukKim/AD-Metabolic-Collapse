# `data/`

Third-party inputs are **not redistributed here**. Each is available from its own source
under its own data use agreement. Place them at the paths below, or point the
corresponding environment variable at wherever you keep them.

| Expected path | Environment variable | Source |
|---|---|---|
| `data/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad` | `H5AD_PATH` | SEA-AD, Allen Institute (`sea-ad.org`) |
| `data/EMORY_CSF_TMT_MS.csv` | `EMORY_CSV` | ADNI (`adni.loni.usc.edu`), Emory TMT-MS CSF proteomics |
| `data/UPENNBIOMK_ROCHE_ELECSYS_09Jan2026.csv` | `ELECSYS_CSV` | ADNI, Roche Elecsys CSF immunoassay |
| `data/ADNIMERGE2/` | `ADNI_RDA` | ADNI, ADNIMERGE2 R package data directory |
| ROSMAP DLPFC snRNA-seq | `ROSMAP_ASTRO`, `ROSMAP_CLIN` | Rush ADRC / Synapse; used only by `R/mediation/46` and `47` |

The files already in this directory (`38_all_genes.csv`, `47_chain_braak.csv`, `sample/`)
are derived summaries small enough to redistribute, not raw data.

Every script takes its paths relative to the repository root. Run them from the root:

```
Rscript R/revision2/91_revision2_S12_all.R
```
