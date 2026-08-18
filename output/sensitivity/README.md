# `output/sensitivity/`

This directory holds the donor-level intermediates that several scripts read rather than
recompute. They are **not committed, and will not be**. `donor_by_gene.rds` is a pair of dense
84 x 36,601 matrices computed directly from the SEA-AD h5ad, so committing it would
redistribute that dataset in derived form. It also adds nothing: anyone running this code
needs the h5ad regardless, and the script that builds this file is in the repository.

| File | Produced by | Read by |
|---|---|---|
| `donor_by_gene.rds` | `R/sensitivity/01_extract_h5ad.R` | `R/revision2/90_revision2_all.R` (C1 consistency gate, S6 pathway screen), `R/audit/78`, `79`, `80` |
| `detection_matched_null_all_genes.csv` | `R/audit/77_detection_matched_null_standalone.R` | `R/audit/04` |

To regenerate:

```
Rscript R/sensitivity/01_extract_h5ad.R      # writes donor_by_gene.rds here
Rscript R/revision2/90_revision2_all.R
Rscript R/revision2/91_revision2_S12_all.R
```

`R/revision2/91_revision2_S12_all.R` does **not** require `donor_by_gene.rds`: it rebuilds
the donor x gene mean and detection matrices itself in the same h5ad pass, and
cross-checks them against `donor_by_gene.rds` when that file is present.
