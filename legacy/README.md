# Original analysis scripts

The original R scripts for the manuscript, by M. A. Zurek-Ost, J. H. Hood and
colleagues (University of Florida). They are kept unchanged for reference and
reproducibility.

| Script | Purpose |
|---|---|
| `FINAL_undesa_role_assessment.R` | Migrant network (Abel & Cohen / UN DESA): REGE, clustering, role assignment, network metrics |
| `FINAL_unhcr_role_assessment.R` | Refugee network (UNHCR): same steps |
| `FINAL_models_undesa_roles.R` | Driver data and role models for migrant roles |
| `FINAL_models_unhcr_roles.R` | Driver data and role models for refugee roles |

## Running them

* Each script starts with `setwd("D:/...")`. Replace it with a folder that
  contains `data/` (the raw files listed in the top-level README) and an empty
  `results/` folder, e.g. this `legacy/` folder.
* Required packages include `igraph`, `sna`, `blockmodeling`, `factoextra`,
  `NbClust`, `clValid`, `xUCINET`, `mclogit`, `tnam` (GitHub:
  `leifeld/tnam`), `readxl`, `tidyverse`, `texreg` and `matrixStats`.
* The role-assessment scripts write `results/UNDESA_full_roles_k3.csv` and
  `results/UNHCR_full_roles_k3.csv`, which the model scripts and
  `../validate_against_original.R` read. The original copies of both files
  (UNHCR version dated 27 Jan 2025) are included in `results/`.

## Known issues

* In the role-assessment scripts, `hcut(x = a9, method = link)` ignores
  `method` (not an `hcut()` argument) and treats the dissimilarity matrix as
  data. The final roles therefore come from Ward.D2 clustering on Euclidean
  distances between rows, not average linkage on REGE distances.
  `../for_colleagues/FIXED_role_assessment_section.R` is a commented drop-in
  replacement for that block.
* The model scripts do not zero-fill missing conflict/disaster values and leave
  covariates unscaled (scaling is commented out), which can drop many
  observations and destabilise the estimates. The pipeline in the repository
  root handles both explicitly (`drivers: zero_fill`, `models: scale`).
