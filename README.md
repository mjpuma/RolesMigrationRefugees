# RolesMigrationRefugees

Structural roles of countries in international **migrant** and **refugee**
mobility networks, identified with **regular equivalence (REGE)**, and their
associations with economic, conflict and environmental drivers.

This repository accompanies a manuscript in preparation, *Structural roles,
temporal dynamics, and characteristics of nation-state actors in international
mobility networks* (Zurek-Ost, Hood, Puma, Muneepeerakul & Johnson). It
contains:

* a **configurable R pipeline** (repository root) that rebuilds the analysis
  and makes it easy to add new migration datasets;
* the **original analysis scripts** (`legacy/`), kept for reference and
  reproducibility;
* a **corrected, line-by-line commented role-assessment section**
  (`for_colleagues/`) that drops into the original scripts.

## Repository layout

```
config.yml                     all settings: datasets, REGE, clustering, drivers, models
run_pipeline.R                 runs every step for every dataset listed in config.yml
validate_against_original.R    compares new role assignments with the original ones
R/
  utils.R                      helpers (config, paths, country-name matching, adjusted Rand index)
  flows.R                      reads any origin-destination flow table into one standard format
  rege.R                       REGE for binary / multi-threshold ties, and valued REGE wrapper
  roles.R                      stacked country-period network, clustering, k diagnostics, role summaries
  drivers.R                    country-year covariate panel and aggregation to each dataset's periods
  models.R                     mixed-effects logit per role (+ optional multinomial model)
  figures.R                    standard figures
src/rege_valued.cpp            White's REGGE on flow values (C++ via Rcpp)
data_ref/                      country reference list and name -> ISO3 crosswalk
data/raw/                      put the raw data here (not tracked; see below)
for_colleagues/                fixed drop-in section for the original role-assessment scripts
legacy/                        original analysis scripts (see legacy/README.md)
```

## Requirements

R ≥ 4.2 and these packages:

```r
install.packages(c("dplyr", "tidyr", "readr", "purrr", "stringr", "ggplot2",
                   "forcats", "scales", "yaml", "readxl", "igraph", "matrixStats",
                   "cluster", "lme4", "mclogit", "nnet", "Rcpp"))
```

`Rcpp` (plus a C++ compiler) is needed only for valued REGE
(`roles: method: valued`). On macOS run `xcode-select --install` once.

## Data

Raw data are **not** included: some files exceed GitHub's size limit, and
several sources do not allow redistribution. Download them into `data/raw/`
using the file names in `config.yml`:

| File (in `data/raw/`) | Source |
|---|---|
| `bilat_mig.csv` | Abel & Cohen (2022) bilateral migration flow estimates, figshare [doi:10.6084/m9.figshare.c.5800838](https://doi.org/10.6084/m9.figshare.c.5800838) (column `da_pb_closed`) |
| `UNHCR_Flow_Data_REF.csv` | UNHCR Forced Displacement Flow Dataset, refugee arrivals (obtained from UNHCR; columns `refflows.1990` … `refflows.2020`) |
| `WPP2022_Demographic_Indicators_Medium.csv` | UN World Population Prospects 2022, [population.un.org/wpp](https://population.un.org/wpp/) |
| `GDP_Global_1990_2021.csv` | World Bank WDI, GDP per capita (current US$), `NY.GDP.PCAP.CD`, DataBank CSV export |
| `spei01-annual-{min,max,mean}{,-cv}.csv` | SPEIbase v2.9 (Beguería et al.), SPEI-01 aggregated to country-year annual min / max / mean and their coefficients of variation |
| `GEDEvent_v22_1.csv` | UCDP Georeferenced Event Dataset v22.1, [ucdp.uu.se/downloads](https://ucdp.uu.se/downloads/) |
| `emdat_public_1990-2020.xlsx` | EM-DAT public export 1990–2020, [emdat.be](https://www.emdat.be/) |

`data_ref/country_list.csv` is the Abel & Cohen country reference list.
`data_ref/country_name_crosswalk.csv` maps country-name spellings used by the
different sources to ISO3 codes; add a row there if a name is reported as
unmatched during a run.

## Running

From the repository root:

```r
source("run_pipeline.R")                 # or: Rscript run_pipeline.R [config.yml]
```

Results go to `output/<dataset>/` (not tracked):

| Step | Code | Output |
|---|---|---|
| 1. Read flows | `R/flows.R` | – |
| 2. Stacked country-period network → REGE → hierarchical clustering | `R/rege.R`, `R/roles.R` | `roles_by_country_period.csv`, `role_profile.csv`, `k_diagnostics.csv`, `role_transitions.csv`, `role_stability_by_country.csv` |
| 3. Driver panel aggregated to each dataset's periods | `R/drivers.R` | `drivers_by_country_period.csv` |
| 4. Role models | `R/models.R` | `binary_role_models.csv`, `multinomial_role_model.csv`, `covariate_missingness.csv` |
| 5. Figures | `R/figures.R` | `figures/*.png` |

To compare with the original role assignments, run
`Rscript validate_against_original.R`. It reads the original outputs in
`legacy/results/` (`UNDESA_full_roles_k3.csv`, `UNHCR_full_roles_k3.csv`;
role labels only, as produced by the legacy scripts).

## Adding a new migration dataset

Add a block under `flows:` in `config.yml` and list it under `run:`. No R code
changes are needed:

```yaml
flows:
  my_new_flows:
    label: "Describe the data"
    file: data/raw/my_flows.csv
    format: long              # or wide (one column per year; see the unhcr entry)
    origin: origin_iso3       # column names in your file (ISO3 codes)
    destination: dest_iso3
    period: year
    value: people
    period_length: 1          # years covered by one record
    period_anchor: end        # "start" if the year is the first year of the interval
    roles:                    # optional per-dataset overrides of the roles settings
      thresholds: [1, 1000]
      k: auto
run: [undesa, unhcr, my_new_flows]
```

Driver data are aggregated automatically over each dataset's own time windows,
so annual, 5-year and irregular periods all work. To add a driver, write a
loader in `R/drivers.R` that returns `iso3, year, <variables>`, add it to
`build_driver_panel()`, and list the variables under `models: covariates:`.

## Key settings

| Setting | Meaning |
|---|---|
| `roles: method` | `binary` (REGE on 0/1 ties, one relation per threshold) or `valued` (White's REGGE on the flow values) |
| `roles: thresholds` | a tie exists when flow ≥ threshold; several values give multi-threshold REGE, in which ties must be matched at each magnitude |
| `roles: rege_iterations` | number of REGE iterations (default 3); report it, since valued results depend on it |
| `roles: cluster_on` | `rege_distance` (clustering on 1 − REGE) or `legacy_profiles` (reproduces the original scripts' clustering) |
| `roles: linkage`, `k`, `k_range` | hierarchical clustering linkage; fixed k or `auto` (maximum average silhouette); range scanned for diagnostics |
| `drivers: zero_fill`, `per_capita`, `aggregation` | how missing event data, population scaling and period aggregation are handled |
| `models: covariates`, `log_transform`, `scale`, `time_effect`, `lag_periods`, `binary_engine`, `multinomial` | model specification |

## Methodological notes

* **Stacked networks.** All periods are placed in one matrix of country-period
  nodes. Ties exist only within a period, but every node is compared with every
  other, so roles share one scale across time.
* **Binary ties in dense networks.** With a threshold of one person, almost
  every country both sends and receives migrants, so binary REGE scores nearly
  all pairs as fully equivalent. Use magnitude thresholds or valued REGE for
  the migrant network.
* **Valued REGE across periods.** In a stacked network, cross-period
  similarities shrink with each iteration, because countries in different
  periods never share partners. Results therefore depend on the number of
  iterations. Running REGE within each period and aligning roles across periods
  structurally is a more robust alternative, planned as a `rege_scope` option.

## Original scripts

`legacy/` holds the original `FINAL_*.R` analysis scripts by the manuscript
authors. They expect their own working directory with `data/` and `results/`
folders; see `legacy/README.md`. `for_colleagues/FIXED_role_assessment_section.R`
is a corrected replacement for the role-assessment block of those scripts, with
every line commented. `for_colleagues/run_fixed_and_compare.R` runs that section
(binary and valued ties, both datasets) and compares the result with the
original roles in `legacy/results/`; it needs the `blockmodeling` package.

## License

The pipeline code in this repository is released under the [MIT License](LICENSE).
The scripts in `legacy/` remain the work of their original authors and are
included for reference. Contact the authors before reusing them. Data are
subject to the terms of their respective providers.
