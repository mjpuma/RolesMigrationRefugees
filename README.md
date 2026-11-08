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
config.yml                     all settings: datasets, REGE, clustering, drivers, models (main analysis)
config_*.yml                   variants: original SPEI-01 drivers, Jordan's SPEI-48 file, interim drought,
                               main specification without conflict controls (see data_prep/SPEI_README.md)
run_pipeline.R                 runs every step for every dataset listed in config.yml
validate_against_original.R    compares new role assignments with the original ones
appendix_figures.R             network diagrams and role-separation plots
data_prep/get_unhcr_flows.R    downloads the current UNHCR flow data
data_prep/build_spei_country.R country-year drought indicators from gridded SPEIbase (area-weighted)
data_prep/SPEI_README.md       drought specification: rationale, timescales, results of each variant
analysis/migrant_volume_robustness.R  migrant groups under other cut-offs, thresholds, Abel & Cohen methods
analysis/migrant_net_roles.R   regular-equivalence roles on net migrant flows (net, backbone, valued)
analysis/migrant_backbone.R    backbone roles: sensitivity to share and method; driver models
analysis/refugee_role_changes.R  drivers of year-to-year changes of refugee role
analysis/manuscript_figures.R  manuscript Figures 1-4 and SI figures (output/manuscript/)
analysis/spei_size_check.R     shows that minimum-over-grid drought indices scale with country size
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
| `spei12.nc`, `spei48.nc` | SPEIbase v2.9 gridded SPEI, 0.5°, [spei.csic.es/spei_database_2_9](https://spei.csic.es/spei_database_2_9) (`nc/spei12.nc`, `nc/spei48.nc`) |
| `ne_10m_admin_0_countries/` | Natural Earth 1:10m admin-0 countries (unzipped), [naturalearthdata.com](https://www.naturalearthdata.com/downloads/10m-cultural-vectors/10m-admin-0-countries/) |
| `new_gdp_percapita_rows.csv` | GDP per capita for five territories missing from WDI (Anguilla, British Virgin Islands, Cook Islands, Montserrat, Palestine); fills gaps only |
| `spei_country_year.csv` | written by `data_prep/build_spei_country.R` (see below) |
| `unhcr_flows_latest.csv` | Current UNHCR forced-displacement flow data, written by `data_prep/get_unhcr_flows.R` (see below) |

### Current UNHCR flow data

`UNHCR_Flow_Data_REF.csv` is an older extract (refugee arrivals, 1990–2020).
To use UNHCR's current release, run once from the repository root:

```sh
Rscript data_prep/get_unhcr_flows.R
```

This installs or updates UNHCR's `refugees` R package (from UNHCR's
R-universe, which is newer than CRAN), takes its `flows` table (one row per
origin × asylum country × year), drops rows without a valid ISO3 code on both
sides (unknown origin, stateless), and writes `data/raw/unhcr_flows_latest.csv`
plus a `_README.txt` recording the package version, download date and the
years covered by each column:

| Column | New arrivals / applications |
|---|---|
| `refugees` | refugee arrivals (group-recognised and temporary protection; what the original analysis used) |
| `asylum_seekers` | new individual asylum applications (broad coverage from 2000) |
| `oip` | other people in need of international protection (from 2018) |
| `returned_refugees` | refugee returns to the country of origin (not used by default) |

Two datasets in `config.yml` read this file:

- `unhcr_latest`: refugee arrivals only, all years in the release (same definition as `unhcr`, newer data).
- `unhcr_displacement`: refugee arrivals + asylum applications + OIP, summed, from 2000.

Add them to `run:` to include them. The driver files listed above end in
2020–2021, so later flow years have missing drivers and are left out of the
models (the run prints which drivers end early) until those files are
updated. Roles are still computed for every year.

`data_ref/country_list.csv` is the Abel & Cohen country reference list.
`data_ref/country_name_crosswalk.csv` maps country-name spellings used by the
different sources to ISO3 codes; add a row there if a name is reported as
unmatched during a run.

## Running

From the repository root, build the drought indicators once (about 20 s), then
run the pipeline and the analyses. Every script takes an optional config file:

```sh
Rscript data_prep/build_spei_country.R           # -> data/raw/spei_country_year.csv
Rscript run_pipeline.R [config.yml]
Rscript analysis/migrant_backbone.R [config.yml]
Rscript analysis/refugee_role_changes.R [config.yml]
Rscript analysis/manuscript_figures.R [config.yml]
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
    value: people             # or several columns, summed: [refugees, asylum_seekers]
    years: [2000, 2024]       # optional: keep only these years
    net: false                # optional: true = one link per pair, net sender -> net receiver
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
