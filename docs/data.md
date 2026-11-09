# Requirements and data

Packages needed and where to get each input file. Back to the [README](../README.md).

## Requirements

R ≥ 4.2 and these packages:

```r
install.packages(c("dplyr", "tidyr", "readr", "purrr", "stringr", "ggplot2",
                   "forcats", "scales", "yaml", "readxl", "igraph", "matrixStats",
                   "cluster", "lme4", "mclogit", "nnet", "Rcpp", "glmmTMB",
                   "patchwork", "ggrepel", "ragg", "sf", "sandwich", "survival", "foreign",
                   "brglm2", "rstanarm", "loo", "bayesplot"))
```

`glmmTMB` fits the main models. `sandwich`, `survival`, `foreign` and `brglm2` are
used by `analysis/regression_robustness.R`; `rstanarm`, `loo` and `bayesplot` only by
`analysis/regression_bayes.R`.

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
| `wdi_agri_share.csv` | World Bank WDI, agriculture, forestry and fishing value added (% of GDP), `NV.AGR.TOTL.ZS`, reshaped to `iso3, year, agri_share` (only for `configs/config_drought_agriculture.yml`) |
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
