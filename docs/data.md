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
| `unhcr_flows_latest.csv` | UNHCR forced displacement flow data, current release; written by `data_prep/get_unhcr_flows.R` (see below). Default refugee data: column `refugees`, 1990–2020 |
| `legacy/UNHCR_Flow_Data_REF.csv` | Older UNHCR extract of refugee arrivals used in the first drafts (dataset `unhcr_legacy`); all of its flows are in the current release unchanged |
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

### UNHCR flow data

The refugee data come from UNHCR's current release. Run once from the repository root:

```sh
Rscript data_prep/get_unhcr_flows.R           # installs UNHCR's `refugees` R package
Rscript data_prep/get_unhcr_flows.R github    # or: reads the data files from UNHCR's GitHub repository
```

This takes the package's `flows` table (one row per origin × asylum country × year),
adds resettlement and returns from its `solutions` table, drops rows without a
country code on both sides (unknown origin, stateless, Tibetan), and writes
`data/raw/unhcr_flows_latest.csv` plus a `_README.txt` recording the package
version, download date and the years covered by each column:

| Column | New arrivals / applications |
|---|---|
| `refugees` | refugee arrivals: group (prima facie) recognitions and temporary protection |
| `asylum_seekers` | new individual asylum applications (20–37 reporting countries before 2000, 137 from 2000) |
| `refugee_like` | people in refugee-like situations (not included in `refugees`) |
| `oip` | other people in need of international protection (from 2018) |
| `resettlement` | refugees resettled to a third country (not used by default) |
| `returned_refugees` | refugee returns to the country of origin (not used by default) |

`refugees` and `asylum_seekers` do not double count: the first is group recognition,
the second individual applications. High-income asylum countries record almost
no `refugees`; their inflows are asylum applications.

Three datasets in `config.yml` read UNHCR data:

- `unhcr` (default): refugee arrivals, 1990–2020.
- `unhcr_displacement`: refugee arrivals + refugee-like + asylum applications + OIP, 2000–2020. Add it to `run:` to include it.
- `unhcr_legacy`: the older extract `data/raw/legacy/UNHCR_Flow_Data_REF.csv` used in the first drafts. All 2,444 of its flows are in the current release with the same values, and the roles and models are identical.

`data_ref/country_list.csv` is the Abel & Cohen country reference list.
`data_ref/country_name_crosswalk.csv` maps country-name spellings used by the
different sources to ISO3 codes; add a row there if a name is reported as
unmatched during a run.
