# Configuration

How to change settings and add a dataset. Back to the [README](../README.md).

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
| `models: covariates`, `log_transform`, `scale`, `time_effect`, `lag_periods`, `binary_engine`, `multinomial` | model specification; `binary_engine` is `glmmTMB` (main analysis), `glmer` or `mblogit` |
| `models: within_between` | covariates replaced by a within-country term (`<var>_within`, deviation from the country's mean) and a between-country term (`<var>_between`, the country's mean) |
| `models: baseline` | variables copied as `<var>_base` = the country's value in its first period with data (time-invariant), e.g. `gdp_pc_base` |
| `models: lags` | variables copied as `<var>_lag1` = value in the previous period (consecutive periods only), e.g. drought in year t-1 for annual refugee data |
| `models: income_groups` | `{var: gdp_pc_base, n: 3}`: countries are split into baseline-income terciles (cut across countries, not rows); adds tercile dummies `income_g2`, `income_g3` |
| `models: slopes_by_group` | variables whose slope is estimated separately in each income group (`<var>_x_g1..3`, replacing the common slope); a Wald test of equal slopes is added to the results (`wald_equal_slopes:<var>`, chi-square in `estimate`, df in `std.error`) |
| `models: interactions` | pairs of covariates whose product of z-scores is added (e.g. drought x baseline agriculture share) |
| `drivers: long_files` | extra long-format files (`iso3, year, <variables>`) joined into the driver panel, e.g. `data/raw/wdi_agri_share.csv` |
| `flows: <dataset>: models:` | per-dataset overrides of any `models:` entry (replaces that entry; used for the annual refugee data) |

Mixed models are fitted by `fit_mixed()` (`R/models.R`) in glmmTMB (Laplace
approximation). A fit that does not converge or has a non-positive-definite
Hessian goes to `fit_glmer_robust()`: lme4 with bobyqa, then Nelder-Mead, then
nloptwrap, and ordinary logistic regression without the country random intercept
if no acceptable mixed fit is found. The `engine` column of every results table
says which was used. All main models in the manuscript converge in glmmTMB.

Alternative configurations shipped with the repository:

| Config | Purpose |
|---|---|
| `config.yml` | main analysis: area-weighted drought indicators, drought slopes by baseline-income tercile, within- and between-country terms for GDP per capita and conflict |
| `configs/config_blended.yml` | as main, with one coefficient per driver (the earlier model, compared in Table S9); writes to `output_blended/` |
| `configs/config_nocontrols_conflict.yml` | as main, without the conflict variables (total rather than direct drought association) |
| `configs/config_drought_agriculture.yml` | drought x baseline agriculture share of GDP instead of income terciles |
| `configs/config_original_spei01.yml`, `configs/config_spei_jordan.yml`, `configs/config_drought_interim.yml` | earlier drought measures, for comparison (Table S1) |
