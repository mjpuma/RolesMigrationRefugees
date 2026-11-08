# Drought indicators for the roles analysis

Three drought specifications are available. Each has its own config file and
output folder, so they can be run and compared side by side:

```sh
Rscript run_pipeline.R config_<name>.yml
Rscript analysis/migrant_backbone.R config_<name>.yml
Rscript analysis/refugee_role_changes.R config_<name>.yml
```

| Config | Drought measure | Status |
|---|---|---|
| `config.yml` (main) | Area-weighted SPEI-12 drought extent and SPEI-48 from `build_spei_country.R`, matched to each dataset's timescale, drought x income | Main analysis since 7 Nov 2026. Needs SPEIbase `spei12.nc`, `spei48.nc` and Natural Earth 1:10m countries in `data/raw/`. |
| `config_nocontrols_conflict.yml` | As main, without conflict variables | Total drought association (drought may act through conflict). |
| `config_original_spei01.yml` | SPEI-01 annual minimum, mean and CV (original files) | First manuscript draft. The minimum scales with country size (see below). |
| `config_spei_jordan.yml` | Jordan's SPEI-48 annual average of monthly minima (`fromJordan/`), standardized within country | Kept as delivered, plus the within-country fix. Adds ~28 small territories. |
| `config_drought_interim.yml` | SPEI-01 annual country mean, previous-year lag for refugees, drought × income | Runs now. Interim measure for the specification below. |

## Why the minimum over grid cells is not used

SPEI is standardized in every grid cell. A country's minimum over its cells
therefore depends on how many cells it has: larger countries have lower minima.
Across countries, the long-run mean of the minimum-based files correlates with
log land area at rho = -0.80 (SPEI-48 avg-min) and -0.90 (SPEI-01 min); the
SPEI-01 mean does not (rho = -0.02). In the models, the minimum-based index
appears as a "drought" effect for roles that differ mainly in country size
(e.g. migrant intermediates). An area-weighted mean, or the area share in
drought, has no such dependence (checked on synthetic grids: rho = 0.02).

## The main specification (`config.yml`)

**Framework.** International moves require both the aspiration to move and
the capability to do so (Carling 2002; de Haas 2021). Drought can raise
aspirations (lost harvests and income) while lowering capabilities (lost
savings and assets), so its net effect on cross-border mobility depends on
income: environmental effects on migration are weaker or reversed in the
poorest countries and strongest in middle-income, agriculture-dependent ones
(Black et al. 2011; Cattaneo and Peri 2016; Hoffmann et al. 2020). The models
therefore include drought x GDP per capita.

**Spatial aggregation.** Each month, the country value is the area-weighted
mean of SPEI over all grid cells that overlap the country (cell area x share of
the cell inside the country, so small islands are covered), and the
area-weighted share of the country with SPEI <= -1 (drought) or <= -1.5
(severe). A population-weighted version is produced if a population raster is
supplied, which targets exposed livelihoods rather than land.

**Timescales.**

| | Refugees (annual flows) | Migrants (5-year periods) |
|---|---|---|
| Response time | Displacement follows within 0-2 years, largely through conflict and loss of livelihoods | International migration responds slowly and cumulatively |
| Drought measure | SPEI-12 drought extent in year t and t-1 | SPEI-12 drought extent averaged over the 5 years (cumulative exposure), plus SPEI-48 at the start of the period (multi-year deficit going in) |
| Interaction | extent(t) x GDP, extent(t-1) x GDP | extent x GDP |

SPEI-12 summarizes the water balance of a year; SPEI-48 summarizes four years,
which matches the 5-year migrant periods. The 48-month scale is therefore a
sensible choice for migrants, provided it is aggregated by area rather than by
minimum.

**Conflict as a pathway.** Drought may affect refugee flows through conflict.
Models with conflict controls estimate the direct association; dropping the
conflict variables gives the total association. Report both.

## Drought by baseline income and by agriculture share (`config.yml`, 8 Nov 2026)

The main configuration now gives each drought indicator a separate slope in
each tercile of baseline GDP per capita (the country's GDP in its first period
with data; terciles cut across countries), with tercile main effects, instead
of the drought x current-GDP product. A Wald test (2 df) asks whether the three
slopes are equal. `config_drought_agriculture.yml` instead interacts drought
with the baseline agriculture share of GDP (World Bank NV.AGR.TOTL.ZS,
`data/raw/wdi_agri_share.csv`). See `R/models.R` for the equations.

* Refugee and migrant (backbone) roles: 1 of 30 tercile slopes has p < 0.05
  (sender-only refugee countries, middle tercile, OR 0.64, p = 0.031); no
  Wald test rejects equal slopes (p = 0.11-0.97).
* Refugee role changes: receiver-only countries start sending more often in
  drought years in the lowest tercile (OR 1.62, p = 0.048; highest tercile
  1.57, p = 0.14; Wald p = 0.42) and when drought extent rises from the
  previous year (1.27, p = 0.047).
* Agriculture: no drought x agriculture interaction for any role (closest
  receiver-only refugees, 1.23, p = 0.09); for countries that start hosting,
  agriculture share 0.61 (p = 0.04) and drought x agriculture 0.64 (p = 0.05).
* Without conflict controls the drought results are the same.
* Current GDP and baseline tercile are correlated (Spearman 0.77 refugees,
  0.87 migrants), so GDP odds ratios have wider intervals in these models;
  joint Wald tests of GDP and the tercile terms are significant for every role.

## Results with the area-weighted SPEI-12 / SPEI-48 indicators (`config.yml`, 7 Nov 2026)

Inputs: SPEIbase v2.9 `spei12.nc`, `spei48.nc`; Natural Earth 1:10m admin-0
(`data/raw/ne_10m_admin_0_countries/`). Output `data/raw/spei_country_year.csv`
covers 237 countries x 1990-2020; no indicator depends on country size
(rho with log area -0.07 to -0.01; `analysis/spei_size_check.R`). Known droughts appear as expected (Syria 2008:
SPEI-12 -1.8, 88% of area in drought; South Africa 2015; Australia 2019).
Natural Earth files go in `data/raw/ne_10m_admin_0_countries/`.
Not covered: French overseas departments (part of France in the admin-0 file)
and several atolls without land cells (Maldives, Marshall Islands, Nauru,
Tuvalu, Tokelau, Bermuda, Cook Islands, Wallis & Futuna).

Models: 1,369 refugee country-years, 1,200 migrant country-periods (vs 1,034).

* Refugee role membership: no drought association for any role (closest:
  two-way countries, OR 1.20 per SD of drought extent, 1.00-1.45, p = 0.055);
  no lagged effect.
* Refugee role changes: receiver-only countries began sending refugees more
  often in drought years (OR 1.33, 1.03-1.72, p = 0.030; rise from the previous
  year 1.27, p = 0.047). Conflict deaths remain the strongest predictor
  (OR 2.29, p = 0.005). The earlier minimum-SPEI result for countries that
  start hosting does not hold (OR 0.98, p = 0.92).
* Migrant (backbone) roles: no association with drought extent or SPEI-48.
* Drought x income: no interaction in any model.
* Results are essentially the same without the conflict variables, so drought
  is not simply acting through conflict at the country-year scale.
* Income and conflict results are unchanged.

With 30+ drought tests, p-values near 0.03-0.05 are weak evidence; the
drought -> start-sending association is the only climate result consistent
across the interim and final measures.

## Interim results (`config_drought_interim.yml`, 7 Nov 2026)

Using the size-neutral SPEI-01 annual mean as a stand-in:

* Refugees: a drier previous year goes with being a two-way country
  (OR 0.72 per SD wetter, p = 0.001) and slightly less with being sender-only
  (OR 1.25, p = 0.03); receiver-only countries are wetter in the current year
  (OR 1.23, p = 0.046). Results are nearly unchanged without conflict controls.
* Refugee role changes: receiver-only countries start sending more often in
  drier years (OR 0.71, p = 0.009).
* Migrants (backbone roles): no drought association.
* Drought x income: no interaction in any model.
* Income and conflict results are unchanged.

## References

Beguería S, Vicente-Serrano SM, Reig F, Latorre B. 2014. Standardized precipitation evapotranspiration index (SPEI) revisited. Int J Climatol. 34:3001–3023.
Black R, et al. 2011. The effect of environmental change on human migration. Glob Environ Change. 21:S3–S11.
Carling J. 2002. Migration in the age of involuntary immobility. J Ethn Migr Stud. 28(1):5–42.
Cattaneo C, Peri G. 2016. The migration response to increasing temperatures. J Dev Econ. 122:127–146.
de Haas H. 2021. A theory of migration: the aspirations-capabilities framework. Comp Migr Stud. 9:8.
Hoffmann R, Dimitrova A, Muttarak R, Crespo Cuaresma J, Peisker J. 2020. A meta-analysis of country-level studies on environmental change and migration. Nat Clim Change. 10:904–912.
Vicente-Serrano SM, Beguería S, López-Moreno JI. 2010. A multiscalar drought index sensitive to global warming: the SPEI. J Clim. 23:1696–1718.
