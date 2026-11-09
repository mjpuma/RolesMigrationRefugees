# RolesMigrationRefugees

Code for a manuscript in preparation by Zurek-Ost, Hood, Puma, Muneepeerakul,
Traff and Johnson on the roles of countries in international **refugee** and **migrant**
networks. The analysis has three steps:

1. Build a network of flows between countries for each year or period.
2. Assign each country a role with **regular equivalence (REGE)**.
3. Relate those roles, and changes of role, to conflict, income, drought and other disasters.

## Find what you need

| I want to... | Go to |
|---|---|
| Run the code | [Run the analysis](#run-the-analysis) below |
| See which file does what | [What is where](#what-is-where) below |
| Get the raw data | [docs/data.md](docs/data.md) |
| Change settings or add a dataset | [docs/configuration.md](docs/configuration.md) |
| Understand the regression model and its outputs | [docs/methods.md](docs/methods.md) |
| See how drought is measured | [data_prep/SPEI_README.md](data_prep/SPEI_README.md) |
| See the original scripts | [legacy/README.md](legacy/README.md) |

## Run the analysis

You need R 4.2 or later, the packages listed in [docs/data.md](docs/data.md), and
the raw data in `data/raw/`. Run everything from the repository root.

**Main analysis only:**

```sh
Rscript data_prep/build_spei_country.R      # once, about 20 seconds
Rscript run_pipeline.R
Rscript analysis/migrant_backbone.R
Rscript analysis/refugee_role_changes.R
```

**Everything in the manuscript and SI.** First run the main analysis and its
variants, then the regression checks and figures, in this order:

```sh
Rscript data_prep/build_spei_country.R

for cfg in config.yml configs/config_*.yml; do
  Rscript run_pipeline.R $cfg
  Rscript analysis/migrant_backbone.R $cfg
  Rscript analysis/refugee_role_changes.R $cfg
done

Rscript analysis/regression_robustness.R           # about 100 minutes on 2 cores
Rscript analysis/regression_robustness_summary.R
Rscript analysis/regression_bayes.R                # about 20 minutes
Rscript analysis/within_between_tests.R
Rscript analysis/glmmTMB_vs_lme4.R
Rscript analysis/regression_evidence.R
Rscript analysis/manuscript_figures.R              # Figs. 3 and 4 need the bootstrap results
```

Results go to `output/` (or `output_<variant>/` for a variant). These folders are not tracked.

## What is where

| Location | Contents |
|---|---|
| `config.yml` | Settings for the main analysis |
| `configs/` | Variants of the main analysis (table below) |
| `run_pipeline.R` | Flows, roles, drivers and role models for each dataset |
| `R/` | Functions used by every script: flows, REGE, roles, drivers, models |
| `analysis/` | One script per analysis (listed below) |
| `data_prep/` | Builds the drought indicators and downloads current UNHCR data |
| `data_ref/` | Country list and country-name crosswalk |
| `data/raw/` | Raw data, not tracked |
| `src/` | C++ code for valued REGE |
| `legacy/`, `for_colleagues/` | Original scripts and a corrected section for them |

Scripts in `analysis/`:

| Purpose | Script |
|---|---|
| Migrant backbone roles and their models | `migrant_backbone.R` |
| Changes of refugee role | `refugee_role_changes.R` |
| Robustness checks and country bootstrap | `regression_robustness.R`, `regression_robustness_summary.R` |
| Comparison of models | `regression_bayes.R`, `within_between_tests.R`, `glmmTMB_vs_lme4.R`, `regression_evidence.R` |
| Manuscript and SI figures | `manuscript_figures.R` |
| Supporting analyses, run separately | `migrant_volume_robustness.R`, `migrant_net_roles.R`, `spei_size_check.R`, `appendix_figures.R`, `validate_against_original.R` |

Config files:

| Config | Purpose |
|---|---|
| `config.yml` | Main analysis |
| `configs/config_blended.yml` | Earlier model with one coefficient per driver |
| `configs/config_nocontrols_conflict.yml` | Main analysis without conflict variables |
| `configs/config_drought_agriculture.yml` | Drought by agriculture share instead of income terciles |
| `configs/config_original_spei01.yml`, `configs/config_spei_jordan.yml`, `configs/config_drought_interim.yml` | Earlier drought measures, for comparison |

## License

The pipeline code is released under the [MIT License](LICENSE). The scripts in
`legacy/` remain the work of their original authors. Contact the authors before
reusing them. Data are subject to the terms of their providers.
