# Methods notes

The regression model, its outputs and notes on the role analysis. Back to the [README](../README.md).

## Regression model and uncertainty

The role and role-change models are mixed-effects logistic regressions with a
country random intercept, fitted by maximum likelihood. GDP per capita, conflict
events and conflict deaths are each separated into a between-country term (the
country's mean) and a within-country term (the deviation from it), following
Bafumi and Gelman (2006). `analysis/within_between_tests.R` tests whether the two
coefficients are equal, which a model with one coefficient per driver assumes.

The confidence intervals in `binary_role_models.csv` and `transition_models.csv`
(`or_low`, `or_high`) are model-based. Given its intercept, the model treats a
country's successive periods as independent, but roles and drivers persist. Under
the placebo test that gives each country another country's full driver history,
within-country terms of the refugee role models exceed |t| = 1.96 in about 20% of
runs instead of 5%. Figs. 3 and 4 and the Results therefore use 95% percentile
intervals from a country cluster bootstrap (`output/robustness/bootstrap.csv`:
whole countries resampled with replacement, 500 samples, non-converged refits
discarded). Fig. 5 and SI Tables S1-S4 remain model-based. For rare EM-DAT disaster
terms, a few events drive the estimate and the percentile interval can be too
narrow or lopsided, so those terms should be read together with the model-based
result and the partially pooled estimate.

Outputs in `output/robustness/` (not tracked):

| File | Content |
|---|---|
| `baseline.csv` | main-model estimates that the checks are compared with |
| `sample_refits.csv`, `spec_refits.csv`, `data_refits.csv` | every refit: subsamples, specifications, alternative data |
| `inference.csv` | pooled-logit standard errors: clustered, two-way, Conley |
| `bootstrap.csv` | country cluster bootstrap of the mixed models |
| `placebo_draws.csv`, `placebo_summary.csv` | placebo runs and randomization p-values |
| `summary_table.csv` | one row per model and driver, all checks combined |
| `within_between_tests.csv`, `within_between_estimates.csv` | tests of equal within and between coefficients |
| `glmmTMB_vs_lme4.csv`, `lme4_restart_check.csv` | comparison of the two fitting programs |
| `bayes/` | Bayesian coefficients, convergence diagnostics, posterior predictive checks |
| `partial_pooling.csv`, `family_tau.csv` | partial pooling by family of drivers |
| `evidence_table.csv`, `Fig_evidence.png` | all estimates for each association (Table S9, Fig. S12) |

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

## Pipeline steps and their outputs

Results go to `output/<dataset>/` (not tracked):

| Step | Code | Output |
|---|---|---|
| 1. Read flows | `R/flows.R` | – |
| 2. Stacked country-period network → REGE → hierarchical clustering | `R/rege.R`, `R/roles.R` | `roles_by_country_period.csv`, `role_profile.csv`, `k_diagnostics.csv`, `role_transitions.csv`, `role_stability_by_country.csv` |
| 3. Driver panel aggregated to each dataset's periods | `R/drivers.R` | `drivers_by_country_period.csv` |
| 4. Role models | `R/models.R` | `binary_role_models.csv`, `multinomial_role_model.csv`, `covariate_missingness.csv` |
| 5. Figures | `R/figures.R` | `figures/*.png` |

To compare with the original role assignments, run
`Rscript analysis/validate_against_original.R`. It reads the original outputs in
`legacy/results/` (`UNDESA_full_roles_k3.csv`, `UNHCR_full_roles_k3.csv`;
role labels only, as produced by the legacy scripts).
