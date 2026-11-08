# =============================================================================
# refugee_role_changes.R — what precedes a change of refugee role?
#
#   Rscript analysis/refugee_role_changes.R     (~1 min; after run_pipeline.R)
#
# Refugee roles change often (27% of country-years). Instead of asking which
# drivers go with a role, this asks which drivers go with a CHANGE of role
# between consecutive years, among countries at risk of that change:
#
#   start sending   receiver only (R3) in year t-1  ->  sender only or two-way in t
#   stop sending    sender only or two-way in t-1   ->  receiver only in t
#   start hosting   sender only (R1) in t-1         ->  two-way or receiver only in t
#
# Drivers: the same standardised country-year drivers as the role models
# (output/unhcr/model_data.csv), measured (a) in year t, the year of the change,
# and (b) as the change from t-1 to t (a "shock" relative to the country's
# previous year). Mixed-effects logit with a country random intercept; years
# without recorded refugee flows break a spell (no transition is counted).
#
# Writes to output/unhcr/role_changes/:
#   transition_counts.csv, transition_models.csv, T1_transition_odds_ratios.png
#
# Model, written out for one type of change. Among country-years at risk
# (country c was in an at-risk role in year t-1 and is observed in year t):
#
#   logit P(z_ct = 1) = a + sum_j b_j x_jct + h (t - tbar)/10 + u_c
#
#   z_ct = 1 if the change happens between t-1 and t
#   x_jct  drivers in year t ("level" spec) or x_jct - x_jc,t-1 ("change" spec)
#   linear time trend in decades (period_c); annual fixed effects would use
#   up too many degrees of freedom with ~70 events
#   u_c ~ N(0, s_u^2) country random intercept
# The level spec uses the same covariates as the refugee role models,
# including drought slopes by income tercile (Wald test of equal slopes).
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(readr); library(ggplot2)
})
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
cfg_file <- commandArgs(trailingOnly = TRUE)   # optional: Rscript <script> my_config.yml
cfg <- read_config(if (length(cfg_file)) cfg_file[1] else "config.yml")
outdir <- dirname(out_path(cfg, "unhcr", "role_changes", "x"))
# the standardized model data of the refugee role models (one row per country-year)
md <- read_csv(out_path(cfg, "unhcr", "model_data.csv"), show_col_types = FALSE)

# Same specification as the refugee role models: read the covariates that
# run_pipeline.R actually used (incl. lags, income-group dummies and
# group-specific drought slopes). Older outputs without that file fall back to
# the covariates listed in the config.
cov_file <- out_path(cfg, "unhcr", "model_covariates.txt")
covs_level <- if (file.exists(cov_file)) readLines(cov_file) else unlist(cfg$models$covariates)
covs_level <- covs_level[covs_level %in% names(md)]
slope_file <- out_path(cfg, "unhcr", "model_slope_groups.csv")
# named list: drought variable -> its tercile slope columns (for the Wald test)
slope_groups <- if (file.exists(slope_file)) split(read_csv(slope_file, show_col_types = FALSE)$column,
                                                   read_csv(slope_file, show_col_types = FALSE)$variable) else list()
# Base drivers only (no lags, products or group dummies) for the change-from-
# previous-year specification: a year-on-year difference of a product is not meaningful.
covs <- covs_level[!grepl("_lag1|_x_|^income_g|_base$", covs_level)]
# a drought variable with tercile-specific slopes is not itself in covs_level
# (its slope columns replace it); add it back so its year-on-year change enters
covs <- unique(c(covs, setdiff(names(slope_groups), grep("_lag1$", names(slope_groups), value = TRUE))))
covs <- covs[covs %in% names(md)]
# drop drivers with no variation in the model data
covs <- covs[vapply(covs, \(v) sd(md[[v]], na.rm = TRUE) > 0, logical(1))]

# previous year's role and the year-on-year change of each base driver
# (d_<var> = value in t minus value in t-1, in SD units of the level)
d <- md |>
  arrange(iso3, period) |>
  group_by(iso3) |>
  mutate(prev_role = lag(role), prev_period = lag(period),
         across(all_of(covs), \(x) x - lag(x), .names = "d_{.col}")) |>
  ungroup() |>
  filter(!is.na(prev_role), period == prev_period + 1)     # consecutive years only

# the three changes: roles at risk in t-1 and roles that count as the event in t
# (1 = sender only, 2 = two-way, 3 = receiver only)
changes <- list(
  `start sending` = list(at_risk = 3,      event = c(1, 2)),
  `stop sending`  = list(at_risk = c(1, 2), event = 3),
  `start hosting` = list(at_risk = 1,      event = c(2, 3)))

# number of country-years at risk and number of events for each change
counts <- imap_dfr(changes, \(ch, nm) d |> filter(prev_role %in% ch$at_risk) |>
  summarise(change = nm, at_risk = n(), events = sum(role %in% ch$event)))
write_csv(counts, file.path(outdir, "transition_counts.csv"))
print(counts)

# Fit one role-change model.
#   dd      at-risk country-years with outcome y
#   vars    covariates
#   slopes  tercile slope groups to test for equality (level spec only)
# Returns one row per term (log-odds estimate, SE, p, engine) plus Wald rows.
fit_one <- function(dd, vars, slopes = list()) {
  # drop terms that do not vary among the countries at risk (e.g. an income group with no events)
  vars <- vars[vapply(vars, \(v) sd(dd[[v]], na.rm = TRUE) > 0, logical(1))]
  f <- as.formula(paste("y ~", paste(c(vars, "period_c"), collapse = " + "), "+ (1 | iso3)"))
  f_glm <- as.formula(paste("y ~", paste(c(vars, "period_c"), collapse = " + ")))
  # same robust fitting as the role models (R/models.R): several optimizers,
  # then ordinary logistic regression if no acceptable mixed fit is found
  fit <- fit_glmer_robust(f, f_glm, dd, label = "role change")
  engine <- attr(fit, "engine")
  co <- summary(fit)$coefficients
  out <- tibble(term = rownames(co), estimate = co[, 1], std.error = co[, 2], p.value = co[, 4], engine = engine)
  # equality of the group-specific drought slopes (Wald chi-square; df in std.error)
  b <- if (inherits(fit, "merMod")) lme4::fixef(fit) else coef(fit)
  V <- as.matrix(vcov(fit))
  for (v in names(slopes)) {
    w <- wald_equal(b, V, slopes[[v]])
    out <- bind_rows(out, tibble(term = paste0("wald_equal_slopes:", v), estimate = w[["chi2"]],
                                 std.error = w[["df"]], p.value = w[["p"]], engine = engine))
  }
  out
}

# both specifications for each change; the change spec uses only rows where
# every d_ variable is available
res <- imap_dfr(changes, function(ch, nm) {
  dd <- d |> filter(prev_role %in% ch$at_risk) |> mutate(y = as.integer(role %in% ch$event))
  bind_rows(
    fit_one(dd, covs_level, slope_groups) |> mutate(spec = "level in year t"),
    fit_one(dd |> filter(if_all(all_of(paste0("d_", covs)), \(x) !is.na(x))), paste0("d_", covs)) |>
      mutate(spec = "change from t-1 to t")) |>
    mutate(change = nm, n = nrow(dd), events = sum(dd$y))
}) |>
  # drop intercept and trend; odds ratios with 95% CIs; flag unstable fits
  filter(!grepl("Intercept|period", term)) |>
  mutate(driver = sub("^d_", "", term), odds_ratio = exp(estimate),
         or_low = exp(estimate - 1.96 * std.error), or_high = exp(estimate + 1.96 * std.error),
         unstable = std.error > 10 | abs(estimate) > 10,
         # Wald rows hold chi-square and df, not log-odds: no odds ratio
         across(c(odds_ratio, or_low, or_high), \(x) ifelse(grepl("^wald_", term), NA, x)),
         unstable = ifelse(grepl("^wald_", term), FALSE, unstable))
write_csv(res, file.path(outdir, "transition_models.csv"))
print(res |> filter(p.value < 0.05, !unstable, !grepl("^wald_", term)) |>
        select(change, spec, driver, odds_ratio, or_low, or_high, p.value) |>
        mutate(across(where(is.double), \(x) signif(x, 3))) |> as.data.frame())

# T1: odds ratios by change (columns) and specification (rows)
theme_set(theme_minimal(10) + theme(plot.title.position = "plot", legend.position = "bottom"))
pT <- res |> filter(!unstable, !grepl("^wald_", term)) |>
  mutate(sig = p.value < 0.05,
         change = factor(paste0(change, " (", events, " of ", n, ")"),
                         unique(paste0(change, " (", events, " of ", n, ")"))),
         spec = factor(spec, c("level in year t", "change from t-1 to t"))) |>
  ggplot(aes(odds_ratio, driver, colour = sig)) +
  geom_vline(xintercept = 1, colour = "grey60") +
  geom_errorbarh(aes(xmin = or_low, xmax = or_high), height = 0) + geom_point() +
  scale_x_log10() + coord_cartesian(xlim = c(0.2, 5)) + facet_grid(spec ~ change) +
  scale_colour_manual(values = c(`FALSE` = "grey55", `TRUE` = "#b34b4b"), labels = c("p >= 0.05", "p < 0.05")) +
  labs(title = "Refugees: what goes with a change of role from one year to the next?",
       subtitle = "Odds ratios per 1 SD (95% CI). Columns: type of change (events of country-years at risk).\nRows: driver level in the year of the change, or its change from the previous year.",
       x = "Odds ratio (log scale)", y = NULL, colour = NULL)
ggsave(file.path(outdir, "T1_transition_odds_ratios.png"), pT, width = 11, height = 9, dpi = 200, bg = "white")
msg("Done: ", outdir)
