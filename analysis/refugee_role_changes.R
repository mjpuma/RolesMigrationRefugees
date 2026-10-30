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
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(readr); library(ggplot2)
})
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
cfg <- read_config("config.yml")
outdir <- dirname(out_path(cfg, "unhcr", "role_changes", "x"))
covs <- unlist(cfg$models$covariates)

md <- read_csv(out_path(cfg, "unhcr", "model_data.csv"), show_col_types = FALSE)
covs <- covs[covs %in% names(md) & vapply(covs, \(v) sd(md[[v]], na.rm = TRUE) > 0, logical(1))]

d <- md |>
  arrange(iso3, period) |>
  group_by(iso3) |>
  mutate(prev_role = lag(role), prev_period = lag(period),
         across(all_of(covs), \(x) x - lag(x), .names = "d_{.col}")) |>
  ungroup() |>
  filter(!is.na(prev_role), period == prev_period + 1)     # consecutive years only

changes <- list(
  `start sending` = list(at_risk = 3,      event = c(1, 2)),
  `stop sending`  = list(at_risk = c(1, 2), event = 3),
  `start hosting` = list(at_risk = 1,      event = c(2, 3)))

counts <- imap_dfr(changes, \(ch, nm) d |> filter(prev_role %in% ch$at_risk) |>
  summarise(change = nm, at_risk = n(), events = sum(role %in% ch$event)))
write_csv(counts, file.path(outdir, "transition_counts.csv"))
print(counts)

fit_one <- function(dd, vars) {
  f <- as.formula(paste("y ~", paste(c(vars, "period_c"), collapse = " + "), "+ (1 | iso3)"))
  f_glm <- as.formula(paste("y ~", paste(c(vars, "period_c"), collapse = " + ")))
  fit <- tryCatch(lme4::glmer(f, data = dd, family = binomial,
                              control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))),
                  error = function(e) NULL)
  engine <- "mixed (country random intercept)"
  co <- if (!is.null(fit)) summary(fit)$coefficients
  # a failed mixed fit shows up as implausibly tiny standard errors: refit as a
  # plain logistic regression and say so
  if (is.null(fit) || any(co[, 2] < 0.01) || length(fit@optinfo$conv$lme4$messages)) {
    fit <- glm(f_glm, data = dd, family = binomial)
    co <- summary(fit)$coefficients
    engine <- "logistic (mixed model did not converge)"
  }
  tibble(term = rownames(co), estimate = co[, 1], std.error = co[, 2], p.value = co[, 4], engine = engine)
}

res <- imap_dfr(changes, function(ch, nm) {
  dd <- d |> filter(prev_role %in% ch$at_risk) |> mutate(y = as.integer(role %in% ch$event))
  bind_rows(
    fit_one(dd, covs) |> mutate(spec = "level in year t"),
    fit_one(dd |> filter(if_all(all_of(paste0("d_", covs)), \(x) !is.na(x))), paste0("d_", covs)) |>
      mutate(spec = "change from t-1 to t")) |>
    mutate(change = nm, n = nrow(dd), events = sum(dd$y))
}) |>
  filter(!grepl("Intercept|period", term)) |>
  mutate(driver = sub("^d_", "", term), odds_ratio = exp(estimate),
         or_low = exp(estimate - 1.96 * std.error), or_high = exp(estimate + 1.96 * std.error),
         unstable = std.error > 10 | abs(estimate) > 10)
write_csv(res, file.path(outdir, "transition_models.csv"))
print(res |> filter(p.value < 0.05, !unstable) |>
        select(change, spec, driver, odds_ratio, or_low, or_high, p.value) |>
        mutate(across(where(is.double), \(x) signif(x, 3))) |> as.data.frame())

theme_set(theme_minimal(10) + theme(plot.title.position = "plot", legend.position = "bottom"))
pT <- res |> filter(!unstable) |>
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
