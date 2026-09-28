# =============================================================================
# models.R — which country characteristics go with which role?
#
# (a) One binary model per role: P(role = r) vs. all other roles, with a
#     random intercept per country (repeated observations over periods).
#     This is what the manuscript reports (odds ratios per role).
# (b) Optionally one multinomial baseline-category logit over all roles
#     (mclogit::mblogit), the model the manuscript names.
#
# Covariates are optionally log1p-transformed (heavy-tailed counts, GDP) and
# z-scored so odds ratios are "per 1 SD", comparable across drivers.
# =============================================================================

prepare_model_data <- function(cfg, roles, drivers) {
  mc <- cfg$models
  covs <- unlist(mc$covariates)
  missing_cols <- setdiff(covs, names(drivers))
  if (length(missing_cols)) stop("Covariates not in driver data: ", paste(missing_cols, collapse = ", "))

  drv <- drivers
  if ((mc$lag_periods %||% 0) > 0) {
    # use the driver values of the previous period(s) for each country
    drv <- drv |> arrange(iso3, period) |> group_by(iso3) |>
      mutate(across(all_of(covs), \(x) lag(x, mc$lag_periods))) |> ungroup()
  }

  d <- roles |>
    select(iso3, country, region, period, role, label) |>
    left_join(drv |> select(iso3, period, all_of(covs)), by = c("iso3", "period"))

  # report missingness BEFORE dropping anything (silent drops were a problem)
  miss <- d |> summarise(across(all_of(covs), \(x) sum(is.na(x)))) |>
    pivot_longer(everything(), names_to = "covariate", values_to = "n_missing") |>
    mutate(share_missing = n_missing / nrow(d))

  logc <- covs[matches_any(covs, unlist(mc$log_transform))]
  # per-capita rates are expressed per 100,000 people before log1p, so small
  # rates are not squashed to ~0
  pc <- unlist(cfg$drivers$per_capita)
  to_log <- function(x, v) log1p(pmax(x, 0) * if (matches_any(v, pc)) 1e5 else 1)
  d <- d |> mutate(across(all_of(logc), \(x) to_log(x, cur_column())))

  n_before <- nrow(d)
  d <- d |> filter(if_all(all_of(covs), \(x) !is.na(x) & is.finite(x)))
  if (isTRUE(mc$scale)) d <- d |> mutate(across(all_of(covs), \(x) as.numeric(scale(x))))
  # drop covariates with no variation (e.g. a disaster type that never occurs)
  keep <- covs[vapply(covs, function(v) sd(d[[v]]) > 0, logical(1))]

  d <- d |> mutate(period_f = factor(period), period_c = (period - mean(period)) / 10)
  list(data = d, covariates = keep, missingness = miss, n_dropped = n_before - nrow(d))
}

model_formula <- function(lhs, covs, time_effect, random = TRUE) {
  rhs <- covs
  if (time_effect == "linear") rhs <- c(rhs, "period_c")
  if (time_effect == "factor") rhs <- c(rhs, "period_f")
  f <- paste(lhs, "~", paste(rhs, collapse = " + "))
  if (random) f <- paste(f, "+ (1 | iso3)")
  as.formula(f)
}

fit_binary_role_models <- function(cfg, md) {
  mc <- cfg$models
  d <- md$data
  map_dfr(sort(unique(d$role)), function(r) {
    d$y <- as.integer(d$role == r)
    if (sum(d$y) < 10 || sum(1 - d$y) < 10) {
      msg("    role ", r, ": fewer than 10 cases -> not modelled")
      return(NULL)
    }
    if (mc$binary_engine == "mblogit") {
      fit <- mclogit::mblogit(model_formula("factor(y)", md$covariates, mc$time_effect, FALSE),
                              random = ~ 1 | iso3, data = d, method = "PQL",
                              estimator = "ML",
                          maxit = 100)
      co <- summary(fit)$coefficients
      est <- tibble(term = sub("^1~", "", rownames(co)), estimate = co[, 1],
                    std.error = co[, 2], p.value = co[, 4])
    } else {
      fit <- tryCatch(
        lme4::glmer(model_formula("y", md$covariates, mc$time_effect), data = d,
                    family = binomial,
                    control = lme4::glmerControl(optimizer = "bobyqa",
                                                 optCtrl = list(maxfun = 2e5))),
        error = function(e) {
          msg("    role ", r, ": glmer failed (", conditionMessage(e),
              ") -> plain logistic regression without random effect")
          glm(model_formula("y", md$covariates, mc$time_effect, FALSE),
              data = d, family = binomial)
        })
      co <- summary(fit)$coefficients
      est <- tibble(term = rownames(co), estimate = co[, 1],
                    std.error = co[, 2], p.value = co[, 4])
      if (inherits(fit, "merMod") && lme4::isSingular(fit))
        msg("    role ", r, ": singular fit (random-effect variance ~0)")
    }
    est |>
      mutate(role = r, n_obs = nrow(d), n_events = sum(d$y),
             odds_ratio = exp(estimate),
             or_low = exp(estimate - 1.96 * std.error),
             or_high = exp(estimate + 1.96 * std.error),
             # flag numerically broken estimates (complete separation etc.)
             unstable = std.error > 10 | abs(estimate) > 10 | std.error < 1e-4)
  })
}

fit_multinomial <- function(cfg, md) {
  d <- md$data |> mutate(role = factor(role))
  f <- model_formula("role", md$covariates, cfg$models$time_effect, FALSE)
  fit <- mclogit::mblogit(f, random = ~ 1 | iso3, data = d, method = "PQL",
                          estimator = "ML",
                          maxit = 100)
  co <- summary(fit)$coefficients
  tibble(equation = sub("~.*", "", rownames(co)), term = sub(".*~", "", rownames(co)),
         estimate = co[, 1], std.error = co[, 2], p.value = co[, 4],
         odds_ratio = exp(co[, 1]))
}

run_models <- function(cfg, dataset, roles, drivers) {
  md <- prepare_model_data(cfg, roles, drivers)
  msg("  model data: ", nrow(md$data), " country-periods (", md$n_dropped,
      " dropped for missing covariates), ", length(md$covariates), " covariates")
  write_csv(md$missingness, out_path(cfg, dataset, "covariate_missingness.csv"))
  write_csv(md$data, out_path(cfg, dataset, "model_data.csv"))

  bin <- fit_binary_role_models(cfg, md)
  write_csv(bin, out_path(cfg, dataset, "binary_role_models.csv"))
  res <- list(binary = bin, data = md)
  if (isTRUE(cfg$models$multinomial)) {
    mn <- tryCatch(fit_multinomial(cfg, md), error = function(e) {
      msg("  multinomial model failed: ", conditionMessage(e)); NULL })
    if (!is.null(mn)) write_csv(mn, out_path(cfg, dataset, "multinomial_role_model.csv"))
    res$multinomial <- mn
  }
  res
}
