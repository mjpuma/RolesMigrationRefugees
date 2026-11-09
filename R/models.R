# =============================================================================
# models.R — which country characteristics go with which role?
#
# (a) One binary model per role: P(role = r) vs. all other roles, with a
#     random intercept per country (repeated observations over periods).
#     This is what the manuscript reports (odds ratios per role). The main
#     config fits it in glmmTMB (binary_engine: glmmTMB; see fit_mixed()).
# (b) Optionally one multinomial baseline-category logit over all roles
#     (mclogit::mblogit), the model the manuscript names.
#
# Covariates are optionally log1p-transformed (heavy-tailed counts, GDP) and
# z-scored so odds ratios are "per 1 SD", comparable across drivers.
#
# Model (a), written out for country c in period t and role r:
#
#   logit P(y_ct = 1) = a + sum_j b_j x_jct + g_t + u_c,   u_c ~ N(0, s_u^2)
#
#   y_ct = 1 if country c is in role r in period t, else 0
#   x_jct  standardized covariate j (mean 0, SD 1 over the model sample)
#   g_t    period fixed effects (time_effect: factor) or a linear trend
#   u_c    country random intercept (absorbs stable country differences)
#   exp(b_j) = odds ratio for a 1-SD higher x_j, other covariates fixed
#
# Within and between countries (config: within_between; Bafumi and Gelman 2006).
# A driver x can differ across countries and over time in one country, and the
# two need not have the same association with a role. Each listed driver is
# replaced by two terms, in units of the driver's overall SD:
#
#   x_ct = xbar_c + (x_ct - xbar_c)
#   ... + b_B xbar_c + b_W (x_ct - xbar_c) + ...
#
#   b_W    "in a year or period when country c has more of x than is usual for it"
#   b_B    "a country that has more of x than other countries, on average"
# With one coefficient, the model assumes b_W = b_B (analysis/within_between_tests.R).
#
# Standard errors. u_c captures stable differences between countries, but given
# u_c the model treats a country's successive periods as independent. Roles and
# drivers persist, so the model-based standard errors of the within-country terms
# are too small. Figs. 3 and 4 and the Results therefore use a country cluster
# bootstrap (analysis/regression_robustness.R, part "bootstrap"); the or_low and
# or_high columns written here remain model-based (exp(b +/- 1.96 SE)).
#
# Drought by income group (config: income_groups, slopes_by_group). The
# common drought slope b_D is replaced by one slope per baseline-income
# tercile k = 1, 2, 3:
#
#   ... + sum_k d_k D_ct 1[G_c = k] + e_2 1[G_c = 2] + e_3 1[G_c = 3] + ...
#
#   D_ct   standardized drought indicator
#   G_c    tercile of the country's GDP per capita in its first period
#          (fixed over time, so it cannot respond to the drought it moderates)
#   d_k    drought log-odds slope for countries in tercile k
#   e_k    tercile main effects (tercile 1, the poorest, is the reference)
# H0: d_1 = d_2 = d_3 is tested with a 2-df Wald test (wald_equal()).
#
# Interactions (config: interactions), e.g. drought x baseline agriculture
# share: the product of two standardized variables is added as a covariate.
# Its coefficient is the change in the drought log-odds slope per 1 SD of
# the moderator; because both inputs are z-scores the product is unitless.
# =============================================================================

# Assemble the model data for one dataset: join roles and period-level drivers,
# transform and standardize covariates, and build the derived terms
# (baselines, lags, income groups, group-specific slopes, interactions).
#   roles    roles_by_country_period table from run_roles()
#   drivers  output of aggregate_drivers() for this dataset's periods
# Returns list(data, covariates, missingness, n_dropped, slope_cols).
prepare_model_data <- function(cfg, roles, drivers) {
  mc <- cfg$models
  covs <- unlist(mc$covariates)
  drv <- drivers

  # Baseline copies of selected drivers: each country's value in its earliest
  # period with data, named <var>_base (time-invariant within a country).
  # Main config: gdp_pc_base (used to form income terciles). Agriculture
  # config: agri_share_base. The baseline is taken from the full driver table,
  # BEFORE rows with missing covariates are dropped.
  for (v in unlist(mc$baseline)) {
    drv <- drv |> arrange(iso3, period) |> group_by(iso3) |>
      mutate("{v}_base" := dplyr::first(.data[[v]][!is.na(.data[[v]])])) |> ungroup()
  }
  # stop early with a clear message if the config names a covariate that the
  # driver data does not have (e.g. a typo, or a source file not loaded)
  missing_cols <- setdiff(covs, names(drv))
  if (length(missing_cols)) stop("Covariates not in driver data: ", paste(missing_cols, collapse = ", "))
  if ((mc$lag_periods %||% 0) > 0) {
    # use the driver values of the previous period(s) for each country
    # (option kept from the original scripts; 0 in all current configs)
    drv <- drv |> arrange(iso3, period) |> group_by(iso3) |>
      mutate(across(all_of(covs), \(x) lag(x, mc$lag_periods))) |> ungroup()
  }

  # Lagged copies of selected drivers: value in the previous period of the same
  # country (consecutive periods only), named <var>_lag1.
  # Refugee data (annual): drought12_ext_lag1 = drought extent in year t-1.
  # `step` is the spacing between periods (1 for annual, 5 for five-year
  # data); if the previous observed period is further back (a gap), the lag
  # is set to NA rather than borrowing an older value.
  lag_vars <- unlist(mc$lags)
  if (length(lag_vars)) {
    step <- min(diff(sort(unique(drv$period))))
    drv <- drv |> arrange(iso3, period) |> group_by(iso3) |>
      mutate(across(all_of(lag_vars), \(x) ifelse(period - lag(period) == step, lag(x), NA),
                    .names = "{.col}_lag1")) |> ungroup()
    covs <- c(covs, paste0(lag_vars, "_lag1"))
  }

  grp <- mc$income_groups                 # e.g. list(var = "gdp_pc_base", n = 3)
  # the grouping variable is carried along (and must be non-missing) even
  # though it is not itself a covariate
  extra <- if (!is.null(grp)) grp$var else character(0)
  d <- roles |>
    select(iso3, country, region, period, role, label) |>
    left_join(drv |> select(iso3, period, all_of(unique(c(covs, extra)))), by = c("iso3", "period"))

  # report missingness BEFORE dropping anything (silent drops were a problem)
  miss <- d |> summarise(across(all_of(covs), \(x) sum(is.na(x)))) |>
    pivot_longer(everything(), names_to = "covariate", values_to = "n_missing") |>
    mutate(share_missing = n_missing / nrow(d))

  # log1p transform for skewed variables named in cfg$models$log_transform
  # (GDP, conflict counts, deaths, affected). A lagged variable follows the
  # rule of its base variable (the "_lag1" suffix is ignored when matching).
  logc <- covs[matches_any(sub("_lag1$", "", covs), unlist(mc$log_transform))]
  # per-capita rates are expressed per 100,000 people before log1p, so small
  # rates are not squashed to ~0
  pc <- unlist(cfg$drivers$per_capita)
  # pmax(x, 0) guards against tiny negative values before the log
  to_log <- function(x, v) log1p(pmax(x, 0) * if (matches_any(v, pc)) 1e5 else 1)
  d <- d |> mutate(across(all_of(logc), \(x) to_log(x, cur_column())))

  # complete-case sample: drop country-periods with any missing (or infinite)
  # covariate; the number dropped is reported in the log
  n_before <- nrow(d)
  d <- d |> filter(if_all(all_of(unique(c(covs, extra))), \(x) !is.na(x) & is.finite(x)))
  # z-score every covariate over the model sample: (x - mean) / SD.
  # The grouping variable (gdp_pc_base) is not scaled; it is only used to cut terciles.
  if (isTRUE(mc$scale)) d <- d |> mutate(across(all_of(covs), \(x) as.numeric(scale(x))))

  # Within/between split (see header): country mean and deviation from it,
  # computed from the standardized values over the country's rows in the model sample
  wb_vars <- intersect(unlist(mc$within_between), covs)
  if (length(wb_vars)) {
    d <- d |> group_by(iso3) |>
      mutate(across(all_of(wb_vars), list(between = mean, within = \(x) x - mean(x)), .names = "{.col}_{.fn}")) |> ungroup()
    covs <- c(setdiff(covs, wb_vars), paste0(wb_vars, "_within"), paste0(wb_vars, "_between"))
  }

  # Income groups: quantile groups (terciles for n = 3) of a baseline variable,
  # cut across COUNTRIES (one value per country), not across country-periods.
  # Cutting across countries gives each tercile about the same number of
  # countries; cutting across rows would weight countries by how many periods
  # they are observed.
  slope_cols <- list()
  if (!is.null(grp)) {
    n_g <- grp$n %||% 3
    base_by_country <- d |> distinct(iso3, .data[[grp$var]]) |> pull(.data[[grp$var]])
    # cut points: the 1/3 and 2/3 quantiles of baseline GDP across countries
    cuts <- quantile(base_by_country, probs = seq_len(n_g - 1) / n_g, names = FALSE)
    # findInterval returns 0 below the first cut, 1 between cuts, 2 above;
    # + 1 gives groups 1..n_g
    d$income_group <- findInterval(d[[grp$var]], cuts) + 1L   # 1 = lowest group
    # group main effects (group 1 is the reference)
    for (k in 2:n_g) { nm <- paste0("income_g", k); d[[nm]] <- as.numeric(d$income_group == k); covs <- c(covs, nm) }
    # separate slopes: <var>_x_g<k> = var * 1(group == k); replaces the common slope.
    # Each row has a non-zero value in exactly one of the n_g columns, so the
    # coefficient of <var>_x_g<k> is the slope of <var> within group k. (This
    # is the same model as "common slope + slope differences", reparametrised
    # so that each group's slope and its CI are read off directly.)
    for (v in unlist(mc$slopes_by_group)) {
      cols <- paste0(v, "_x_g", seq_len(n_g))
      for (k in seq_len(n_g)) d[[cols[k]]] <- d[[v]] * as.numeric(d$income_group == k)
      covs <- c(setdiff(covs, v), cols)
      slope_cols[[v]] <- cols
    }
    attr(slope_cols, "cuts") <- cuts
  }
  # Interactions as products of the (scaled) covariates, named <a>_x_<b>.
  # With scaled inputs, the main effects are effects at the mean of the other.
  for (pr in mc$interactions %||% list()) {
    pr <- unlist(pr); nm <- paste(pr, collapse = "_x_")
    d[[nm]] <- d[[pr[1]]] * d[[pr[2]]]
    covs <- c(covs, nm)
  }
  # drop covariates with no variation (e.g. a disaster type that never occurs)
  keep <- covs[vapply(covs, function(v) sd(d[[v]]) > 0, logical(1))]

  # time terms: period_f = one fixed effect per period; period_c = centred
  # period in decades (used for a linear trend and in the role-change models)
  d <- d |> mutate(period_f = factor(period), period_c = (period - mean(period)) / 10)
  list(data = d, covariates = keep, missingness = miss, n_dropped = n_before - nrow(d),
       slope_cols = slope_cols)
}

# Build a model formula, e.g. y ~ gdp_pc + ... + period_f + (1 | iso3)
#   lhs          response ("y" or "role")
#   covs         covariate names
#   time_effect  "none", "linear" (period_c) or "factor" (period_f)
#   random       add a country random intercept (1 | iso3)?
model_formula <- function(lhs, covs, time_effect, random = TRUE) {
  rhs <- covs
  if (time_effect == "linear") rhs <- c(rhs, "period_c")
  if (time_effect == "factor") rhs <- c(rhs, "period_f")
  f <- paste(lhs, "~", paste(rhs, collapse = " + "))
  if (random) f <- paste(f, "+ (1 | iso3)")
  as.formula(f)
}

# Fit a mixed logistic model, trying several optimizers (bobyqa, Nelder-Mead,
# nloptwrap; Laplace approximation, up to 200,000 function evaluations).
#   f_mixed  formula with (1 | iso3);  f_glm  the same without it
# Acceptance, in order:
#   1. "clean": lme4 reports no convergence problem and every fixed-effect
#      standard error is plausible (between 0.001 and 10 on the log-odds
#      scale; a degenerate Hessian shows up as near-zero standard errors).
#      The first optimizer giving a clean fit is used.
#   2. "optimizers agree": no clean fit, but at least two optimizers give
#      fits with plausible standard errors, no warning about a singular or
#      degenerate Hessian, and practically the same estimates: every fixed
#      effect differs between optimizers by less than 0.1 of its standard
#      error, and the random-effect SD by less than 0.05. lme4's gradient
#      check is known to flag such fits although they have converged
#      (lme4 ?convergence; Bolker's allFit advice); the fit with the smallest
#      scaled gradient is used.
#   3. otherwise: ordinary logistic regression without the country random
#      intercept, with a message.
# The rule used is stored in attr(fit, "engine") and reported in the results.
# Warnings are suppressed during fitting because the checks below read the
# convergence messages directly.
fit_glmer_robust <- function(f_mixed, f_glm, d, label = "") {
  se_ok <- function(fit) {
    se <- sqrt(diag(as.matrix(vcov(fit))))
    all(is.finite(se)) && all(se > 1e-3 & se < 10)
  }
  msgs <- function(fit) fit@optinfo$conv$lme4$messages
  # scaled gradient (lme4 ?convergence): solve(Hessian, gradient); Inf if unavailable
  relgrad <- function(fit) tryCatch(max(abs(with(fit@optinfo$derivs, solve(Hessian, gradient)))),
                                    error = function(e) Inf)
  fits <- list()
  for (opt in c("bobyqa", "Nelder_Mead", "nloptwrap")) {
    fit <- tryCatch(suppressWarnings(
      lme4::glmer(f_mixed, data = d, family = binomial,
                  control = lme4::glmerControl(optimizer = opt, optCtrl = list(maxfun = 2e5)))),
      error = function(e) NULL)
    if (is.null(fit)) next                                    # glmer stopped with an error
    if (!length(msgs(fit)) && se_ok(fit)) {                   # rule 1: clean fit
      attr(fit, "engine") <- paste0("mixed (", opt, ")"); return(fit)
    }
    # keep candidates for rule 2: plausible SEs and no singular/degenerate Hessian
    if (se_ok(fit) && !any(grepl("degenerate|singular|unidentifiable", msgs(fit)))) fits[[opt]] <- fit
  }
  if (length(fits) >= 2) {
    # spread of each estimate across the candidate fits: fixed effects in units
    # of their standard error (largest SE among the fits), theta on its own scale
    b  <- do.call(rbind, lapply(fits, lme4::fixef))
    se <- do.call(rbind, lapply(fits, function(f) sqrt(diag(as.matrix(vcov(f))))))
    th <- vapply(fits, function(f) lme4::getME(f, "theta")[[1]], numeric(1))
    spread <- max(apply(b, 2, function(x) diff(range(x))) / apply(se, 2, max))
    if (spread < 0.1 && diff(range(th)) < 0.05) {             # rule 2: optimizers agree
      best <- names(fits)[which.min(vapply(fits, relgrad, numeric(1)))]
      fit <- fits[[best]]
      msg("    ", label, ": gradient warning, but ", length(fits), " optimizers agree (largest difference ",
          signif(spread, 2), " SE) -> mixed model (", best, ")")
      attr(fit, "engine") <- paste0("mixed (", best, "; optimizers agree)")
      return(fit)
    }
  }
  msg("    ", label, ": no acceptable mixed-model fit -> logistic regression without random intercept")
  fit <- glm(f_glm, data = d, family = binomial)
  attr(fit, "engine") <- "logistic (mixed model did not converge)"
  fit
}

# Main fitting function: glmmTMB (same Laplace likelihood as lme4, much faster and
# more stable for these models); a fit that does not converge goes through
# fit_glmer_robust(). The engine used is stored in attr(fit, "engine").
fit_mixed <- function(f_mixed, f_glm, d, label = "") {
  fit <- tryCatch(suppressWarnings(suppressMessages(
    glmmTMB::glmmTMB(f_mixed, data = d, family = binomial, control = glmmTMB::glmmTMBControl(parallel = 1)))),
    error = function(e) NULL)
  if (!is.null(fit) && fit$fit$convergence == 0 && isTRUE(fit$sdr$pdHess)) {
    attr(fit, "engine") <- "mixed (glmmTMB)"
    return(fit)
  }
  fit_glmer_robust(f_mixed, f_glm, d, label = label)
}

# Coefficient table, estimates and covariance for any of the fit classes used here
co_table <- function(fit) if (inherits(fit, "glmmTMB")) summary(fit)$coefficients$cond else summary(fit)$coefficients
b_vec <- function(fit) if (inherits(fit, "glmmTMB")) glmmTMB::fixef(fit)$cond else if (inherits(fit, "merMod")) lme4::fixef(fit) else coef(fit)
v_mat <- function(fit) if (inherits(fit, "glmmTMB")) as.matrix(vcov(fit)$cond) else as.matrix(vcov(fit))

# Wald test that several coefficients are equal (e.g. the drought slope in each
# income group): H0: b_1 = b_2 = ... = b_G, tested with G-1 contrasts b_k - b_1.
#   b     named vector of fixed-effect estimates
#   V     their covariance matrix
#   cols  names of the coefficients to compare
# With contrast matrix L (one row per contrast: -1 at b_1, +1 at b_k),
#   chi2 = (L b)' (L V L')^-1 (L b),  df = G - 1,
# and p = P(chi-square_df > chi2). Returns c(chi2, df, p); NA if fewer than two
# of the coefficients are in the model or L V L' cannot be inverted.
wald_equal <- function(b, V, cols) {
  cols <- intersect(cols, names(b))
  if (length(cols) < 2) return(c(chi2 = NA_real_, df = NA_real_, p = NA_real_))
  L <- matrix(0, length(cols) - 1, length(b), dimnames = list(NULL, names(b)))
  for (k in 2:length(cols)) { L[k - 1, cols[1]] <- -1; L[k - 1, cols[k]] <- 1 }
  Lb <- L %*% b
  chi2 <- tryCatch(as.numeric(t(Lb) %*% solve(L %*% V %*% t(L)) %*% Lb), error = function(e) NA_real_)
  c(chi2 = chi2, df = nrow(L), p = pchisq(chi2, nrow(L), lower.tail = FALSE))
}

# One binary model per role (role r vs. all other roles).
# Output: one row per role and term with the log-odds estimate, SE, p-value,
# odds ratio and 95% CI (exp(b +/- 1.96 SE)), the engine used, and extra rows
# "wald_equal_slopes:<var>" for each group-specific slope (estimate = chi2,
# std.error = df, p.value = Wald p; no odds ratio). Roles with fewer than 10
# members (or fewer than 10 non-members) are not modelled.
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
      # option reproducing the original scripts: mclogit with PQL estimation
      fit <- mclogit::mblogit(model_formula("factor(y)", md$covariates, mc$time_effect, FALSE),
                              random = ~ 1 | iso3, data = d, method = "PQL",
                              estimator = "ML",
                          maxit = 100)
      co <- summary(fit)$coefficients
      # mblogit prefixes terms with the equation ("1~gdp_pc"); strip it
      est <- tibble(term = sub("^1~", "", rownames(co)), estimate = co[, 1],
                    std.error = co[, 2], p.value = co[, 4])
    } else {
      # main analysis: glmmTMB (fit_mixed), or lme4 through fit_glmer_robust() when
      # the config says binary_engine: glmer
      fit_fun <- if (identical(mc$binary_engine, "glmmTMB")) fit_mixed else fit_glmer_robust
      fit <- fit_fun(model_formula("y", md$covariates, mc$time_effect),
                     model_formula("y", md$covariates, mc$time_effect, FALSE), d,
                     label = paste("role", r))
      engine <- attr(fit, "engine")
      co <- co_table(fit)
      # columns of co: Estimate, Std. Error, z value, Pr(>|z|)
      est <- tibble(term = rownames(co), estimate = co[, 1],
                    std.error = co[, 2], p.value = co[, 4], engine = engine)
      if (inherits(fit, "merMod") && lme4::isSingular(fit))
        msg("    role ", r, ": singular fit (random-effect variance ~0)")
      # equality of group-specific slopes (one test per variable)
      b <- b_vec(fit)
      V <- v_mat(fit)
      for (v in names(md$slope_cols)) {
        w <- wald_equal(b, V, md$slope_cols[[v]])
        est <- bind_rows(est, tibble(term = paste0("wald_equal_slopes:", v), estimate = w[["chi2"]],
                                     std.error = w[["df"]], p.value = w[["p"]], engine = engine))
      }
    }
    est |>
      mutate(role = r, n_obs = nrow(d), n_events = sum(d$y),
             odds_ratio = exp(estimate),
             or_low = exp(estimate - 1.96 * std.error),
             or_high = exp(estimate + 1.96 * std.error),
             # flag numerically broken estimates (complete separation etc.)
             unstable = std.error > 10 | abs(estimate) > 10 | std.error < 1e-4) |>
      # Wald rows carry chi2 in `estimate` and df in `std.error`: no odds ratio
      mutate(across(c(odds_ratio, or_low, or_high), \(x) ifelse(grepl("^wald_", term), NA, x)),
             unstable = ifelse(grepl("^wald_", term), FALSE, unstable))
  })
}

# Multinomial baseline-category logit over all roles (role 1 = baseline),
# with a country random intercept, fitted by PQL (mclogit::mblogit).
# Optional (models: multinomial: true); reported in the SI only.
# Output: one row per equation (role vs. role 1) and term.
fit_multinomial <- function(cfg, md) {
  d <- md$data |> mutate(role = factor(role))
  f <- model_formula("role", md$covariates, cfg$models$time_effect, FALSE)
  fit <- mclogit::mblogit(f, random = ~ 1 | iso3, data = d, method = "PQL",
                          estimator = "ML",
                          maxit = 100)
  co <- summary(fit)$coefficients
  # row names look like "2~gdp_pc": equation before "~", term after
  tibble(equation = sub("~.*", "", rownames(co)), term = sub(".*~", "", rownames(co)),
         estimate = co[, 1], std.error = co[, 2], p.value = co[, 4],
         odds_ratio = exp(co[, 1]))
}

# Model step for one dataset. Writes to <output_dir>/<dataset>/:
#   covariate_missingness.csv   missing values per covariate before dropping
#   model_data.csv              the exact data the models were fitted to
#   model_covariates.txt        the final covariate list (one per line)
#   model_slope_groups.csv      which columns hold each group-specific slope
#   binary_role_models.csv      results of fit_binary_role_models()
#   multinomial_role_model.csv  (optional) results of fit_multinomial()
run_models <- function(cfg, dataset, roles, drivers) {
  # a dataset may override model settings (e.g. lags for annual refugee data)
  # (top-level replacement, so a dataset list such as interactions replaces the default)
  ov <- cfg$flows[[dataset]]$models %||% list()
  cfg$models[names(ov)] <- ov
  md <- prepare_model_data(cfg, roles, drivers)
  msg("  model data: ", nrow(md$data), " country-periods (", md$n_dropped,
      " dropped for missing covariates), ", length(md$covariates), " covariates")
  write_csv(md$missingness, out_path(cfg, dataset, "covariate_missingness.csv"))
  write_csv(md$data, out_path(cfg, dataset, "model_data.csv"))
  # the exact covariate set and slope groups, so later scripts (role changes) use the same specification
  writeLines(md$covariates, out_path(cfg, dataset, "model_covariates.txt"))
  if (length(md$slope_cols))
    write_csv(tibble(variable = rep(names(md$slope_cols), lengths(md$slope_cols)),
                     column = unlist(md$slope_cols)), out_path(cfg, dataset, "model_slope_groups.csv"))

  bin <- fit_binary_role_models(cfg, md)
  write_csv(bin, out_path(cfg, dataset, "binary_role_models.csv"))
  res <- list(binary = bin, data = md)
  if (isTRUE(cfg$models$multinomial)) {
    # a failure of the optional multinomial model is reported, not fatal
    mn <- tryCatch(fit_multinomial(cfg, md), error = function(e) {
      msg("  multinomial model failed: ", conditionMessage(e)); NULL })
    if (!is.null(mn)) write_csv(mn, out_path(cfg, dataset, "multinomial_role_model.csv"))
    res$multinomial <- mn
  }
  res
}
