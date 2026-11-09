# =============================================================================
# regression_bayes.R — Bayesian multilevel versions of the role regressions
#
#   Rscript analysis/regression_bayes.R [config.yml]
#   (after run_pipeline.R, migrant_backbone.R, refugee_role_changes.R and
#    regression_robustness.R; about 20 minutes on 2 cores)
#   Environment: BAYES_MODELS (ids, comma-separated; default all), BAYES_ITER (default 1000)
#
# The analysis follows the approach of Gelman and coauthors:
#   - every model has a country intercept (priors make all of them estimable,
#     including those fitted as plain logistic regression in the paper);
#   - role models: a linear time trend plus varying intercepts for periods,
#     instead of a separate unpooled effect for each period;
#   - GDP and conflict enter as a within-country and a between-country
#     component (Bafumi and Gelman 2006), as in the main model:
#         x_ct = xbar_c + (x_ct - xbar_c);
#     the coefficient on (x_ct - xbar_c) answers "when a country has more of x
#     than usual", the coefficient on xbar_c "when a country has more of x than
#     other countries";
#   - weakly informative priors on predictors in SD units: normal(0, 1) on the
#     log-odds scale, so a 95% prior interval for an odds ratio per SD of
#     0.14 to 7; normal(0, 2.5) is the prior sensitivity check;
#   - convergence (R-hat, effective sample size, divergences) and posterior
#     predictive checks are recorded for every model.
#
# Specifications (spec):
#   main      same predictors as the main model in the paper
#   prior2.5  main specification with normal(0, 2.5) priors
#   wb        only when the config has one coefficient per driver (config_blended.yml):
#             the same model with the drivers split into within and between terms
#
# These models are a check on the main model. They share its likelihood, and
# therefore also treat a country's successive years as independent given its intercept.
#
# Writes to <output_dir>/robustness/bayes/:
#   coefficients.csv   posterior median, 50% and 95% intervals, Pr(b > 0) for each coefficient
#   diagnostics.csv    max R-hat, min effective sample size, divergences
#   ppc.csv            posterior predictive checks
#   loo.csv            approximate leave-one-out comparison of main and wb
#                      (empty for the main config, which has no wb specification)
# =============================================================================
# Reuse the model definitions: run regression_robustness.R up to its "baseline" section.
source_head <- readLines("analysis/regression_robustness.R")
source_head <- source_head[seq_len(grep("^# ---- baseline", source_head)[1] - 1)]
eval(parse(text = source_head))
suppressPackageStartupMessages(library(rstanarm))
built <- build_models(); models <- built$models
base <- read_csv(file.path(outdir, "baseline.csv"), show_col_types = FALSE)
bdir <- file.path(outdir, "bayes"); dir.create(bdir, showWarnings = FALSE)
ids <- strsplit(Sys.getenv("BAYES_MODELS", paste(names(models), collapse = ",")), ",")[[1]]
iter <- as.integer(Sys.getenv("BAYES_ITER", "1000"))
split_vars <- c("gdp_pc", "conflict_events", "conflict_deaths")

# formula: predictors + time + country intercept (+ period intercept for role models)
fml <- function(m, covs) {
  # a linear trend in every model: GDP is in current US$ and rises over time, so
  # period intercepts shrunk towards a constant could not absorb that trend alone
  rhs <- c(covs, "period_c")
  res <- c("(1 | iso3)", if (m$time == "factor") "(1 | period_f)")
  as.formula(paste("y ~", paste(c(rhs, res), collapse = " + ")))
}
# Fit one model in rstanarm. scale = SD of the normal prior on the predictors.
# adapt_delta is raised to 0.95 to limit divergent transitions.
fit_b <- function(f, d, scale, seed) {
  suppressWarnings(stan_glmer(f, data = d, family = binomial(), prior = normal(0, scale), prior_intercept = normal(0, 2.5),
    prior_aux = NULL, chains = 4, cores = cores, iter = iter, seed = seed, adapt_delta = 0.95, refresh = 0, QR = FALSE))
}
# Posterior summaries of the predictors (no country or period intercepts, no trend)
tidy_b <- function(fit, m, spec) {
  a <- as.matrix(fit); keep <- !grepl("^b\\[|^Sigma\\[|Intercept|^period_c$", colnames(a)); a <- a[, keep, drop = FALSE]
  tibble(model = m$id, spec = spec, term = colnames(a), median = apply(a, 2, median),
         q2.5 = apply(a, 2, quantile, .025), q25 = apply(a, 2, quantile, .25), q75 = apply(a, 2, quantile, .75), q97.5 = apply(a, 2, quantile, .975),
         pr_pos = colMeans(a > 0))
}
# Convergence diagnostics for one fit
diag_b <- function(fit, m, spec, secs) {
  s <- summary(fit); np <- bayesplot::nuts_params(fit)
  tibble(model = m$id, spec = spec, max_rhat = max(s[, "Rhat"], na.rm = TRUE), min_neff = min(s[, "n_eff"], na.rm = TRUE),
         divergences = sum(np$Value[np$Parameter == "divergent__"]), minutes = secs / 60)
}
# Posterior predictive checks: observed statistic vs. simulated replications
ppc_b <- function(fit, m, spec) {
  yrep <- posterior_predict(fit, draws = 400); d <- m$d
  stat <- list(
    `share of cases` = function(y) mean(y),
    `share of countries never in the outcome` = function(y) mean(tapply(y, d$iso3, sum) == 0),
    `share of countries always in the outcome` = function(y) mean(tapply(y, d$iso3, mean) == 1),
    `SD across countries of the outcome rate` = function(y) sd(tapply(y, d$iso3, mean)),
    `SD across periods of the outcome rate` = function(y) sd(tapply(y, d$period, mean)))
  imap_dfr(stat, function(f, nm) { r <- apply(yrep, 1, f); o <- f(d$y)
    tibble(model = m$id, spec = spec, stat = nm, observed = o, rep_median = median(r), rep_low = quantile(r, .025), rep_high = quantile(r, .975), ppp = mean(r >= o)) })
}
# Within / between predictors, for configs with one coefficient per driver. With the
# main config, m$covs already holds the _within and _between terms and this is not used.
wb_data <- function(m) {
  v <- intersect(split_vars, m$covs)
  d <- m$d |> group_by(iso3) |> mutate(across(all_of(v), list(between = mean, within = \(x) x - mean(x)), .names = "{.col}__{.fn}")) |> ungroup()
  list(d = d, covs = c(setdiff(m$covs, v), paste0(v, "__within"), paste0(v, "__between")))
}
out <- list(); dg <- list(); pp <- list(); loo_r <- list()   # results, diagnostics, predictive checks, LOO
for (id in ids) {
  m <- models[[id]]; m$d$period_f <- droplevels(factor(m$d$period)); covs <- usable(m$covs, m$d)
  is_level <- !grepl("\\(change\\)", id)
  t0 <- Sys.time(); f1 <- fit_b(fml(m, covs), m$d, 1, cfg$seed); sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  out[[paste(id, "main")]] <- tidy_b(f1, m, "main"); dg[[paste(id, "main")]] <- diag_b(f1, m, "main", sec); pp[[paste(id, "main")]] <- ppc_b(f1, m, "main")
  msg("  ", id, " main: ", round(sec / 60, 1), " min; max R-hat ", round(dg[[paste(id, "main")]]$max_rhat, 3), "; divergences ", dg[[paste(id, "main")]]$divergences)
  if (length(intersect(split_vars, m$covs))) {
    w <- wb_data(m); covw <- usable(w$covs, w$d)
    m2 <- m; m2$d <- w$d
    t0 <- Sys.time(); f2 <- fit_b(fml(m2, covw), w$d, 1, cfg$seed); sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    out[[paste(id, "wb")]] <- tidy_b(f2, m2, "wb"); dg[[paste(id, "wb")]] <- diag_b(f2, m2, "wb", sec); pp[[paste(id, "wb")]] <- ppc_b(f2, m2, "wb")
    msg("  ", id, " wb: ", round(sec / 60, 1), " min; max R-hat ", round(dg[[paste(id, "wb")]]$max_rhat, 3))
    l <- tryCatch(loo::loo_compare(rstanarm::loo(f1), rstanarm::loo(f2)), error = function(e) NULL)
    if (!is.null(l)) loo_r[[id]] <- tibble(model = id, spec = rownames(l), elpd_diff = l[, "elpd_diff"], se_diff = l[, "se_diff"]) |>
        mutate(spec = recode(spec, model1 = "main", model2 = "wb", f1 = "main", f2 = "wb"))
  }
  t0 <- Sys.time(); f3 <- fit_b(fml(m, covs), m$d, 2.5, cfg$seed + 1); sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  out[[paste(id, "p25")]] <- tidy_b(f3, m, "prior2.5"); dg[[paste(id, "p25")]] <- diag_b(f3, m, "prior2.5", sec)
  write_csv(bind_rows(out), file.path(bdir, "coefficients.csv")); write_csv(bind_rows(dg), file.path(bdir, "diagnostics.csv"))
  write_csv(bind_rows(pp), file.path(bdir, "ppc.csv")); write_csv(bind_rows(loo_r), file.path(bdir, "loo.csv"))
}
msg("Done: ", bdir)
