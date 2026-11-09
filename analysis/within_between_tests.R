# =============================================================================
# within_between_tests.R — should GDP and conflict be split into within- and
# between-country parts?
#
#   Rscript analysis/within_between_tests.R [config.yml]
#   (after run_pipeline.R, migrant_backbone.R, refugee_role_changes.R)
#
# The model with one coefficient per driver assumes that a driver has the same
# association with a role within a country over time as between countries:
#     b_W = b_B    (Bafumi and Gelman 2006)
# For the main role and role-change models (drivers measured in the year of the
# change), this script fits both versions by maximum likelihood (Laplace, glmmTMB)
# and reports a Wald test of b_W = b_B for each driver and jointly, and AIC.
#
# Writes to <output_dir>/robustness/:
#   within_between_estimates.csv   odds ratios within and between, test for each driver
#   within_between_tests.csv       joint test and AIC for each model
# =============================================================================
# The Wald tests use model-based standard errors, which are too small for within-
# country terms (see part "bootstrap" of regression_robustness.R), so they are approximate.
# Reuse the model definitions: run regression_robustness.R up to its "baseline" section.
src <- readLines("analysis/regression_robustness.R")
src <- src[seq_len(grep("^# ---- baseline", src)[1] - 1)]
eval(parse(text = src))
built <- build_models(); models <- built$models
ests <- list(); tests <- list()
for (id in names(models)[!grepl("\\(change\\)", names(models))]) {
  m <- models[[id]]; d <- m$d; d$period_f <- droplevels(factor(d$period))
  cov_split <- usable(m$covs, d)
  vars <- sub("_within$", "", grep("_within$", cov_split, value = TRUE))
  # the same model with one coefficient per driver (the standardized driver itself)
  cov_one <- c(setdiff(cov_split, c(paste0(vars, "_within"), paste0(vars, "_between"))), vars)
  time <- time_terms(m$time)
  ctl <- glmmTMB::glmmTMBControl(parallel = 1)
  fit <- function(cv) glmmTMB::glmmTMB(as.formula(paste("y ~", paste(c(cv, time), collapse = " + "), "+ (1 | iso3)")),
                                       data = d, family = binomial, control = ctl)
  f_one <- fit(cov_one); f_split <- fit(cov_split)
  est <- glmmTMB::fixef(f_split)$cond; V <- as.matrix(vcov(f_split)$cond); se <- sqrt(diag(V))
  for (v in vars) {
    w <- paste0(v, "_within"); bt <- paste0(v, "_between")
    dif <- est[[bt]] - est[[w]]; se_dif <- sqrt(V[bt, bt] + V[w, w] - 2 * V[bt, w])
    ests[[paste(id, v)]] <- tibble(model = id, driver = v,
      or_within = exp(est[[w]]), within_low = exp(est[[w]] - 1.96 * se[[w]]), within_high = exp(est[[w]] + 1.96 * se[[w]]),
      or_between = exp(est[[bt]]), between_low = exp(est[[bt]] - 1.96 * se[[bt]]), between_high = exp(est[[bt]] + 1.96 * se[[bt]]),
      p_equal = 2 * pnorm(-abs(dif / se_dif)))
  }
  # joint Wald test of b_B = b_W for all split drivers: one row of contrasts per driver,
  # chi2 = (L b)' (L V L')^-1 (L b) with one degree of freedom per driver
  L <- matrix(0, length(vars), length(est), dimnames = list(NULL, names(est)))
  for (i in seq_along(vars)) { L[i, paste0(vars[i], "_between")] <- 1; L[i, paste0(vars[i], "_within")] <- -1 }
  chi <- as.numeric(t(L %*% est) %*% solve(L %*% V %*% t(L)) %*% (L %*% est))
  tests[[id]] <- tibble(model = id, n = nrow(d), events = sum(d$y),
    converged_one = f_one$fit$convergence == 0 && isTRUE(f_one$sdr$pdHess),
    converged_split = f_split$fit$convergence == 0 && isTRUE(f_split$sdr$pdHess),
    aic_one = AIC(f_one), aic_split = AIC(f_split), aic_difference = AIC(f_split) - AIC(f_one),
    joint_chi2 = chi, joint_df = length(vars), joint_p = pchisq(chi, length(vars), lower.tail = FALSE))
}
write_csv(bind_rows(ests), file.path(outdir, "within_between_estimates.csv"))
write_csv(bind_rows(tests), file.path(outdir, "within_between_tests.csv"))
print(as.data.frame(bind_rows(tests) |> mutate(across(where(is.double), \(x) round(x, 3)))))
msg("Done: ", outdir)
