# =============================================================================
# regression_evidence.R — one evidence table for each association, no thresholds
#
#   Rscript analysis/regression_evidence.R [config.yml]
#   (after regression_robustness.R, regression_robustness_summary.R, regression_bayes.R,
#    and after running the pipeline with config_blended.yml, which writes output_blended/)
#
# Brings together, as estimates with uncertainty:
#   - the main model of the paper (mixed-effects logistic regression with
#     within-country and between-country terms; Bafumi and Gelman 2006)
#   - the earlier model with one coefficient per driver (config_blended.yml)
#   - the Bayesian multilevel model (weakly informative priors)
#   - partial pooling across coefficients of the same family of drivers
#     (Gelman, Hill and Yajima 2012), in place of a multiple-testing correction
#   - the range over all alternative analyses (multiverse; Steegen et al. 2016)
#   - sign and magnitude errors (Gelman and Carlin 2014)
#   - the country cluster bootstrap interval and the placebo p-value of the main model
#
# Partial pooling. For the K coefficients of one family (conflict, income, income
# terciles, SPEI drought, disasters), with estimate y_k and standard error s_k:
#     y_k ~ N(theta_k, s_k^2),   theta_k ~ N(0, tau^2),   tau ~ half-Cauchy(0, 1)
# The centre is 0 because the role models are contrasts whose signs differ. tau is
# the typical size of a real effect in the family, learned from the family's K
# estimates; a small tau pulls noisy estimates strongly towards zero. tau is
# integrated over a grid (Gelman et al., Bayesian Data Analysis, ch. 5). Estimates
# from different models share data and are treated as independent here, so the
# tau intervals are somewhat too narrow.
#
# Sign and magnitude errors. For a true effect of size D and standard error s:
#   power     probability that the estimate is more than 1.96 s from zero
#   type S    probability that such an estimate has the wrong sign
#   type M    average exaggeration: E(|estimate| given significance) / D
# D is set to the family's tau, the typical real effect in that family.
#
# Writes to <output_dir>/robustness/: partial_pooling.csv, family_tau.csv,
#   evidence_table.csv, Fig_evidence.png
# =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(purrr); library(readr); library(ggplot2); library(stringr) })
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
cfg_file <- commandArgs(trailingOnly = TRUE)
cfg <- read_config(if (length(cfg_file)) cfg_file[1] else "config.yml")
outdir <- dirname(out_path(cfg, "robustness", "x"))
rd <- function(f) read_csv(file.path(outdir, f), show_col_types = FALSE)
set.seed(cfg$seed)

base <- rd("baseline.csv") |> filter(!grepl("^wald_equal|^joint_income", term))
bayes <- rd("bayes/coefficients.csv")

# ---- partial pooling ---------------------------------------------------------
family_of <- function(t) case_when(
  grepl("drought12|spei48", t) ~ "SPEI drought",
  grepl("conflict.*between", t) ~ "conflict, between countries",
  grepl("conflict", t) ~ "conflict, within country",
  grepl("^income_g", t) ~ "income tercile",
  grepl("gdp_pc.*between", t) ~ "GDP per capita, between countries",
  grepl("gdp_pc", t) ~ "GDP per capita, within country",
  TRUE ~ "disasters (EM-DAT)")
base <- base |> mutate(family = family_of(term))
tau_grid <- seq(0.001, 2, length.out = 2000)
pool_family <- function(d) {
  # log posterior of tau: half-Cauchy prior times the marginal likelihood of the estimates
  lp <- vapply(tau_grid, function(tau) dcauchy(tau, 0, 1, log = TRUE) + sum(dnorm(d$estimate, 0, sqrt(d$se^2 + tau^2), log = TRUE)), 0)
  w <- exp(lp - max(lp)); w <- w / sum(w)
  tau_draw <- sample(tau_grid, 4000, TRUE, w)
  theta <- vapply(seq_len(nrow(d)), function(k) {
    v <- 1 / (1 / d$se[k]^2 + 1 / tau_draw^2)
    rnorm(length(tau_draw), v * d$estimate[k] / d$se[k]^2, sqrt(v))
  }, numeric(length(tau_draw)))
  list(tab = d |> mutate(pool_median = apply(theta, 2, median), pool_low = apply(theta, 2, quantile, .025), pool_high = apply(theta, 2, quantile, .975),
                         pool_pr_pos = colMeans(theta > 0)),
       tau = tibble(family = d$family[1], k = nrow(d), tau_median = median(tau_draw), tau_low = quantile(tau_draw, .025), tau_high = quantile(tau_draw, .975)))
}
pooled <- map(split(base, base$family), pool_family)
pp <- bind_rows(map(pooled, "tab")); taus <- bind_rows(map(pooled, "tau"))
write_csv(pp, file.path(outdir, "partial_pooling.csv")); write_csv(taus, file.path(outdir, "family_tau.csv"))

# ---- sign and magnitude errors -----------------------------------------------
# D = assumed true effect, s = standard error; n = simulated estimates for type M
retro <- function(D, s, n = 1e5) {
  z <- qnorm(.975); l <- D / s
  hi <- 1 - pnorm(z - l); lo <- pnorm(-z - l); pow <- hi + lo
  est <- D + s * rnorm(n); sig <- abs(est) > z * s
  c(power = pow, type_s = lo / pow, type_m = mean(abs(est[sig])) / D)
}
rt <- pp |> left_join(taus |> select(family, tau_median), by = "family") |> rowwise() |>
  mutate(r = list(retro(tau_median, se))) |> ungroup() |>
  mutate(power = map_dbl(r, 1), type_s = map_dbl(r, 2), type_m = map_dbl(r, 3)) |> select(-r)

# ---- range over alternative analyses -----------------------------------------
# Analyses on the log-odds scale that enter the range (the same outcome and SD units).
# The "wb" specification exists only for configs with one coefficient per driver.
logit_specs <- c("adaptive quadrature (10 points)", "decade effects", "linear time trend", "drivers lagged one period",
                 "no disaster deaths or affected", "no conflict variables", "no SPEI drought variables",
                 "no current GDP (terciles only)", "one drought slope, no income terciles", "Firth logit (no random intercept)",
                 "conditional logit (within country)")
paths <- bind_rows(
  rd("spec_refits.csv") |> filter(spec %in% logit_specs) |> transmute(model, term, analysis = spec, b = estimate),
  rd("data_refits.csv") |> transmute(model, term, analysis = version, b = estimate),
  bayes |> filter(spec != "wb") |> transmute(model, term, analysis = paste("Bayes", spec), b = median))
sum_paths <- paths |> inner_join(base |> select(model, term, b0 = estimate), by = c("model", "term")) |>
  group_by(model, term) |> summarise(n_analyses = n(), or_low_path = exp(min(b)), or_high_path = exp(max(b)),
                                     share_same_sign = mean(sign(b) == sign(b0)), .groups = "drop")

# ---- evidence table ----------------------------------------------------------
# The associations stated in the Results, flagged as `claim` in the table (Table S9, Fig. S12).
claims <- bind_rows(
  tidyr::expand_grid(model = c("R1", "R2", "R3", "M1", "M2"),
                     term = c("conflict_events_within", "conflict_events_between", "conflict_deaths_within", "conflict_deaths_between",
                              "gdp_pc_within", "gdp_pc_between")),
  tribble(~model, ~term,
    "start sending (level)", "conflict_deaths_within", "start sending (level)", "conflict_deaths_between", "start sending (level)", "drought12_ext_x_g1",
    "start sending (change)", "d_conflict_deaths", "start sending (change)", "d_drought12_ext",
    "stop sending (level)", "conflict_deaths_between",
    "start hosting (level)", "conflict_deaths_within", "start hosting (level)", "conflict_deaths_between", "start hosting (level)", "income_g3"))
b1 <- bayes |> filter(spec == "main") |> transmute(model, term, bayes_m = median, bayes_l = q2.5, bayes_h = q97.5, pr_pos)
# the earlier model with one coefficient per driver (run config_blended.yml first)
bl_dir <- file.path(cfg$root, "output_blended")
blended <- bind_rows(
  read_csv(file.path(bl_dir, "unhcr", "binary_role_models.csv"), show_col_types = FALSE) |> transmute(model = paste0("R", role), term, blended = estimate),
  read_csv(file.path(bl_dir, "undesa_backbone", "binary_role_models.csv"), show_col_types = FALSE) |> transmute(model = paste0("M", role), term, blended = estimate),
  read_csv(file.path(bl_dir, "unhcr", "role_changes", "transition_models.csv"), show_col_types = FALSE) |>
    transmute(model = paste0(change, ifelse(spec == "level in year t", " (level)", " (change)")), term, blended = estimate))
ev <- rt |> transmute(model, term, family, b0 = estimate, se0 = se, p0 = p, pool_median, pool_low, pool_high, pool_pr_pos, tau_median, power, type_s, type_m) |>
  left_join(b1, by = c("model", "term")) |> mutate(base_term = sub("_(within|between)$", "", term)) |>
  left_join(blended, by = c("model", "base_term" = "term")) |>
  left_join(sum_paths, by = c("model", "term")) |>
  # uncertainty of the main model that does not rely on its standard errors
  left_join(rd("bootstrap.csv") |> select(model, term, n_boot, boot_low, boot_high), by = c("model", "term")) |>
  left_join(rd("placebo_summary.csv") |> filter(test == "country reassignment") |> transmute(model, term, p_placebo = p_perm, placebo_reject = rate_reject_05), by = c("model", "term")) |>
  mutate(claim = paste(model, term) %in% paste(claims$model, claims$term))
write_csv(ev, file.path(outdir, "evidence_table.csv"))

# ---- figure ------------------------------------------------------------------
model_lab <- c(R1 = "R1 sender only", R2 = "R2 two-way", R3 = "R3 receiver only", M1 = "M1 source", M2 = "M2 intermediate",
               `start sending (level)` = "Start sending (level)", `start sending (change)` = "Start sending (change)", `stop sending (level)` = "Stop sending (level)", `start hosting (level)` = "Start hosting (level)")
term_lab <- c(gdp_pc_within = "GDP per capita, within", gdp_pc_between = "GDP per capita, between", income_g3 = "Income tercile: high",
              conflict_events_within = "Conflict events, within", conflict_events_between = "Conflict events, between",
              conflict_deaths_within = "Conflict deaths, within", conflict_deaths_between = "Conflict deaths, between",
              drought12_ext_x_g1 = "Drought extent, low income", d_conflict_deaths = "Rise in conflict deaths", d_drought12_ext = "Rise in drought extent")
long <- ev |> filter(claim) |> transmute(model, term,
    `Main model_m` = b0, `Main model_l` = b0 - 1.96 * se0, `Main model_h` = b0 + 1.96 * se0,
    `Country bootstrap_m` = b0, `Country bootstrap_l` = boot_low, `Country bootstrap_h` = boot_high,
    `Bayesian multilevel_m` = bayes_m, `Bayesian multilevel_l` = bayes_l, `Bayesian multilevel_h` = bayes_h,
    `Partially pooled_m` = pool_median, `Partially pooled_l` = pool_low, `Partially pooled_h` = pool_high) |>
  pivot_longer(-c(model, term), names_to = c("est", ".value"), names_sep = "_") |> filter(!is.na(m)) |>
  mutate(est = factor(est, rev(c("Main model", "Country bootstrap", "Bayesian multilevel", "Partially pooled"))),
         panel = factor(model_lab[model], unname(model_lab)), label = factor(term_lab[term], rev(unname(term_lab))))
cols <- c("Main model" = "#3b3b3b", "Country bootstrap" = "#c4602d", "Bayesian multilevel" = "#2f6fb5", "Partially pooled" = "#1aa08a")
p <- ggplot(long, aes(y = label, x = exp(m), xmin = exp(l), xmax = exp(h), colour = est, group = est)) +
  geom_vline(xintercept = 1, colour = "grey60") +
  geom_pointrange(position = position_dodge(width = 0.8), size = 0.18, linewidth = 0.45) +
  scale_x_log10(breaks = c(1/32, 1/8, 1/2, 2, 8, 32), labels = c("1/32", "1/8", "1/2", "2", "8", "32")) +
  coord_cartesian(xlim = c(0.025, 40)) +
  scale_colour_manual(values = cols, breaks = names(cols)) +
  facet_wrap(~panel, scales = "free_y", ncol = 3) +
  labs(x = "Odds ratio per SD (log scale) with 95% interval (intervals beyond the axis are cut at its edge)", y = NULL, colour = NULL) +
  theme_minimal(9) + theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave(file.path(outdir, "Fig_evidence.png"), p, width = 9.5, height = 9, dpi = 300, bg = "white")
print(taus |> mutate(across(where(is.double), \(x) round(x, 3))) |> as.data.frame())
msg("Done: ", outdir)
