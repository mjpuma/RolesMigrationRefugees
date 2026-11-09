# =============================================================================
# regression_robustness.R — are the role and role-change regressions robust?
#
#   Rscript analysis/regression_robustness.R [config.yml]
#   (after run_pipeline.R, analysis/migrant_backbone.R, analysis/refugee_role_changes.R;
#    part "data" also reads the outputs of the alternative config files, so run those
#    first; about 100 minutes on 2 cores for all parts)
#
#   Environment variables:
#     ROBUST_PARTS   any of: sample spec data inference placebo bootstrap   (default: all)
#     ROBUST_CORES   number of cores (default 2)
#     ROBUST_NPERM   number of placebo runs (default 500)
#     ROBUST_NBOOT   number of bootstrap samples (default 500)
#
#   Packages beyond the main pipeline: glmmTMB, sandwich, survival, foreign, brglm2 (Firth logit only).
#
# There was no written protocol for this paper, so the checks follow the
# "Summary of robustness checks" table in the supplement of Esmaili, Puma et al.,
# "El Nino amplified food insecurity in early modern Europe":
#
#   El Nino paper                         here
#   ------------------------------------  -----------------------------------------
#   Sample composition                    part "sample"
#     leave out each location             leave out each country
#     leave out each famine episode       leave out each event (role changes)
#     drop 20-year blocks                 drop 5-year blocks / each 5-year period
#                                         also: each subregion, each continent,
#                                         countries under 1 million people, 2020
#   Model specification                   part "spec"
#     parametric alternatives             adaptive quadrature; linear probability
#                                         model; conditional logit; Firth logit;
#                                         multinomial logit
#     distributed lags                    drivers lagged one period
#     nonlinearity                        quadratic and cubic terms
#     controls                            drop each group of drivers; time effects
#   Alternative data                      part "data"
#     other price / harvest datasets      other drought measures (alternative
#                                         configs); other role definitions
#   Inference                             part "inference"
#     Conley, spatial HAC, two-way        cluster by country, by subregion, two-way
#                                         (country x period), Conley spatial-serial;
#                                         part "bootstrap": country cluster bootstrap
#                                         of the mixed model itself
#   Identification                        part "placebo"
#     shuffle ENSO, N = 500               reassign each country's whole driver
#                                         history to another country, N = 500;
#                                         shuffle driver years within country
#
# The models examined (same data and fitting rules as the paper):
#   R1 R2 R3   refugee roles             (output/unhcr/model_data.csv)
#   M1 M2      migrant backbone roles    (M3 has 12 cases and is not interpreted)
#   start sending / stop sending / start hosting, each with drivers measured
#              in the year of the change ("level") and as the change from the
#              previous year ("change")
#
# Conventions used throughout:
#   - The main model is a mixed-effects logistic regression with within-country
#     and between-country terms for GDP per capita and conflict (R/models.R).
#   - Refits on a subsample keep the standardization, the country means behind the
#     within/between terms and the income terciles of the full sample, so that all
#     odds ratios stay "per the same SD" and each refit estimates the same quantity.
#   - A main model that fell back to ordinary logistic regression (no acceptable
#     mixed fit) would be refitted the same way here (m$engine_glm). None does at
#     present; M3, with 12 cases, is not examined.
#   - Part "bootstrap" is not only a check: its percentile intervals are the 95%
#     intervals shown in Figs. 3 and 4 and quoted in the Results.
#
# Writes to <output_dir>/robustness/:
#   baseline.csv               the estimates all checks are compared with
#   sample_refits.csv          every refit of part "sample"
#   spec_refits.csv            every refit of part "spec"
#   nonlinearity.csv           joint tests of quadratic / cubic terms
#   vif.csv                    variance inflation factors
#   multinomial.csv            multinomial logit (one row per equation and term)
#   data_refits.csv            alternative drought measures and role definitions
#   inference.csv              standard errors under each assumption
#   placebo_draws.csv          t statistics of every placebo run
#   placebo_summary.csv        randomization p-values
#   bootstrap.csv              country cluster bootstrap of the mixed models: SE, 95%
#                              percentile interval (log odds), refits used (n_boot)
#   summary_table.csv          one row per model and driver, all checks combined
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(readr); library(ggplot2)
  library(survival)      # conditional logit; strata() in the formula
})
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
cfg_file <- commandArgs(trailingOnly = TRUE)
cfg <- read_config(if (length(cfg_file)) cfg_file[1] else "config.yml")
set.seed(cfg$seed)

parts <- strsplit(Sys.getenv("ROBUST_PARTS", "sample spec data inference placebo bootstrap"), " +")[[1]]
cores <- as.integer(Sys.getenv("ROBUST_CORES", "2"))
nperm <- as.integer(Sys.getenv("ROBUST_NPERM", "500"))
outdir <- dirname(out_path(cfg, "robustness", "x"))

# Parallel map over a list. An element that throws an error (for example a refit
# that does not converge) returns NULL and is dropped, so one failure does not
# stop the others; callers report how many elements were used.
pmap_safe <- function(x, f) {
  res <- parallel::mclapply(x, function(el) tryCatch(f(el), error = function(e) NULL),
                            mc.cores = cores, mc.preschedule = TRUE)
  res[!vapply(res, function(r) is.null(r) || inherits(r, "try-error"), logical(1))]
}

# ---- models ------------------------------------------------------------------

# One model to examine:
#   id, group    short name, and "refugee role" / "migrant role" / "role change"
#   d            data with outcome y (one row per country-period at risk)
#   pool         all rows with driver values (placebo donors, lags)
#   covs, slopes covariates; tercile slope columns for each drought variable
#   time         "factor" (period effects) or "linear" (trend), as in the paper
make_model <- function(id, group, d, pool, covs, slopes, time) {
  d$period_f <- factor(d$period); d$iso3 <- as.character(d$iso3)
  pool$iso3 <- as.character(pool$iso3)
  list(id = id, group = group, d = d, pool = pool, covs = covs, slopes = slopes, time = time,
       step = min(diff(sort(unique(pool$period)))))
}

build_models <- function() {
  areas <- read_csv(in_root(cfg, cfg$country_list), show_col_types = FALSE) |>
    transmute(iso3 = iso3c, area = un_area, subregion = un_region)
  add_area <- function(d) d |> select(-any_of(c("area", "subregion"))) |> left_join(areas, by = "iso3") |>
    mutate(area = coalesce(area, "Other"), subregion = coalesce(subregion, "Other"))
  out <- list(); mds <- list()

  # refugee roles: the model data written by run_pipeline.R
  md <- read_csv(out_path(cfg, "unhcr", "model_data.csv"), show_col_types = FALSE) |> add_area() |>
    mutate(period_f = factor(period))        # the CSV stores it as a number
  covs_r <- readLines(out_path(cfg, "unhcr", "model_covariates.txt"))
  sg <- read_csv(out_path(cfg, "unhcr", "model_slope_groups.csv"), show_col_types = FALSE)
  slopes_r <- split(sg$column, sg$variable)
  mds$refugee <- list(data = md, covariates = covs_r, slope_cols = slopes_r)
  for (r in sort(unique(md$role)))
    out[[paste0("R", r)]] <- make_model(paste0("R", r), "refugee role",
                                        md |> mutate(y = as.integer(role == r)), md, covs_r, slopes_r, "factor")

  # migrant backbone roles: same call as analysis/migrant_backbone.R
  roles_bb <- read_csv(out_path(cfg, "undesa_backbone", "backbone_roles.csv"), show_col_types = FALSE)
  drivers_m <- read_csv(out_path(cfg, "undesa", "drivers_by_country_period.csv"), show_col_types = FALSE)
  mdb <- prepare_model_data(cfg, roles_bb, drivers_m)
  mdb$data <- add_area(mdb$data)
  mds$migrant <- mdb
  for (r in sort(unique(mdb$data$role))) {
    y <- as.integer(mdb$data$role == r)
    # M3 (terminal destinations, 12 cases) is not interpreted in the paper: too few to test
    if (sum(y) >= 30 && sum(1 - y) >= 30)
      out[[paste0("M", r)]] <- make_model(paste0("M", r), "migrant role",
                                          mdb$data |> mutate(y = y), mdb$data, mdb$covariates, mdb$slope_cols, "factor")
  }

  # refugee role changes: same construction as analysis/refugee_role_changes.R
  covs_level <- covs_r[covs_r %in% names(md)]
  covs_base <- covs_level[!grepl("_lag1|_x_|^income_g|_base$", covs_level)]
  covs_base <- unique(sub("_(within|between)$", "", covs_base))
  covs_base <- unique(c(covs_base, setdiff(names(slopes_r), grep("_lag1$", names(slopes_r), value = TRUE))))
  covs_base <- covs_base[covs_base %in% names(md)]
  covs_base <- covs_base[vapply(covs_base, \(v) sd(md[[v]], na.rm = TRUE) > 0, logical(1))]
  d <- md |> arrange(iso3, period) |> group_by(iso3) |>
    mutate(prev_role = lag(role), prev_period = lag(period),
           across(all_of(covs_base), \(x) x - lag(x), .names = "d_{.col}")) |> ungroup() |>
    filter(!is.na(prev_role), period == prev_period + 1)
  changes <- list(`start sending` = list(a = 3, e = c(1, 2)), `stop sending` = list(a = c(1, 2), e = 3),
                  `start hosting` = list(a = 1, e = c(2, 3)))
  dvars <- paste0("d_", covs_base)
  for (nm in names(changes)) {
    ch <- changes[[nm]]
    dd <- d |> filter(prev_role %in% ch$a) |> mutate(y = as.integer(role %in% ch$e))
    out[[paste0(nm, " (level)")]] <- make_model(paste0(nm, " (level)"), "role change", dd, d, covs_level, slopes_r, "linear")
    dc <- dd |> filter(if_all(all_of(dvars), \(x) !is.na(x)))
    out[[paste0(nm, " (change)")]] <- make_model(paste0(nm, " (change)"), "role change", dc,
                                                 d |> filter(if_all(all_of(dvars), \(x) !is.na(x))), dvars, list(), "linear")
  }
  list(models = out, mds = mds)
}

# ---- fitting -----------------------------------------------------------------

# covariates that can be estimated in this data: numeric with variation (or a factor with > 1 level)
usable <- function(covs, d) covs[vapply(covs, function(v) {
  if (!v %in% names(d)) return(FALSE)
  x <- d[[v]]
  if (is.factor(x)) return(nlevels(droplevels(x)) > 1)
  s <- sd(x, na.rm = TRUE); !is.na(s) && s > 0
}, logical(1))]

# formula term for the time effect: period fixed effects, a linear trend, or none
time_terms <- function(time) switch(time, factor = "period_f", linear = "period_c", decade = "decade_f", none = character(0))

# Fit model m (optionally on other data, covariates or time effects).
#   fast = TRUE   glmmTMB, the same engine as the paper's main fits (R/models.R,
#                 fit_mixed), with fit_glmer_robust() if it does not converge
#   fast = FALSE  lme4 through fit_glmer_robust() (analysis/glmmTMB_vs_lme4.R)
#   nAGQ > 1      lme4 with adaptive quadrature (a specification check)
# If the paper's main fit was ordinary logistic regression, m$engine_glm = TRUE
# and every refit of that model is ordinary logistic too.
fit_m <- function(m, d = m$d, covs = m$covs, time = m$time, nAGQ = 1, fast = TRUE) {
  d$period_f <- droplevels(factor(d$period))
  covs <- usable(covs, d)
  tt <- if (time == "decade") { "decade_f" } else time_terms(time)
  rhs <- c(covs, tt); if ("decade_f" %in% rhs) d$decade_f <- droplevels(factor(floor((pmax(d$period, 1991) - 1) / 10) * 10))  # 1990s (incl. 1990), 2000s, 2010s (incl. 2020)
  f_g <- as.formula(paste("y ~", paste(rhs, collapse = " + ")))
  f_m <- as.formula(paste("y ~", paste(rhs, collapse = " + "), "+ (1 | iso3)"))
  if (isTRUE(m$engine_glm)) {
    fit <- glm(f_g, data = d, family = binomial); attr(fit, "engine") <- "logistic"; return(fit)
  }
  if (nAGQ > 1) {
    fit <- suppressWarnings(lme4::glmer(f_m, data = d, family = binomial, nAGQ = nAGQ,
                                        control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))))
    attr(fit, "engine") <- paste0("mixed (AGQ ", nAGQ, ")"); return(fit)
  }
  if (!fast) return(suppressMessages(fit_glmer_robust(f_m, f_g, d, label = m$id)))
  # refits: same model and Laplace likelihood as the main fit, in glmmTMB, which is
  # about 70 times faster than lme4 here (comparison: analysis/glmmTMB_vs_lme4.R);
  # a fit that does not converge goes through fit_glmer_robust()
  fit <- tryCatch(suppressWarnings(suppressMessages(
    glmmTMB::glmmTMB(f_m, data = d, family = binomial, control = glmmTMB::glmmTMBControl(parallel = 1)))),
    error = function(e) NULL)
  if (!is.null(fit) && fit$fit$convergence == 0 && isTRUE(fit$sdr$pdHess)) {
    attr(fit, "engine") <- "mixed (glmmTMB)"; return(fit)
  }
  suppressMessages(fit_glmer_robust(f_m, f_g, d, label = m$id))
}

# Coefficient table, estimates and covariance for any of the three fit classes
# (glmmTMB, lme4 merMod, glm). The same helpers are in R/models.R.
co_table <- function(fit) if (inherits(fit, "glmmTMB")) summary(fit)$coefficients$cond else summary(fit)$coefficients
b_vec <- function(fit) if (inherits(fit, "glmmTMB")) glmmTMB::fixef(fit)$cond else if (inherits(fit, "merMod")) lme4::fixef(fit) else coef(fit)
v_mat <- function(fit) if (inherits(fit, "glmmTMB")) as.matrix(vcov(fit)$cond) else as.matrix(vcov(fit))

# Joint Wald test that a set of coefficients are all zero.
wald_zero <- function(b, V, cols) {
  cols <- intersect(cols, names(b)); if (!length(cols)) return(c(chi2 = NA, df = NA, p = NA))
  chi2 <- tryCatch(as.numeric(t(b[cols]) %*% solve(V[cols, cols, drop = FALSE]) %*% b[cols]), error = function(e) NA_real_)
  c(chi2 = chi2, df = length(cols), p = pchisq(chi2, length(cols), lower.tail = FALSE))
}

# One row per covariate, plus Wald rows used in the paper:
#   wald_equal_slopes:<var>   equal drought slopes across income terciles
#   joint_income              GDP and the tercile terms jointly zero
tidy_m <- function(fit, m) {
  co <- co_table(fit)
  b <- b_vec(fit)
  b <- b[!is.na(b)]; V <- v_mat(fit)[names(b), names(b)]
  out <- tibble(term = rownames(co), estimate = co[, 1], se = co[, 2], p = co[, 4]) |>
    filter(!grepl("Intercept|^period_|^decade_", term))
  for (v in names(m$slopes)) {
    w <- wald_equal(b, V, m$slopes[[v]])
    out <- bind_rows(out, tibble(term = paste0("wald_equal_slopes:", v), estimate = w[["chi2"]], se = w[["df"]], p = w[["p"]]))
  }
  if ("gdp_pc" %in% names(b) && any(grepl("^income_g", names(b)))) {
    w <- wald_zero(b, V, c("gdp_pc", grep("^income_g", names(b), value = TRUE)))
    out <- bind_rows(out, tibble(term = "joint_income", estimate = w[["chi2"]], se = w[["df"]], p = w[["p"]]))
  }
  out |> mutate(engine = attr(fit, "engine"), singular = inherits(fit, "merMod") && lme4::isSingular(fit))
}

# fit_m() followed by tidy_m(), with the number of rows and cases actually used
fit_tidy <- function(m, ...) {
  d <- list(...)$d %||% m$d
  tidy_m(fit_m(m, ...), m) |> mutate(n = nrow(d), events = sum(d$y))
}

# ---- baseline ----------------------------------------------------------------
built <- build_models(); models <- built$models; mds <- built$mds
msg("Models: ", paste(names(models), collapse = ", "))
baseline <- map_dfr(names(models), function(id) {
  m <- models[[id]]
  f0 <- fit_m(m)                                    # same fitting rule as the paper (R/models.R, fit_mixed)
  models[[id]]$engine_glm <<- inherits(f0, "glm")   # TRUE only if the main fit fell back to ordinary logistic regression
  tidy_m(f0, m) |> mutate(model = id, group = m$group, n = nrow(m$d), events = sum(m$d$y))
})
write_csv(baseline, file.path(outdir, "baseline.csv"))
msg("Baseline fitted. Logistic (no random intercept): ",
    paste(names(models)[vapply(models, \(m) isTRUE(m$engine_glm), TRUE)], collapse = ", "))

# Baseline must reproduce the paper's results exactly (else the checks would be about another model).
check_repro <- function() {
  a <- bind_rows(
    read_csv(out_path(cfg, "unhcr", "binary_role_models.csv"), show_col_types = FALSE) |> mutate(model = paste0("R", role)),
    read_csv(out_path(cfg, "undesa_backbone", "binary_role_models.csv"), show_col_types = FALSE) |> mutate(model = paste0("M", role)),
    read_csv(out_path(cfg, "unhcr", "role_changes", "transition_models.csv"), show_col_types = FALSE) |>
      mutate(model = paste0(change, ifelse(spec == "level in year t", " (level)", " (change)")))) |>
    select(model, term, ref = estimate)
  x <- baseline |> inner_join(a, by = c("model", "term")) |> filter(!grepl("^wald", term))
  msg("Reproduction check: ", nrow(x), " coefficients; largest difference from the paper's output = ",
      signif(max(abs(x$estimate - x$ref)), 3))
  if (max(abs(x$estimate - x$ref)) > 1e-3) warning("Baseline does not reproduce the pipeline output")
}
check_repro()

# ---- A. sample composition ---------------------------------------------------
# Refit after removing each unit in turn. Each split = list(type, unit, keep).
splits_for <- function(m) {
  d <- m$d; sp <- list()
  add <- function(type, unit, keep) sp[[length(sp) + 1]] <<- list(type = type, unit = unit, keep = keep)
  for (u in unique(d$iso3)) add("country", u, d$iso3 != u)
  for (u in unique(d$subregion)) add("subregion", u, d$subregion != u)
  for (u in unique(d$area)) add("continent", u, d$area != u)
  if (m$step == 1) {                                    # annual refugee data: 5-year blocks
    blk <- 5 * ceiling(pmax(d$period, 1991) / 5)         # 1990 joins the first block
    for (u in sort(unique(blk))) add("time block", sprintf("%d-%d", if (u == 1995) 1990 else u - 4, u), blk != u)
    add("year 2020", "2020", d$period != 2020)
  } else for (u in sort(unique(d$period))) add("time block", sprintf("%d-%d", u - 4, u), d$period != u)
  # small states
  small <- pops$iso3[pops$pop < 1e6]
  add("small states", "under 1 million people", !(d$iso3 %in% small))
  # each event separately (the analogue of leaving out one famine)
  if (m$group == "role change") for (i in which(d$y == 1))
    add("event", paste(d$iso3[i], d$period[i]), seq_len(nrow(d)) != i)
  sp
}

if ("sample" %in% parts) {
  pops <- read_csv(out_path(cfg, "undesa", "drivers_by_country_period.csv"), show_col_types = FALSE) |>
    group_by(iso3) |> summarise(pop = median(pop, na.rm = TRUE))
  jobs <- flatten(imap(models, function(m, id) map(splits_for(m), \(s) c(s, id = id))))
  msg("Sample composition: ", length(jobs), " refits")
  t0 <- Sys.time()
  res <- pmap_safe(jobs, function(j) {
    m <- models[[j$id]]; d <- m$d[j$keep, ]
    if (sum(d$y) < 5 || sum(1 - d$y) < 5) return(NULL)
    fit_tidy(m, d = d) |> mutate(model = j$id, group = m$group, type = j$type, unit = j$unit)
  })
  write_csv(bind_rows(res), file.path(outdir, "sample_refits.csv"))
  msg("  done: ", length(res), " of ", length(jobs), " refits in ",
      round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")
}

# ---- B. model specification --------------------------------------------------
is_dis  <- function(v) grepl("(^|d_)(drought|flood|storm|landslide)_(affected|deaths)$", v)
is_conf <- function(v) grepl("^(d_)?conflict_", v)
is_spei <- function(v) grepl("(^|d_)(drought12|spei48)", v)

# Lag every time-varying driver by one observed period (year, or 5-year period).
lag_data <- function(m) {
  vars <- m$covs[!grepl("^income_g", m$covs)]
  lagged <- m$pool |> arrange(iso3, period) |> group_by(iso3) |>
    mutate(across(all_of(vars), \(x) ifelse(period - lag(period) == m$step, lag(x), NA))) |> ungroup() |>
    select(iso3, period, all_of(vars))
  m$d |> select(-all_of(vars)) |> left_join(lagged, by = c("iso3", "period")) |>
    filter(if_all(all_of(vars), \(x) !is.na(x)))
}

# common drought slope: the slope columns add up to the drought variable itself
common_slope <- function(m) c(setdiff(m$covs, c(unlist(m$slopes), grep("^income_g", m$covs, value = TRUE))), names(m$slopes))

# linear probability model with country fixed effects; SEs clustered by country
fit_lpm <- function(m, d = m$d, covs = m$covs, time = m$time, B = NULL) {
  d$period_f <- droplevels(factor(d$period))
  covs <- usable(covs, d)
  # terms that do not vary within country are absorbed by the fixed effects
  covs <- covs[vapply(covs, \(v) any(tapply(d[[v]], d$iso3, function(x) length(unique(x)) > 1)), TRUE)]
  f <- as.formula(paste("y ~", paste(c(covs, time_terms(time), "factor(iso3)"), collapse = " + ")))
  fit <- lm(f, data = d)
  b <- coef(fit); b <- b[!is.na(b)]
  V <- sandwich::vcovCL(fit, cluster = d$iso3, type = "HC1")
  se <- sqrt(diag(V))[names(b)]
  tibble(term = names(b), estimate = b, se = se, p = 2 * pnorm(-abs(b / se))) |>
    filter(!grepl("Intercept|^period_|factor\\(iso3\\)", term))
}

# conditional (country fixed-effects) logit: uses only within-country change
fit_clogit <- function(m, d = m$d, covs = m$covs, time = m$time) {
  d$period_f <- droplevels(factor(d$period))
  covs <- usable(covs, d)
  covs <- covs[vapply(covs, \(v) any(tapply(d[[v]], d$iso3, function(x) length(unique(x)) > 1)), TRUE)]
  f <- as.formula(paste("y ~", paste(c(covs, time_terms(time), "strata(iso3)"), collapse = " + ")))
  fit <- suppressWarnings(survival::clogit(f, data = d, method = "exact"))
  co <- summary(fit)$coefficients
  tibble(term = rownames(co), estimate = co[, "coef"], se = co[, "se(coef)"], p = co[, "Pr(>|z|)"]) |>
    filter(!grepl("^period_", term))
}

# Firth (bias-reduced) logit, no random intercept: for rare outcomes
fit_firth <- function(m) {
  d <- m$d; d$period_f <- droplevels(factor(d$period)); covs <- usable(m$covs, d)
  f <- as.formula(paste("y ~", paste(c(covs, time_terms(m$time)), collapse = " + ")))
  fit <- glm(f, data = d, family = binomial, method = brglm2::brglmFit, type = "AS_mean")
  co <- summary(fit)$coefficients
  tibble(term = rownames(co), estimate = co[, 1], se = co[, 2], p = co[, 4]) |> filter(!grepl("Intercept|^period_", term))
}

# label one table of estimates with its model, specification and sample size
add_info <- function(tb, m, spec, scale = "log-odds", d = m$d) tb |> mutate(model = m$id, group = m$group, spec = spec, scale = scale,
                                                                             n = nrow(d), events = sum(d$y))

if ("spec" %in% parts) {
  jobs <- flatten(imap(models, function(m, id) {
    is_level <- !grepl("\\(change\\)", id)
    sp <- list(
      `adaptive quadrature (10 points)` = \(m) add_info(fit_tidy(m, nAGQ = 10), m, "adaptive quadrature (10 points)"),
      `decade effects` = \(m) add_info(fit_tidy(m, time = "decade"), m, "decade effects"),
      `no disaster deaths or affected` = \(m) add_info(fit_tidy(m, covs = m$covs[!is_dis(m$covs)]), m, "no disaster deaths or affected"),
      `no conflict variables` = \(m) add_info(fit_tidy(m, covs = m$covs[!is_conf(m$covs)]), m, "no conflict variables"),
      `no SPEI drought variables` = \(m) add_info(fit_tidy(m, covs = m$covs[!is_spei(m$covs)]), m, "no SPEI drought variables"),
      `linear probability, country fixed effects` = \(m) add_info(fit_lpm(m), m, "linear probability, country fixed effects", "probability"),
      `conditional logit (within country)` = \(m) add_info(fit_clogit(m), m, "conditional logit (within country)"))
    if (m$time == "factor") sp$`linear time trend` <- \(m) add_info(fit_tidy(m, time = "linear"), m, "linear time trend")
    if (m$time == "linear" ) sp$`year effects` <- NULL
    if (is_level) {
      sp$`drivers lagged one period` <- \(m) { d <- lag_data(m); add_info(fit_tidy(m, d = d), m, "drivers lagged one period", d = d) }
      if (length(m$slopes)) sp$`one drought slope, no income terciles` <- \(m) add_info(fit_tidy(m, covs = common_slope(m)), m, "one drought slope, no income terciles")
      if (any(grepl("^gdp_pc", m$covs))) sp$`no current GDP (terciles only)` <- \(m) add_info(fit_tidy(m, covs = m$covs[!grepl("^gdp_pc", m$covs)]), m, "no current GDP (terciles only)")
    }
    if (sum(m$d$y) < 100 || isTRUE(m$engine_glm)) sp$`Firth logit (no random intercept)` <- \(m) add_info(fit_firth(m), m, "Firth logit (no random intercept)")
    imap(sp, \(f, nm) list(id = id, spec = nm, f = f))
  }))
  msg("Specifications: ", length(jobs), " fits")
  res <- pmap_safe(jobs, function(j) j$f(models[[j$id]]))
  sr <- bind_rows(res)
  write_csv(sr, file.path(outdir, "spec_refits.csv"))
  msg("  done: ", length(res), " of ", length(jobs))
  missing_jobs <- setdiff(paste(map_chr(jobs, "id"), map_chr(jobs, "spec")), unique(paste(sr$model, sr$spec)))
  if (length(missing_jobs)) msg("  FAILED: ", paste(missing_jobs, collapse = "; "))

  # Nonlinearity: add quadratic (and cubic) terms of one driver at a time
  poly_vars <- c("gdp_pc_within", "conflict_events_within", "conflict_deaths_within", "drought12_ext", "spei48_dec")
  jobs <- flatten(imap(models[!grepl("\\(change\\)", names(models))], function(m, id) {
    vs <- intersect(poly_vars, names(m$d))
    vs <- vs[vs %in% c(m$covs, names(m$slopes))]
    flatten(map(vs, \(v) map(2:3, \(dg) list(id = id, v = v, dg = dg))))
  }))
  res <- pmap_safe(jobs, function(j) {
    m <- models[[j$id]]; d <- m$d; extra <- character(0)
    for (k in 2:j$dg) { nm <- paste0(j$v, "_pow", k); d[[nm]] <- d[[j$v]]^k; extra <- c(extra, nm) }
    fit <- fit_m(m, d = d, covs = c(m$covs, extra))
    b <- b_vec(fit)
    b <- b[!is.na(b)]; V <- v_mat(fit)[names(b), names(b)]
    w <- wald_zero(b, V, extra)
    tibble(model = j$id, driver = j$v, degree = j$dg, chi2 = w[["chi2"]], df = w[["df"]], p = w[["p"]])
  })
  write_csv(bind_rows(res), file.path(outdir, "nonlinearity.csv"))

  # Variance inflation factors
  vif <- imap_dfr(models, function(m, id) {
    x <- m$d[usable(m$covs, m$d)]
    tibble(model = id, term = names(x), vif = diag(solve(cor(as.matrix(x)))))
  })
  write_csv(vif, file.path(outdir, "vif.csv"))

  # Multinomial logit over all roles (role 1 = baseline), as in the first draft
  mn <- bind_rows(
    tryCatch(fit_multinomial(cfg, mds$refugee) |> mutate(dataset = "refugee"), error = function(e) { msg("multinomial refugee: ", conditionMessage(e)); NULL }),
    tryCatch(fit_multinomial(cfg, mds$migrant) |> mutate(dataset = "migrant backbone"), error = function(e) { msg("multinomial migrant: ", conditionMessage(e)); NULL }))
  write_csv(mn, file.path(outdir, "multinomial.csv"))
}

# ---- C. alternative data -----------------------------------------------------
if ("data" %in% parts) {
  out <- list()
  # (1) alternative drought measures: results of the other config files, when they have been run
  alts <- c(output_original_spei01 = "SPEI-01 min/mean/CV (first draft)",
            output_spei_jordan = "SPEI-48 average minimum, within-country z",
            output_drought_interim = "interim drought measure",
            output_nocontrols_conflict = "no conflict controls",
            output_drought_agriculture = "drought x agriculture share")
  for (o in names(alts)) {
    p <- function(...) file.path(cfg$root, o, ...)
    if (file.exists(p("unhcr", "binary_role_models.csv")))
      out[[paste(o, "R")]] <- read_csv(p("unhcr", "binary_role_models.csv"), show_col_types = FALSE) |>
        transmute(model = paste0("R", role), term, estimate, se = std.error, p = p.value, n = n_obs, events = n_events, engine, version = alts[[o]], kind = "drought measure")
    if (file.exists(p("undesa_backbone", "binary_role_models.csv")))
      out[[paste(o, "M")]] <- read_csv(p("undesa_backbone", "binary_role_models.csv"), show_col_types = FALSE) |>
        transmute(model = paste0("M", role), term, estimate, se = std.error, p = p.value, n = n_obs, events = n_events, engine, version = alts[[o]], kind = "drought measure")
    if (file.exists(p("unhcr", "role_changes", "transition_models.csv")))
      out[[paste(o, "C")]] <- read_csv(p("unhcr", "role_changes", "transition_models.csv"), show_col_types = FALSE) |>
        transmute(model = paste0(change, ifelse(spec == "level in year t", " (level)", " (change)")), term, estimate, se = std.error, p = p.value, n, events, engine, version = alts[[o]], kind = "drought measure")
    if (!file.exists(p("unhcr", "binary_role_models.csv"))) msg("  not run yet: ", o)
  }

  # (2) other definitions of the roles
  drivers_r <- read_csv(out_path(cfg, "unhcr", "drivers_by_country_period.csv"), show_col_types = FALSE)
  drivers_m <- read_csv(out_path(cfg, "undesa", "drivers_by_country_period.csv"), show_col_types = FALSE)
  altdir <- file.path(outdir, "alt_roles")
  cfg_ref <- cfg; ov <- cfg$flows$unhcr$models %||% list(); cfg_ref$models[names(ov)] <- ov

  # refugees: minimum tie size of 1 (paper), 10 or 100 people
  flows_r <- read_flows(cfg, "unhcr")
  for (thr in c(10, 100)) {
    c2 <- cfg_ref; c2$flows$unhcr$roles <- list(thresholds = list(thr)); c2$output_dir <- file.path(cfg$output_dir, "robustness", "alt_roles", paste0("thr", thr))
    rr <- run_roles(c2, "unhcr", flows_r)$roles
    # role numbers can differ between runs: map through the labels
    lab <- rr |> distinct(role, label) |>
      mutate(code = case_when(grepl("sender only", label) ~ 1L, grepl("two-way", label) ~ 2L, grepl("receiver only", label) ~ 3L, TRUE ~ NA_integer_))
    msg("  refugee tie threshold ", thr, ": ", paste(lab$label, collapse = " | "))
    dup <- anyDuplicated(lab$code[!is.na(lab$code)]) > 0
    rr <- rr |> left_join(lab |> select(role, code), by = "role") |> mutate(role = code) |> select(-code) |> filter(!is.na(role))
    bin <- fit_binary_role_models(c2, prepare_model_data(c2, rr, drivers_r))
    out[[paste("thr", thr)]] <- bin |> transmute(model = paste0("R", role), term, estimate, se = std.error, p = p.value, n = n_obs, events = n_events,
                                                 engine, version = paste0("refugee ties of at least ", thr, if (dup) " (labels not unique)" else ""), kind = "role definition")
  }

  # migrants: backbone share 5%, 20% and the other Abel & Cohen estimate
  label_bb <- function(roles_df) {
    p <- roles_df |> group_by(role) |> summarise(s = mean(in_degree == 0), r = mean(out_degree == 0), .groups = "drop")
    setNames(case_when(p$s > 0.5 ~ 1L, p$r > 0.5 ~ 3L, TRUE ~ 2L), p$role)
  }
  rc <- cfg$roles; rc$k <- 3
  for (v in list(list("da_pb_closed", 0.05, "pseudo-Bayesian, 5% backbone"), list("da_pb_closed", 0.20, "pseudo-Bayesian, 20% backbone"),
                 list("da_min_closed", 0.10, "minimisation, 10% backbone"))) {
    c2 <- cfg; c2$flows$undesa$value <- v[[1]]; c2$flows$undesa$net <- TRUE
    fl <- read_flows(c2, "undesa")
    bb <- fl |> group_by(origin, period) |> filter(flow >= v[[2]] * sum(flow)) |> ungroup()
    sim <- rege(build_stacked_network(bb, 1), rc$rege_iterations, verbose = FALSE)$similarity
    cl <- cluster_roles(1 - sim, rc)
    k3 <- cutree(cl$hclust, 3)
    roles <- tibble(node = names(k3), role = unname(k3)) |> separate(node, c("iso3", "period"), sep = "_", convert = TRUE) |>
      left_join(node_metrics(bb, 1) |> select(iso3, period, in_degree, out_degree), by = c("iso3", "period"))
    code <- label_bb(roles)
    dup <- anyDuplicated(code) > 0
    # if two clusters get the same position, keep the larger of them for that position
    sizes <- table(roles$role)
    roles <- roles |> mutate(role_new = code[as.character(role)], size = as.integer(sizes[as.character(role)])) |>
      group_by(role_new) |> filter(!dup | size == max(size)) |> ungroup() |> mutate(role = role_new) |> select(-role_new, -size)
    roles <- roles |> left_join(read_csv(in_root(cfg, cfg$country_list), show_col_types = FALSE) |>
                                  select(iso3 = iso3c, country = name, region = un_region), by = "iso3") |> mutate(label = as.character(role))
    bin <- fit_binary_role_models(cfg, prepare_model_data(cfg, roles, drivers_m))
    out[[v[[3]]]] <- bin |> transmute(model = paste0("M", role), term, estimate, se = std.error, p = p.value, n = n_obs, events = n_events,
                                      engine, version = paste0(v[[3]], if (dup) " (two clusters in one position: larger kept)" else ""), kind = "role definition")
  }
  write_csv(bind_rows(out), file.path(outdir, "data_refits.csv"))
}

# ---- D. inference ------------------------------------------------------------
# Spatial-serial HAC (Conley 1999; Hsiang 2010), for a logistic regression.
#   s_it   score of row (country i, period t)
#   Meat = sum_t sum_ij K(dist_ij) s_it s_jt'                      (uniform kernel within `cutoff` km)
#        + sum_i sum_{l=1..L} (1 - l/(L+1)) (s_it s_i,t-l' + s_i,t-l s_it')   (Bartlett over `L` lags)
#   Var  = Bread Meat Bread
vcov_conley <- function(fit, d, D, cutoff, L) {
  s <- sandwich::estfun(fit); br <- summary(fit)$cov.unscaled
  k <- ncol(s); M <- matrix(0, k, k)
  idx <- match(d$iso3, rownames(D))
  for (t in unique(d$period)) {
    i <- which(d$period == t); w <- D[idx[i], idx[i], drop = FALSE] <= cutoff
    w[is.na(w)] <- FALSE; diag(w) <- TRUE
    M <- M + crossprod(s[i, , drop = FALSE], w %*% s[i, , drop = FALSE])
  }
  step <- min(diff(sort(unique(d$period))))
  for (u in unique(d$iso3)) {
    i <- which(d$iso3 == u); if (length(i) < 2) next
    for (a in seq_along(i)) for (b in seq_along(i)) {
      l <- (d$period[i[a]] - d$period[i[b]]) / step
      if (l > 0 && l <= L) {
        g <- tcrossprod(s[i[a], ], s[i[b], ]) * (1 - l / (L + 1))
        M <- M + g + t(g)
      }
    }
  }
  br %*% M %*% br
}

if ("inference" %in% parts) {
  # country coordinates (Natural Earth label points; successor state for historical codes)
  shp <- in_root(cfg, "data/raw/ne_10m_admin_0_countries/ne_10m_admin_0_countries.dbf")
  geo <- foreign::read.dbf(shp, as.is = TRUE) |>
    transmute(iso3 = ifelse(ISO_A3_EH == "-99", ADM0_A3, ISO_A3_EH), lon = LABEL_X, lat = LABEL_Y) |> distinct(iso3, .keep_all = TRUE)
  geo <- bind_rows(geo, geo |> filter(iso3 == "SDN") |> mutate(iso3 = "SUD"), geo |> filter(iso3 == "SRB") |> mutate(iso3 = "SCG"))
  isos <- sort(unique(unlist(map(models, \(m) m$d$iso3))))
  miss <- setdiff(isos, geo$iso3); if (length(miss)) msg("  no coordinates for: ", paste(miss, collapse = ", "))
  g <- geo |> filter(iso3 %in% isos)
  rad <- pi / 180
  D <- outer(seq_len(nrow(g)), seq_len(nrow(g)), Vectorize(function(i, j) {
    a <- sin((g$lat[j] - g$lat[i]) * rad / 2)^2 + cos(g$lat[i] * rad) * cos(g$lat[j] * rad) * sin((g$lon[j] - g$lon[i]) * rad / 2)^2
    6371 * 2 * asin(pmin(1, sqrt(a)))
  })); dimnames(D) <- list(g$iso3, g$iso3)

  inf <- imap_dfr(models, function(m, id) {
    d <- m$d; d$period_f <- droplevels(factor(d$period)); covs <- usable(m$covs, d)
    fit <- glm(as.formula(paste("y ~", paste(c(covs, time_terms(m$time)), collapse = " + "))), data = d, family = binomial)
    b <- coef(fit); ok <- !is.na(b); terms <- names(b)[ok & !grepl("Intercept|^period_", names(b))]
    L <- if (m$step == 1) 5 else 1
    vs <- list(`logit, naive` = vcov(fit),
               `clustered by country` = sandwich::vcovCL(fit, cluster = d$iso3, type = "HC1"),
               `clustered by subregion` = sandwich::vcovCL(fit, cluster = d$subregion, type = "HC1"),
               `two-way: country and period` = sandwich::vcovCL(fit, cluster = d[c("iso3", "period")], type = "HC1"),
               `Conley 1,000 km` = vcov_conley(fit, d, D, 1000, L),
               `Conley 2,000 km` = vcov_conley(fit, d, D, 2000, L))
    imap_dfr(vs, function(V, nm) {
      se <- sqrt(pmax(diag(V)[terms], 0))
      tibble(model = id, term = terms, estimate = b[terms], se = se, p = 2 * pnorm(-abs(b[terms] / se)), inference = nm)
    })
  })
  write_csv(inf, file.path(outdir, "inference.csv"))
}

# ---- E. placebo --------------------------------------------------------------
# Studentized permutation test: compare the observed t statistic with its
# distribution when the drivers carry no information about the outcome.
#
#  (1) "country reassignment" (for all models). Each country is given the whole
#      driver history of another country drawn at random: the drivers keep their
#      correlations with each other and their persistence over time, while the
#      link with the outcome is broken. Donors come from the full driver panel
#      (all countries with complete drivers in every period of the model), so
#      no rows are lost. Draws in which the model does not converge are discarded
#      and counted.
#  (2) "shuffled years within country" (linear probability model with country
#      fixed effects). Keeps each country's own drivers but randomizes their
#      timing: this tests the within-country timing of the drivers.
#
# The observed t statistic is computed with the same estimator as the placebo runs.
# the driver columns that are permuted in the placebo runs (all usable covariates)
xcols <- function(m) unique(c(usable(m$covs, m$d)))
# Fit the mixed model in glmmTMB and return term, log-odds estimate and z statistic.
# Used by the placebo and bootstrap runs. An error (no convergence, non-positive-
# definite Hessian, or a non-finite z) discards that run through pmap_safe().
glmm_fast <- function(m, d, covs) {
  d$period_f <- droplevels(factor(d$period)); covs <- usable(covs, d)
  rhs <- c(covs, time_terms(m$time))
  if (isTRUE(m$engine_glm)) {
    co <- summary(glm(as.formula(paste("y ~", paste(rhs, collapse = " + "))), data = d, family = binomial))$coefficients
  } else {
    fit <- suppressMessages(suppressWarnings(glmmTMB::glmmTMB(
      as.formula(paste("y ~", paste(rhs, collapse = " + "), "+ (1 | iso3)")), data = d, family = binomial,
      control = glmmTMB::glmmTMBControl(parallel = 1))))
    if (fit$fit$convergence != 0 || !isTRUE(fit$sdr$pdHess)) stop("not converged")   # draw is discarded
    co <- summary(fit)$coefficients$cond
  }
  if (any(!is.finite(co[, 3]))) stop("not converged")
  co <- co[!grepl("Intercept|^period_", rownames(co)), , drop = FALSE]
  tibble(term = rownames(co), est = co[, 1], z = co[, 3])
}
# Full driver panel, transformed like the model data. The role column is a
# placeholder: only the drivers are used. A t statistic does not depend on the
# scale of a covariate, so standardizing over this larger panel is harmless.
panel_pool <- function(c0, drivers) {
  grid <- drivers |> distinct(iso3, period) |> mutate(country = iso3, region = "x", role = 1L, label = "x")
  prepare_model_data(c0, grid, drivers)$data |> mutate(iso3 = as.character(iso3))
}
# adds to m: P (panel of donors) and donors (countries complete in every period of the model)
add_donors <- function(m, P) {
  X <- xcols(m)
  dv <- X[grepl("^d_", X)]
  if (length(dv)) P <- P |> arrange(iso3, period) |> group_by(iso3) |>
    mutate(across(all_of(sub("^d_", "", dv)), \(x) ifelse(period - lag(period) == 1, x - lag(x), NA), .names = "d_{.col}")) |> ungroup()
  P <- P[stats::complete.cases(P[X]), ]
  per <- unique(m$d$period)
  cnt <- P |> filter(period %in% per) |> count(iso3)
  m$P <- P; m$donors <- cnt$iso3[cnt$n == length(per)]
  m
}
# one draw of test (1): each country gets the complete driver history of a donor
permute_countries <- function(m) {
  X <- xcols(m)
  ctry <- unique(m$d$iso3)
  donor <- setNames(sample(m$donors, length(ctry), replace = length(m$donors) < length(ctry)), ctry)
  i <- match(paste(donor[m$d$iso3], m$d$period), paste(m$P$iso3, m$P$period))
  d <- m$d; d[X] <- m$P[i, X]
  d
}
# one draw of test (2): the rows of drivers are shuffled within each country
permute_within <- function(m) {
  X <- xcols(m); d <- m$d
  idx <- ave(seq_len(nrow(d)), d$iso3, FUN = function(x) if (length(x) > 1) sample(x) else x)
  d[X] <- d[idx, X]
  d
}

if ("placebo" %in% parts) {
  draws <- list()
  cfg_r <- cfg; ov <- cfg$flows$unhcr$models %||% list(); cfg_r$models[names(ov)] <- ov
  P_r <- panel_pool(cfg_r, read_csv(out_path(cfg, "unhcr", "drivers_by_country_period.csv"), show_col_types = FALSE))
  P_m <- panel_pool(cfg, read_csv(out_path(cfg, "undesa", "drivers_by_country_period.csv"), show_col_types = FALSE))
  for (id in names(models)) {
    m <- add_donors(models[[id]], if (grepl("^M", id)) P_m else P_r)
    msg("  ", id, ": ", length(m$donors), " donor countries for ", length(unique(m$d$iso3)), " countries")
    t0 <- Sys.time()
    obs <- glmm_fast(m, m$d, m$covs)
    res <- pmap_safe(as.list(seq_len(nperm)), function(k) {
      set.seed(cfg$seed + k)
      d <- permute_countries(m)
      glmm_fast(m, d, xcols(m)) |> mutate(draw = k, n_used = nrow(d))
    })
    a <- bind_rows(res) |> mutate(model = id, test = "country reassignment")
    obs_l <- fit_lpm(m) |> transmute(term, est = estimate, z = estimate / se)
    res2 <- pmap_safe(as.list(seq_len(nperm)), function(k) {
      set.seed(cfg$seed + 100000 + k)
      fit_lpm(m, d = permute_within(m)) |> transmute(term, est = estimate, z = estimate / se, draw = k, n_used = nrow(m$d))
    })
    b <- bind_rows(res2) |> mutate(model = id, test = "shuffled years within country")
    draws[[id]] <- list(draws = bind_rows(a, b),
                        obs = bind_rows(obs |> mutate(test = "country reassignment"), obs_l |> mutate(test = "shuffled years within country")) |> mutate(model = id))
    msg("  placebo ", id, ": ", length(res), " + ", length(res2), " draws in ",
        round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min; mean rows used ", round(mean(a$n_used)), " of ", nrow(m$d))
  }
  all_d <- bind_rows(map(draws, "draws")); all_o <- bind_rows(map(draws, "obs"))
  all_d <- all_d |> filter(is.finite(z), abs(z) < 50)              # separated fits are not informative
  write_csv(all_d |> mutate(across(c(est, z), \(x) signif(x, 4))), file.path(outdir, "placebo_draws.csv"))
  ps <- all_d |> inner_join(all_o |> rename(obs_est = est, obs_z = z), by = c("model", "term", "test")) |>
    group_by(model, term, test, obs_est, obs_z) |>
    summarise(n_draws = n(), mean_placebo_z = mean(z), sd_placebo_z = sd(z),
              # share of placebo runs with a t statistic at least as extreme as the observed one
              p_perm = (1 + sum(abs(z) >= abs(obs_z[1]))) / (n() + 1),
              rate_reject_05 = mean(abs(z) > 1.96), .groups = "drop")
  write_csv(ps, file.path(outdir, "placebo_summary.csv"))
}
# ---- D2. cluster bootstrap of the mixed model ---------------------------------
# The standard errors in part D belong to a pooled logistic regression, which is
# not the model in the paper (it has no country intercept, so its coefficients
# differ). This part keeps the paper's model: resample whole countries with
# replacement (each resampled copy gets its own random intercept), refit, and
# read the uncertainty from the spread of the estimates.
#
# Why it matters. The random intercept treats a country's successive years as
# independent once its average is accounted for. Roles and drivers persist from
# year to year, so the model-based standard errors of the within-country terms are
# too small: under the country-reassignment placebo, which keeps that persistence,
# these terms exceed |t| = 1.96 in about 20% of runs instead of 5%. Resampling whole
# countries keeps each country's time series intact and needs no such assumption.
#
# Results (bootstrap.csv, log-odds scale):
#   boot_se                standard deviation of the refitted estimates
#   boot_low, boot_high    2.5% and 97.5% percentiles: the 95% intervals in Figs. 3
#                          and 4 and in the Results, around the main model's estimate
#   n_boot                 refits used for this term (those that converged)
# The standardization and the within/between terms are those of the full sample.
# Caution: for rare EM-DAT disaster terms a few events drive the estimate, and the
# percentile interval can be narrower than the model-based one or lopsided.
if ("bootstrap" %in% parts) {
  nboot <- as.integer(Sys.getenv("ROBUST_NBOOT", "500"))
  boot <- imap_dfr(models, function(m, id) {
    idx <- split(seq_len(nrow(m$d)), m$d$iso3)
    res <- pmap_safe(as.list(seq_len(nboot)), function(k) {
      set.seed(cfg$seed + 200000 + k)
      pick <- sample(names(idx), replace = TRUE)
      rows <- unlist(idx[pick], use.names = FALSE)
      d <- m$d[rows, ]
      d$iso3 <- rep(paste0("b", seq_along(pick)), lengths(idx[pick]))    # a new id for each copy
      glmm_fast(m, d, m$covs) |> mutate(draw = k)
    })
    msg("  bootstrap ", id, ": ", length(res), " of ", nboot, " refits used")
    # |log odds| >= 20 means separation in that resample: no usable estimate of the term
    bind_rows(res) |> filter(is.finite(est), abs(est) < 20) |> mutate(model = id)
  })
  bs <- boot |> group_by(model, term) |>
    summarise(n_boot = n(), boot_se = sd(est), boot_low = quantile(est, 0.025), boot_high = quantile(est, 0.975), .groups = "drop") |>
    inner_join(baseline |> select(model, term, b0 = estimate), by = c("model", "term")) |>
    # z and p_boot use the bootstrap SE and a normal approximation (Fig. S14)
    mutate(z = b0 / boot_se, p_boot = 2 * pnorm(-abs(z)), percentile_excludes_0 = boot_low > 0 | boot_high < 0)
  write_csv(bs, file.path(outdir, "bootstrap.csv"))
}
msg("Done: ", outdir)
