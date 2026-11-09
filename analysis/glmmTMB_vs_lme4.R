# =============================================================================
# glmmTMB_vs_lme4.R — do glmmTMB and lme4 give the same main models?
#
#   Rscript analysis/glmmTMB_vs_lme4.R [config.yml]
#   (after regression_robustness.R, which writes robustness/baseline.csv; a few minutes)
#
# 1. Refits every main model with lme4 through fit_glmer_robust() and compares the
#    coefficients with glmmTMB.
# 2. For the models where lme4 from its default start does not converge cleanly,
#    restarts lme4 from the glmmTMB estimates and compares the log-likelihoods.
#
# Writes to <output_dir>/robustness/: glmmTMB_vs_lme4.csv, lme4_restart_check.csv
# =============================================================================
src <- readLines("analysis/regression_robustness.R")
eval(parse(text = src[seq_len(grep("^# ---- baseline", src)[1] - 1)]))
built <- build_models(); models <- built$models
base <- read_csv(file.path(outdir, "baseline.csv"), show_col_types = FALSE)
res <- bind_rows(parallel::mclapply(names(models), function(id) {
  m <- models[[id]]; if (grepl("logistic", unique(base$engine[base$model == id]))) return(NULL)
  f <- fit_m(m, fast = FALSE); t <- tidy_m(f, m) |> filter(!grepl("^wald|^joint", term))
  x <- t |> inner_join(base |> filter(model == id) |> select(term, b = estimate, s = se), by = "term")
  tibble(model = id, eng = unique(t$engine), maxb = max(abs(x$estimate - x$b)), maxse = max(abs(x$se - x$s)), maxz = max(abs(x$estimate / x$se - x$b / x$s)))
}, mc.cores = 2, mc.preschedule = FALSE))
write_csv(res, file.path(outdir, "glmmTMB_vs_lme4.csv"))
ids <- res$model[grepl("logistic", res$eng)]
chk <- map_dfr(ids, function(id) {
  m <- models[[id]]; d <- m$d; d$period_f <- droplevels(factor(d$period)); covs <- usable(m$covs, d)
  f <- as.formula(paste("y ~", paste(c(covs, time_terms(m$time)), collapse = " + "), "+ (1 | iso3)"))
  ctl <- lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
  t <- glmmTMB::glmmTMB(f, data = d, family = binomial); b <- glmmTMB::fixef(t)$cond
  g <- suppressWarnings(lme4::glmer(f, data = d, family = binomial, control = ctl))
  g2 <- suppressWarnings(lme4::glmer(f, data = d, family = binomial, control = ctl,
                                     start = list(theta = sqrt(as.numeric(glmmTMB::VarCorr(t)$cond$iso3)), fixef = b)))
  tibble(model = id, ll_glmmTMB = as.numeric(logLik(t)), ll_lme4 = as.numeric(logLik(g)), ll_lme4_restart = as.numeric(logLik(g2)),
         maxb = max(abs(lme4::fixef(g)[names(b)] - b)), maxb_restart = max(abs(lme4::fixef(g2)[names(b)] - b)),
         messages = paste(g@optinfo$conv$lme4$messages, collapse = " | "))
})
write_csv(chk, file.path(outdir, "lme4_restart_check.csv"))
print(as.data.frame(res)); print(as.data.frame(chk))
