# =============================================================================
# regression_robustness_summary.R — combine the checks and draw the figures
#
#   Rscript analysis/regression_robustness_summary.R [config.yml]
#   (after analysis/regression_robustness.R)
#
# Reads the CSV files in <output_dir>/robustness/ and writes there:
#   summary_table.csv    one row per model and driver: every check
#   sample_summary.csv   range of each estimate over each type of exclusion
#   Fig_robust_sample.png, Fig_robust_checks.png, Fig_robust_placebo.png
#
# No pass/fail classification is made: each check is reported as an estimate and
# its uncertainty. analysis/regression_evidence.R brings the results together with
# the Bayesian multilevel models.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(readr); library(ggplot2); library(stringr)
})
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
cfg_file <- commandArgs(trailingOnly = TRUE)
cfg <- read_config(if (length(cfg_file)) cfg_file[1] else "config.yml")
outdir <- dirname(out_path(cfg, "robustness", "x"))
rd <- function(f) read_csv(file.path(outdir, f), show_col_types = FALSE)

base <- rd("baseline.csv")
sample_r <- rd("sample_refits.csv")
spec_r <- rd("spec_refits.csv")
data_r <- rd("data_refits.csv")
inf_r <- rd("inference.csv")
plc <- rd("placebo_summary.csv")
draws <- rd("placebo_draws.csv")
bs <- rd("bootstrap.csv")

# rows of the results tables that hold test statistics, not coefficients
skip_term <- function(t) grepl("^wald_equal|^joint_income", t)
b0 <- base |> filter(!skip_term(term)) |> select(model, group, term, b0 = estimate, se0 = se, p0 = p, n, events, engine)

# ---- sample composition ------------------------------------------------------
# for each model, driver and type of exclusion: range of the odds ratio, share of
# refits with the baseline sign, share with p < 0.05, and the unit whose removal
# moves the estimate most
sample_sum <- sample_r |> filter(!skip_term(term)) |> inner_join(b0 |> select(model, term, b0), by = c("model", "term")) |>
  group_by(model, term, type, b0) |>
  summarise(refits = n(), or_min = exp(min(estimate)), or_max = exp(max(estimate)),
            share_same_sign = mean(sign(estimate) == sign(b0[1])), share_p05 = mean(p < 0.05),
            most_influential = unit[which.max(abs(estimate - b0[1]))],
            or_when_removed = exp(estimate[which.max(abs(estimate - b0[1]))]),
            refits_logistic = sum(engine != "mixed (glmmTMB)" & !grepl("^mixed", engine)), .groups = "drop") |>
  select(-b0)
write_csv(sample_sum, file.path(outdir, "sample_summary.csv"))

# ---- specification -----------------------------------------------------------
# specifications on the log-odds scale; the linear probability and conditional
# logit models have other scales or samples and are summarised separately below
logit_specs <- c("adaptive quadrature (10 points)", "decade effects", "linear time trend", "drivers lagged one period",
                 "no disaster deaths or affected", "no conflict variables", "no SPEI drought variables",
                 "no current GDP (terciles only)", "one drought slope, no income terciles", "Firth logit (no random intercept)")
spec_sum <- spec_r |> filter(spec %in% logit_specs, !skip_term(term)) |> inner_join(b0 |> select(model, term, b0), by = c("model", "term")) |>
  group_by(model, term) |>
  summarise(spec_n = n(), spec_or_min = exp(min(estimate)), spec_or_max = exp(max(estimate)),
            spec_same = mean(sign(estimate) == sign(b0)), spec_p05 = mean(p < 0.05), .groups = "drop")
within_sum <- spec_r |> filter(grepl("conditional logit|linear probability", spec), !skip_term(term)) |>
  mutate(k = ifelse(grepl("conditional", spec), "clogit", "lpm")) |>
  select(model, term, k, estimate, p) |>
  pivot_wider(names_from = k, values_from = c(estimate, p), names_glue = "{k}_{.value}")

# ---- alternative data --------------------------------------------------------
data_sum <- data_r |> filter(!skip_term(term), !grepl("Intercept|period", term)) |>
  inner_join(b0 |> select(model, term, b0), by = c("model", "term")) |>
  group_by(model, term) |>
  summarise(data_n = n(), data_or_min = exp(min(estimate)), data_or_max = exp(max(estimate)),
            data_same = mean(sign(estimate) == sign(b0)), data_p05 = mean(p < 0.05), .groups = "drop")

# ---- inference ---------------------------------------------------------------
inf_wide <- inf_r |> select(model, term, inference, p) |>
  pivot_wider(names_from = inference, values_from = p, names_prefix = "p ")
inf_need <- c("p clustered by country", "p two-way: country and period", "p Conley 1,000 km", "p Conley 2,000 km")

# ---- combine -----------------------------------------------------------------
wide <- function(type_) sample_sum |> filter(type == type_) |> select(-type)
s_c <- wide("country") |> rename_with(\(x) paste0("loco_", x), -c(model, term))
s_other <- sample_sum |> filter(type %in% c("subregion", "continent", "time block", "small states", "year 2020")) |>
  group_by(model, term) |> summarise(other_same = weighted.mean(share_same_sign, refits), .groups = "drop")
s_block <- wide("time block") |> transmute(model, term, block_or_min = or_min, block_or_max = or_max, block_same = share_same_sign, block_p05 = share_p05)
s_region <- wide("subregion") |> transmute(model, term, region_or_min = or_min, region_or_max = or_max, region_same = share_same_sign, region_p05 = share_p05,
                                           region_most = most_influential)
s_event <- wide("event") |> transmute(model, term, event_same = share_same_sign, event_p05 = share_p05)
pl <- plc |> filter(test == "country reassignment") |> transmute(model, term, p_placebo = p_perm)
pl2 <- plc |> filter(test == "shuffled years within country") |> transmute(model, term, p_placebo_within = p_perm)

tab <- b0 |> left_join(s_c, by = c("model", "term")) |> left_join(s_other, by = c("model", "term")) |>
  left_join(s_block, by = c("model", "term")) |> left_join(s_region, by = c("model", "term")) |>
  left_join(s_event, by = c("model", "term")) |> left_join(spec_sum, by = c("model", "term")) |>
  left_join(within_sum, by = c("model", "term")) |> left_join(data_sum, by = c("model", "term")) |>
  left_join(inf_wide, by = c("model", "term")) |> left_join(bs |> select(model, term, n_boot, boot_se, boot_low, boot_high, p_boot), by = c("model", "term")) |>
  left_join(pl, by = c("model", "term")) |> left_join(pl2, by = c("model", "term")) |>
  mutate(or0 = exp(b0), ci_low = exp(b0 - 1.96 * se0), ci_high = exp(b0 + 1.96 * se0))

# Benjamini-Hochberg across all drivers in the role and role-change models
tab$fdr_q <- p.adjust(tab$p0, method = "BH")
write_csv(tab, file.path(outdir, "summary_table.csv"))

# ---- figures -----------------------------------------------------------------
model_lab <- c(R1 = "R1 sender only", R2 = "R2 two-way", R3 = "R3 receiver only", M1 = "M1 source", M2 = "M2 intermediate",
               `start sending (level)` = "Start sending (level)", `start sending (change)` = "Start sending (change)",
               `stop sending (level)` = "Stop sending (level)", `stop sending (change)` = "Stop sending (change)",
               `start hosting (level)` = "Start hosting (level)", `start hosting (change)` = "Start hosting (change)")
term_lab <- function(t) {
  r <- c(gdp_pc_within = "GDP per capita, within", gdp_pc_between = "GDP per capita, between",
         conflict_events_within = "Conflict events, within", conflict_events_between = "Conflict events, between",
         conflict_deaths_within = "Conflict deaths, within", conflict_deaths_between = "Conflict deaths, between",
         gdp_pc = "GDP per capita", income_g2 = "Income tercile: middle", income_g3 = "Income tercile: high",
         conflict_events = "Conflict events", conflict_deaths = "Conflict deaths",
         drought12_ext = "Drought extent", d_drought12_ext = "Rise in drought extent",
         d_conflict_events = "Rise in conflict events", d_conflict_deaths = "Rise in conflict deaths", d_gdp_pc = "Change in GDP",
         drought12_ext_x_g1 = "Drought extent, low income", drought12_ext_x_g2 = "Drought extent, middle income",
         drought12_ext_x_g3 = "Drought extent, high income",
         drought12_ext_lag1_x_g1 = "Drought extent t-1, low income", drought12_ext_lag1_x_g2 = "Drought extent t-1, middle income",
         drought12_ext_lag1_x_g3 = "Drought extent t-1, high income",
         spei48_dec_x_g1 = "SPEI-48, low income", spei48_dec_x_g2 = "SPEI-48, middle income", spei48_dec_x_g3 = "SPEI-48, high income")
  out <- r[t]; out[is.na(out)] <- t[is.na(out)]; unname(out)
}
# The associations stated in the Results: these are the rows of Figs. S13-S15,
# and the same list defines the rows of Table S9 in regression_evidence.R.
claims <- bind_rows(
  tidyr::expand_grid(model = c("R1", "R2", "R3", "M1", "M2"),
                     term = c("conflict_events_within", "conflict_events_between", "conflict_deaths_within", "conflict_deaths_between",
                              "gdp_pc_within", "gdp_pc_between")),
  tribble(~model, ~term,
    "start sending (level)", "conflict_deaths_within", "start sending (level)", "conflict_deaths_between", "start sending (level)", "drought12_ext_x_g1",
    "start sending (change)", "d_conflict_deaths", "start sending (change)", "d_drought12_ext",
    "stop sending (level)", "conflict_deaths_between",
    "start hosting (level)", "conflict_deaths_within", "start hosting (level)", "conflict_deaths_between", "start hosting (level)", "income_g3"))
ml <- function(m) factor(model_lab[m], levels = unname(model_lab))
theme_set(theme_minimal(9) + theme(plot.title.position = "plot", legend.position = "bottom", panel.grid.minor = element_blank()))
save_fig <- function(p, file, w, h) ggsave(file.path(outdir, file), p, width = w, height = h, dpi = 300, bg = "white")

# Fig A: odds ratio with its 95% CI, and the range over each type of exclusion
A <- claims |> inner_join(tab, by = c("model", "term")) |>
  transmute(model, term, or0, ci_low, ci_high, label = term_lab(term), m = ml(model))
rng <- claims |> inner_join(sample_sum |> filter(type %in% c("country", "subregion", "time block")), by = c("model", "term")) |>
  mutate(check = recode(type, country = "Leave out a country", subregion = "Leave out a subregion", `time block` = "Leave out a time block"),
         check = factor(check, c("Leave out a country", "Leave out a subregion", "Leave out a time block"))) |>
  left_join(A |> select(model, term, label, m), by = c("model", "term"))
pA <- ggplot(A, aes(y = label)) +
  geom_vline(xintercept = 1, colour = "grey60") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0, colour = "grey75", linewidth = 2.6) +
  geom_linerange(data = rng, aes(xmin = or_min, xmax = or_max, colour = check), linewidth = 0.9,
                 position = position_dodge(width = 0.55)) +
  geom_point(aes(x = or0), size = 1.6, colour = "black") +
  scale_x_log10(breaks = c(0.25, 1, 4, 16), labels = c("1/4", "1", "4", "16")) +
  scale_colour_manual(values = c("Leave out a country" = "#2f6fb5", "Leave out a subregion" = "#c0504a", "Leave out a time block" = "#1aa08a")) +
  facet_wrap(~m, scales = "free_y", ncol = 3) +
  labs(x = "Odds ratio per SD (log scale). Black point: full sample; grey bar: its 95% CI; coloured lines: range of estimates", y = NULL, colour = NULL)
save_fig(pA, "Fig_robust_sample.png", 9.5, 9)

# Fig B: z statistic of each association under every check. Each cell is the
# estimate divided by the standard error of that analysis; the bootstrap column
# uses the bootstrap SE and the placebo columns the observed t statistic.
zb <- base |> filter(!skip_term(term)) |> transmute(model, term, check = "Baseline (mixed model)", z = estimate / se, sig = p < 0.05, cat = "Baseline")
zs <- spec_r |> filter(!skip_term(term)) |> transmute(model, term, check = spec, z = estimate / se, sig = p < 0.05, cat = "Specification")
zi <- inf_r |> filter(inference != "logit, naive") |> transmute(model, term, check = inference, z = estimate / se, sig = p < 0.05, cat = "Inference (pooled logit)")
zbs <- bs |> transmute(model, term, check = "Country bootstrap", z = z, sig = p_boot < 0.05, cat = "Inference (mixed model)")
zp <- plc |> transmute(model, term, check = paste0("Placebo: ", test), z = obs_z, sig = p_perm < 0.05, cat = "Placebo")
Z <- bind_rows(zb, zs, zbs, zi, zp) |> inner_join(claims, by = c("model", "term")) |>
  mutate(rowlab = paste0(model_lab[model], ": ", term_lab(term)),
         rowlab = factor(rowlab, rev(unique(paste0(model_lab[claims$model], ": ", term_lab(claims$term))))),
         cat = recode(cat, Baseline = "Base", `Inference (mixed model)` = "Mixed", `Inference (pooled logit)` = "Pooled logit inference"),
         cat = factor(cat, c("Base", "Specification", "Mixed", "Pooled logit inference", "Placebo")),
         check = recode(check, `adaptive quadrature (10 points)` = "Adaptive quadrature", `no disaster deaths or affected` = "No disaster variables",
                        `no SPEI drought variables` = "No SPEI drought", `no conflict variables` = "No conflict variables",
                        `one drought slope, no income terciles` = "One drought slope", `no current GDP (terciles only)` = "No current GDP",
                        `drivers lagged one period` = "Drivers lagged", `decade effects` = "Decade effects", `linear time trend` = "Linear trend",
                        `Firth logit (no random intercept)` = "Firth logit", `conditional logit (within country)` = "Conditional logit",
                        `linear probability, country fixed effects` = "Linear prob., country FE",
                        `clustered by country` = "Cluster: country", `clustered by subregion` = "Cluster: subregion",
                        `two-way: country and period` = "Two-way cluster",
                        `Placebo: country reassignment` = "Country reassignment", `Placebo: shuffled years within country` = "Years within country"))
order_chk <- c("Baseline (mixed model)", "Adaptive quadrature", "Decade effects", "Linear trend", "Drivers lagged", "No disaster variables",
               "No conflict variables", "No SPEI drought", "No current GDP", "One drought slope", "Firth logit", "Conditional logit",
               "Linear prob., country FE", "Country bootstrap", "Cluster: country", "Cluster: subregion", "Two-way cluster", "Conley 1,000 km", "Conley 2,000 km",
               "Country reassignment", "Years within country")
Z <- Z |> mutate(check = factor(check, order_chk))
pB <- ggplot(Z, aes(check, rowlab, fill = pmax(pmin(z, 4), -4))) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.1f", z)), size = 1.8, colour = "black") +
  scale_fill_gradient2(low = "#2f6fb5", mid = "#f2f0ec", high = "#c0504a", midpoint = 0, limits = c(-4, 4), name = "z statistic (capped at 4)") +
  facet_grid(~cat, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = NULL, caption = "Number in each cell: z statistic. Blank: the term is not in that model.") +
  theme(axis.text.x = element_text(angle = 50, hjust = 1), panel.grid = element_blank(), strip.text = element_text(face = "bold"))
save_fig(pB, "Fig_robust_checks.png", 9.5, 9)

# Fig C: placebo distribution of the t statistic, with the observed value
# (associations with an observed |t| of at least 1.96 are shown; all are in placebo_summary.csv)
# From here on `claims` is the reduced list of associations shown in this figure.
shown <- plc |> filter(test == "country reassignment", abs(obs_z) >= 1.96) |> inner_join(claims, by = c("model", "term")) |> select(model, term)
claims_all <- claims; claims <- shown |> semi_join(claims, by = c("model", "term"))
D <- draws |> filter(test == "country reassignment") |> inner_join(claims, by = c("model", "term")) |>
  mutate(lab = paste0(model_lab[model], "\n", term_lab(term)))
O <- plc |> filter(test == "country reassignment") |> inner_join(claims, by = c("model", "term")) |>
  mutate(lab = paste0(model_lab[model], "\n", term_lab(term)))
lv <- unique(paste0(model_lab[claims$model], "\n", term_lab(claims$term)))
pC <- ggplot(D |> mutate(lab = factor(lab, lv)), aes(z)) +
  geom_histogram(binwidth = 0.25, fill = "#8aa4c1", colour = "white", linewidth = 0.1) +
  geom_vline(data = O |> mutate(lab = factor(lab, lv)), aes(xintercept = obs_z), colour = "#c0504a", linewidth = 0.7) +
  geom_text(data = O |> mutate(lab = factor(lab, lv)), aes(x = Inf, y = Inf, label = sprintf("p = %.3f", p_perm)),
            hjust = 1.1, vjust = 1.4, size = 2.3) +
  coord_cartesian(xlim = c(-6, 6)) +
  facet_wrap(~lab, ncol = 5, scales = "free_y") +
  labs(x = "t statistic in placebo runs (bars) and in the data (red line)", y = "Placebo runs") +
  theme(strip.text = element_text(size = 6.8))
save_fig(pC, "Fig_robust_placebo.png", 9.5, 6.5)

msg("Done: ", outdir)
