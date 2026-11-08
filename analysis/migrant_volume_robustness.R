# =============================================================================
# migrant_volume_robustness.R — do the migrant groups hold under other cut-offs?
#
#   Rscript analysis/migrant_volume_robustness.R      (~15 min; after run_pipeline.R)
#
# The migrant roles found by REGE track migration volume (in + out flow). This
# script asks how robust that grouping is, in three ways:
#   A. Volume tiers: group country-periods directly by total volume, using
#      several different cut-off rules (fixed breaks, quantiles, natural breaks).
#   B. REGE roles: re-run multi-threshold REGE with other tie thresholds.
#   C. Flow estimates: re-run with other Abel & Cohen (2022) estimation methods
#      (the pipeline default is da_pb_closed).
# Agreement between any two groupings = adjusted Rand index (ARI: 1 = identical,
# 0 = no better than chance) and % of country-periods placed in the same tier.
#
# Writes to output/undesa/robustness/:
#   volume_tiers.csv              tier of every country-period under each rule
#   volume_tier_profile.csv       description of the decade tiers (the main grouping)
#   rege_threshold_runs.csv       k, role sizes, agreement for each REGE variant
#   agreement_ari.csv             ARI between every pair of groupings
#   estimate_methods_usa_deu.csv  USA and Germany flows under each Abel method
#   R1_volume_cutoffs.png, R2_agreement.png, R3_tier_switches.png
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(readr); library(ggplot2)
})
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
cfg_file <- commandArgs(trailingOnly = TRUE)   # optional: Rscript <script> my_config.yml
cfg <- read_config(if (length(cfg_file)) cfg_file[1] else "config.yml")
ds <- "undesa"
outdir <- dirname(out_path(cfg, ds, "robustness", "x"))
rc <- modifyList(cfg$roles, cfg$flows[[ds]]$roles)
base_thr <- unlist(rc$thresholds)

flows <- read_flows(cfg, ds)
cl_list <- read_csv(in_root(cfg, cfg$country_list), show_col_types = FALSE) |>
  select(iso3 = iso3c, country = name)
nm <- node_metrics(flows, 1) |>
  mutate(volume = in_strength + out_strength, node = node_id(iso3, period)) |>
  left_join(cl_list, by = "iso3")

# ---- helper: REGE + clustering for one flow table and one threshold set ------
run_rege <- function(fl, thr, k_fixed = NULL) {
  A <- build_stacked_network(fl, thr)
  sim <- rege(A, rc$rege_iterations, verbose = FALSE)$similarity
  rr <- rc; if (!is.null(k_fixed)) rr$k <- k_fixed
  cl <- cluster_roles(1 - sim, rr)
  cl$role[match(nm$node, names(cl$role))]          # aligned with nm rows
}
ari <- function(a, b) { ok <- !is.na(a) & !is.na(b); adjusted_rand(a[ok], b[ok]) }
same <- function(a, b) mean(a == b, na.rm = TRUE)

# =============================================================================
# A. Volume tiers under different cut-off rules (5-year flows, people)
# =============================================================================
lv <- log10(nm$volume)
jenks <- function(k) {
  br <- classInt::classIntervals(lv, k, style = "fisher")$brks
  cut(lv, br, include.lowest = TRUE, labels = FALSE)
}
q <- function(k) cut(lv, quantile(lv, seq(0, 1, length.out = k + 1)), include.lowest = TRUE, labels = FALSE)
fixed <- function(b) cut(nm$volume, c(0, b, Inf), labels = FALSE)

tiers <- tibble(
  node = nm$node,
  `fixed 10k/100k/1M (4 tiers)`  = fixed(c(1e4, 1e5, 1e6)),
  `fixed 30k/300k/3M (4 tiers)`  = fixed(c(3e4, 3e5, 3e6)),
  `fixed 10k/100k (3 tiers)`     = fixed(c(1e4, 1e5)),
  `fixed 100k/1M (3 tiers)`      = fixed(c(1e5, 1e6)),
  `terciles`                     = q(3),
  `quartiles`                    = q(4),
  `natural breaks, 3`            = jenks(3),
  `natural breaks, 4`            = jenks(4)
)
breaks_used <- list(
  `fixed 10k/100k/1M (4 tiers)` = c(1e4, 1e5, 1e6),
  `fixed 30k/300k/3M (4 tiers)` = c(3e4, 3e5, 3e6),
  `terciles`  = 10^quantile(lv, c(1/3, 2/3)),
  `quartiles` = 10^quantile(lv, c(.25, .5, .75)),
  `natural breaks, 3` = 10^classInt::classIntervals(lv, 3, style = "fisher")$brks[2:3],
  `natural breaks, 4` = 10^classInt::classIntervals(lv, 4, style = "fisher")$brks[2:4])

# Description of the main grouping (decade breaks)
tier_names <- c("1 low (< 10k)", "2 medium (10k-100k)", "3 high (100k-1M)", "4 very high (> 1M)")
main <- nm |> mutate(tier = tier_names[tiers$`fixed 10k/100k/1M (4 tiers)`])
profile <- main |>
  group_by(tier) |>
  summarise(country_periods = n(), countries = n_distinct(iso3),
            median_inflow = median(in_strength), median_outflow = median(out_strength),
            share_net_receiving = mean(net_flow > 0),
            median_in_degree = median(in_degree), median_out_degree = median(out_degree),
            .groups = "drop") |>
  left_join(main |> filter(period == max(period)) |> group_by(tier) |>
              slice_max(volume, n = 6, with_ties = FALSE) |>
              summarise(largest_2020 = paste(coalesce(country, iso3), collapse = "; ")), by = "tier")
write_csv(profile, file.path(outdir, "volume_tier_profile.csv"))
msg("Volume tiers (10k/100k/1M):"); print(as.data.frame(profile[, 1:6]))

# how often does a country-period change tier when the cut-offs move?
stab <- tiers |> select(-node) |> as.matrix()
msg("country-periods in the same tier under the two 4-tier fixed rules: ",
    round(100 * same(stab[, 1], stab[, 2]), 1), "%")

# =============================================================================
# B. REGE roles with other tie thresholds (same flows)
# =============================================================================
thr_sets <- list(
  `REGE 1/1k/10k/100k (base)` = c(1, 1e3, 1e4, 1e5),
  `REGE 1/10k/100k`           = c(1, 1e4, 1e5),
  `REGE 1/1k/10k`             = c(1, 1e3, 1e4),
  `REGE 1/100/1k/10k/100k`    = c(1, 1e2, 1e3, 1e4, 1e5),
  `REGE 1/3k/30k/300k`        = c(1, 3e3, 3e4, 3e5),
  `REGE 1/1k/10k/100k/1M`     = c(1, 1e3, 1e4, 1e5, 1e6))
roles_auto <- list(); roles_k5 <- list()
for (nmx in names(thr_sets)) {
  msg("  ", nmx)
  A <- build_stacked_network(flows, thr_sets[[nmx]])
  sim <- rege(A, rc$rege_iterations, verbose = FALSE)$similarity
  cla <- cluster_roles(1 - sim, rc)
  roles_auto[[nmx]] <- cla$role[match(nm$node, names(cla$role))]
  roles_k5[[nmx]]   <- cutree(cla$hclust, 5)[match(nm$node, names(cla$role))]
}

# =============================================================================
# C. Other Abel & Cohen estimation methods (base thresholds)
# =============================================================================
methods <- c(da_min_closed = "demographic accounting, minimisation, closed",
             da_min_open   = "demographic accounting, minimisation, open",
             sd_drop_neg   = "stock differencing, negatives dropped")
est_roles <- list(); est_vol <- list(); usa_deu <- list()
for (m in c("da_pb_closed", names(methods))) {
  c2 <- cfg; c2$flows[[ds]]$value <- m
  fl <- read_flows(c2, ds)
  e <- fl |> filter(flow >= 1)
  v <- full_join(e |> group_by(iso3 = origin, period) |> summarise(out = sum(flow), .groups = "drop"),
                 e |> group_by(iso3 = destination, period) |> summarise(inn = sum(flow), .groups = "drop"),
                 by = c("iso3", "period")) |>
    mutate(across(c(out, inn), \(x) coalesce(x, 0)), node = node_id(iso3, period))
  vol_m <- (v$out + v$inn)[match(nm$node, v$node)]
  est_vol[[m]] <- cut(vol_m, c(0, 1e4, 1e5, 1e6, Inf), labels = FALSE)
  usa_deu[[m]] <- v |> filter(iso3 %in% c("USA", "DEU")) |>
    transmute(method = m, iso3, period, outflow = out, inflow = inn, net = inn - out)
  if (m != "da_pb_closed") {
    msg("  REGE base thresholds, ", m)
    est_roles[[m]] <- run_rege(fl, base_thr)
  }
}
write_csv(bind_rows(usa_deu), file.path(outdir, "estimate_methods_usa_deu.csv"))

# =============================================================================
# Agreement
# =============================================================================
groupings <- c(
  as.list(tiers |> select(-node)),
  setNames(roles_auto, paste0(names(roles_auto), ", k auto")),
  setNames(roles_k5, paste0(names(roles_k5), ", k = 5")),
  setNames(est_vol[names(methods)], paste0("tiers 10k/100k/1M, ", names(methods))),
  setNames(est_roles, paste0("REGE base, ", names(est_roles))))
G <- names(groupings)
ARI <- outer(G, G, Vectorize(function(a, b) ari(groupings[[a]], groupings[[b]])))
dimnames(ARI) <- list(G, G)
write_csv(as_tibble(ARI, rownames = "grouping"), file.path(outdir, "agreement_ari.csv"))

base_nm <- "REGE 1/1k/10k/100k (base), k auto"
runs <- tibble(grouping = G,
               k = map_int(groupings, \(x) n_distinct(na.omit(x))),
               sizes = map_chr(groupings, \(x) paste(table(x), collapse = "/")),
               ari_vs_base_rege = ARI[, base_nm],
               ari_vs_decade_tiers = ARI[, "fixed 10k/100k/1M (4 tiers)"])
write_csv(runs, file.path(outdir, "rege_threshold_runs.csv"))
write_csv(bind_cols(tiers, as_tibble(roles_auto)), file.path(outdir, "volume_tiers.csv"))
print(as.data.frame(runs |> mutate(across(where(is.numeric), \(x) round(x, 2)))))

# =============================================================================
# Figures
# =============================================================================
theme_set(theme_minimal(11) + theme(plot.title.position = "plot", legend.position = "bottom"))

cuts <- imap_dfr(breaks_used, \(b, n) tibble(rule = n, cut = b))
pR1 <- ggplot(nm, aes(volume)) +
  geom_histogram(bins = 60, fill = "grey70", colour = "white") +
  geom_vline(data = cuts, aes(xintercept = cut), colour = "#c0392b") +
  facet_wrap(~rule, ncol = 2) +
  scale_x_log10(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
  labs(title = "Migrants: total volume (in + out, 5-year flows) of every country-period",
       subtitle = "Red lines: cut-offs under each rule. The distribution has no deep gaps, so tiers are set by the cut-offs.",
       x = "People moving in or out (log scale)", y = "Country-periods")
ggsave(file.path(outdir, "R1_volume_cutoffs.png"), pR1, width = 9, height = 8, dpi = 200, bg = "white")

show <- c(names(tiers)[-1], paste0(names(thr_sets), ", k auto"),
          paste0("tiers 10k/100k/1M, ", names(methods)), paste0("REGE base, ", names(methods)))
Al <- as_tibble(ARI[show, show], rownames = "a") |> pivot_longer(-a, names_to = "b", values_to = "ari") |>
  mutate(a = factor(a, show), b = factor(b, rev(show)))
pR2 <- ggplot(Al, aes(a, b, fill = ari)) + geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", ari)), size = 2.3) +
  scale_fill_gradient(low = "white", high = "#1f3a5f", limits = c(0, 1)) +
  labs(title = "Agreement between migrant groupings (adjusted Rand index)",
       subtitle = "Volume tiers, REGE roles with other tie thresholds, and other Abel & Cohen flow estimates",
       x = NULL, y = NULL, fill = "ARI") +
  theme(axis.text.x = element_text(angle = 50, hjust = 1, size = 7), axis.text.y = element_text(size = 7),
        panel.grid = element_blank())
ggsave(file.path(outdir, "R2_agreement.png"), pR2, width = 10, height = 9, dpi = 200, bg = "white")

# which countries change tier when the decade cut-offs shift by half a decade?
sw <- main |> mutate(alt = tiers$`fixed 30k/300k/3M (4 tiers)`, base = tiers$`fixed 10k/100k/1M (4 tiers)`) |>
  count(base, alt) |> group_by(base) |> mutate(share = n / sum(n))
pR3 <- ggplot(sw, aes(factor(alt), factor(base, 4:1), fill = share)) + geom_tile(colour = "white") +
  geom_text(aes(label = n), size = 3.5) +
  scale_fill_gradient(low = "#f2f5f9", high = "#1f3a5f", labels = scales::percent) +
  labs(title = "Tier membership when cut-offs move from 10k/100k/1M to 30k/300k/3M",
       x = "Tier with cut-offs 30k / 300k / 3M", y = "Tier with cut-offs 10k / 100k / 1M",
       fill = "Share of row") + theme(panel.grid = element_blank())
ggsave(file.path(outdir, "R3_tier_switches.png"), pR3, width = 6.5, height = 5, dpi = 200, bg = "white")
msg("Done: ", outdir)
