# =============================================================================
# migrant_backbone.R — backbone migrant roles: sensitivity and drivers
#
#   Rscript analysis/migrant_backbone.R     (~3 min; after run_pipeline.R)
#
# The backbone keeps, for every country and period, only the net migrant flows
# that carry at least `share` of its total net outflow ("its main net
# destinations"). REGE + clustering on that network gives positions such as
# sources (not a main destination for anyone), intermediates (a main
# destination for some, sending onward) and terminal destinations.
#
#   1. Sensitivity: share = 5%, 10%, 20% x Abel & Cohen pseudo-Bayesian and
#      minimisation estimates; agreement with the reference (10%, pseudo-Bayesian)
#   2. Drivers: mixed-effects logit per backbone role (reference version), same
#      drivers and specification as the main pipeline
#
# Writes to output/undesa_backbone/:
#   backbone_sensitivity.csv, backbone_roles.csv, binary_role_models.csv,
#   B1_backbone_sensitivity.png, B2_backbone_odds_ratios.png,
#   B3_backbone_timeline_africa.png
#
# Backbone definition, written out. Let N_ij,t >= 0 be the net flow from
# country i to country j in period t (only the net sender of each pair has a
# link; see read_flows(net = TRUE)). Country i's total net outflow is
# O_i,t = sum_j N_ij,t. The backbone keeps link i -> j if
#
#      N_ij,t >= s * O_i,t        (s = 0.10 in the reference version)
#
# i.e. links carrying at least 10% of the origin's net outflow. A country
# keeps at most 1/s = 10 outgoing links; a country whose net outflow is spread
# thinly (no destination with >= 10%) keeps none and appears only through its
# backbone inflows, if any.
# Roles are then found on the backbone exactly as in the main pipeline
# (binary REGE, average-linkage clustering), with k = 3 roles reported.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(readr); library(ggplot2)
})
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
cfg_file <- commandArgs(trailingOnly = TRUE)   # optional: Rscript <script> my_config.yml
cfg <- read_config(if (length(cfg_file)) cfg_file[1] else "config.yml")
ds <- "undesa_backbone"                       # name of the output subfolder
outdir <- dirname(out_path(cfg, ds, "x"))
rc <- cfg$roles; rc$k <- "auto"              # roles settings; k chosen by silhouette for the diagnostics
shares <- c(0.05, 0.10, 0.20)                 # backbone thresholds s compared in the sensitivity analysis
# two Abel & Cohen estimation methods (column names in the flow file -> labels)
methods <- c(da_pb_closed = "pseudo-Bayesian", da_min_closed = "minimisation")

# Net migrant flows from one Abel & Cohen estimate: a copy of the config with
# that flow column and net = TRUE, read with the standard read_flows().
read_net <- function(method) {
  c2 <- cfg; c2$flows$undesa$value <- method; c2$flows$undesa$net <- TRUE
  read_flows(c2, "undesa")
}
# Keep the links that carry at least `share` of their origin's total net
# outflow in that period (N_ij,t >= s * O_i,t; see header).
backbone <- function(fl, share) fl |> group_by(origin, period) |>
  filter(flow >= share * sum(flow)) |> ungroup()

# Label each role from the in/out pattern of its members on the backbone
#   s   share of members with no backbone inflow  (nobody's main destination)
#   r   share of members with no backbone outflow (no main net destination)
#   hub median number of countries for which members are a main destination
# > 50% with no inflow -> source; > 50% with no outflow -> terminal;
# otherwise intermediate.
label_backbone <- function(roles_df) {
  p <- roles_df |> group_by(role) |>
    summarise(s = mean(in_degree == 0), r = mean(out_degree == 0),
              hub = median(in_degree), .groups = "drop")
  lab <- case_when(p$s > 0.5 ~ "source: not a main destination for anyone",
                   p$r > 0.5 ~ "terminal destination: receives, no main outflow",
                   TRUE ~ "intermediate: main destination for some, sends onward")
  # several intermediate roles -> distinguish by how many countries they serve
  dup <- lab %in% lab[duplicated(lab)]
  lab[dup] <- paste0(lab[dup], " (median ", round(p$hub[dup]), " countries)")
  setNames(lab, p$role)
}

# ---- 1. sensitivity -------------------------------------------------------------
# net flow networks for both estimation methods
nets <- map(names(methods), read_net) |> setNames(names(methods))
runs <- list()
for (m in names(methods)) for (s in shares) {
  key <- sprintf("%s, %d%%", methods[[m]], round(100 * s))
  msg("  backbone ", key)
  bb <- backbone(nets[[m]], s)
  # binary REGE on the backbone (one threshold: any backbone link), then clustering
  sim <- rege(build_stacked_network(bb, 1), rc$rege_iterations, verbose = FALSE)$similarity
  cl <- cluster_roles(1 - sim, rc)
  # keep the automatic-k solution and the k = 3 solution used in the paper
  runs[[key]] <- list(m = m, s = s, flows = bb, hc = cl$hclust, auto = cl$role, k = cl$k,
                      k3 = cutree(cl$hclust, 3), diag = cl$diagnostics)
}
ref_key <- "pseudo-Bayesian, 10%"             # reference version reported in the paper
nodes <- names(runs[[ref_key]]$k3)
# ARI between two role vectors on the reference version's nodes; nodes absent
# from one version (no backbone link at that threshold) are skipped
ari <- function(a, b) { a <- a[nodes]; b <- b[nodes]; ok <- !is.na(a) & !is.na(b); adjusted_rand(a[ok], b[ok]) }
# one row per version: number of backbone links, automatic k and silhouette,
# silhouette and role sizes at k = 3, and agreement with the reference
sens <- imap_dfr(runs, function(r, key) tibble(
  version = key, method = methods[[r$m]], share = r$s,
  ties = nrow(r$flows), k_auto = r$k,
  silhouette_auto = r$diag$avg_silhouette[r$diag$k == r$k],
  silhouette_k3 = r$diag$avg_silhouette[r$diag$k == 3],
  sizes_k3 = paste(sort(table(r$k3), decreasing = TRUE), collapse = "/"),
  ari_vs_reference_k3 = ari(r$k3, runs[[ref_key]]$k3)))
write_csv(sens, file.path(outdir, "backbone_sensitivity.csv"))
print(as.data.frame(sens |> mutate(across(where(is.double), \(x) round(x, 3)))))

# ---- reference roles with labels and metrics ---------------------------------
ref <- runs[[ref_key]]
cl_list <- read_csv(in_root(cfg, cfg$country_list), show_col_types = FALSE) |>
  select(iso3 = iso3c, country = name, region = un_region, area = un_area)
# in/out degree on the BACKBONE (used for labelling), not on the full network
bb_metrics <- node_metrics(ref$flows, 1) |>
  select(iso3, period, bb_in_degree = in_degree, bb_out_degree = out_degree)
roles <- tibble(node = names(ref$k3), role = unname(ref$k3)) |>
  separate(node, c("iso3", "period"), sep = "_", convert = TRUE) |>
  left_join(bb_metrics, by = c("iso3", "period")) |>
  rename(in_degree = bb_in_degree, out_degree = bb_out_degree)
labs_bb <- label_backbone(roles)
# order roles: sources, intermediates, terminal
ord <- order(match(sub(":.*", "", labs_bb), c("source", "intermediate", "terminal destination")))
recode <- setNames(seq_along(ord), names(labs_bb)[ord])
# renumber roles 1 = source (M1), 2 = intermediate (M2), 3 = terminal (M3);
# add net flow and volume from the FULL net network (descriptive table)
roles <- roles |>
  mutate(label = labs_bb[as.character(role)], role = unname(recode[as.character(role)])) |>
  left_join(node_metrics(nets$da_pb_closed, 1) |> select(iso3, period, net_flow, in_strength, out_strength),
            by = c("iso3", "period")) |>
  left_join(cl_list, by = "iso3")
stab <- role_stability(roles)
write_csv(roles, file.path(outdir, "backbone_roles.csv"))
msg("  reference roles: ", paste(names(table(roles$label)), table(roles$label), sep = " = ", collapse = "; "))
msg("  persistence between periods: ", round(100 * stab$overall_persistence, 1), "%")
print(roles |> group_by(role, label) |>
        summarise(n = n(), share_net_receiving = mean(net_flow > 0),
                  median_volume = median(in_strength + out_strength), .groups = "drop"))
write_csv(stab$transition_matrix, file.path(outdir, "role_transitions.csv"))

# ---- 2. drivers -------------------------------------------------------------------
# Same drivers (5-year aggregates written by run_pipeline.R) and the same
# model specification as the main migrant models: prepare_model_data() and
# fit_binary_role_models() from R/models.R, incl. income-tercile drought slopes
# and Wald tests when the config asks for them.
drivers <- read_csv(out_path(cfg, "undesa", "drivers_by_country_period.csv"), show_col_types = FALSE)
md <- prepare_model_data(cfg, roles, drivers)
msg("  model data: ", nrow(md$data), " country-periods")
bin <- fit_binary_role_models(cfg, md) |>
  left_join(distinct(roles, role, label), by = "role")
write_csv(bin, file.path(outdir, "binary_role_models.csv"))
print(bin |> filter(!grepl("Intercept|period", term), p.value < 0.05) |>
        select(role, term, odds_ratio, or_low, or_high, p.value, n_events, unstable))

# ---- figures ------------------------------------------------------------------------
theme_set(theme_minimal(11) + theme(plot.title.position = "plot", legend.position = "bottom"))
# B1: agreement of each version with the reference
S <- sens |> mutate(version = factor(version, version))
pB1 <- ggplot(S, aes(ari_vs_reference_k3, version, fill = method)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%.2f  (sizes %s)", ari_vs_reference_k3, sizes_k3)), hjust = -0.05, size = 3) +
  scale_x_continuous(limits = c(0, 1.6), breaks = seq(0, 1, 0.25)) +
  scale_fill_manual(values = c(`pseudo-Bayesian` = "#1f3a5f", minimisation = "#8aa4c1")) +
  labs(title = "Backbone roles: sensitivity to the backbone share and the Abel & Cohen method",
       subtitle = "Agreement (adjusted Rand index, 3 roles) with the reference version (pseudo-Bayesian, 10%).",
       x = "Adjusted Rand index", y = NULL, fill = NULL)
ggsave(file.path(outdir, "B1_backbone_sensitivity.png"), pB1, width = 9, height = 4, dpi = 200, bg = "white")

# B2: odds ratios per backbone role (unstable estimates and Wald rows have no
# odds ratio and are not drawn: unstable is TRUE or odds_ratio is NA)
pB2 <- bin |> filter(!grepl("Intercept|period", term), !unstable) |>
  mutate(sig = p.value < 0.05, facet = paste0(role, ": ", sub(":.*", "", label), " (n = ", n_events, ")")) |>
  ggplot(aes(odds_ratio, term, colour = sig)) +
  geom_vline(xintercept = 1, colour = "grey60") +
  geom_errorbarh(aes(xmin = or_low, xmax = or_high), height = 0) + geom_point() +
  scale_x_log10() + facet_wrap(~facet) +
  scale_colour_manual(values = c(`FALSE` = "grey55", `TRUE` = "#b34b4b"), labels = c("p >= 0.05", "p < 0.05")) +
  labs(title = "Migrants, backbone roles: drivers of role membership",
       x = "Odds ratio per 1 SD (95% CI, log scale)", y = NULL, colour = NULL)
ggsave(file.path(outdir, "B2_backbone_odds_ratios.png"), pB2, width = 10, height = 7, dpi = 200, bg = "white")

# B3: role of each African country over time
rl <- roles |> mutate(role_lab = factor(paste0(role, ": ", sub(":.*", "", label))))
pB3 <- rl |> filter(area == "Africa") |>
  ggplot(aes(period, forcats::fct_rev(coalesce(country, iso3)), fill = role_lab)) +
  geom_tile(colour = "white") +
  scale_fill_manual(values = c("#7b5ea7", "#e08a2c", "#333333", "#3a9a5b")[seq_len(n_distinct(rl$role_lab))]) +
  labs(title = "Migrants, backbone roles: African countries", x = NULL, y = NULL, fill = "Role") +
  theme(panel.grid = element_blank(), axis.text.y = element_text(size = 7))
ggsave(file.path(outdir, "B3_backbone_timeline_africa.png"), pB3, width = 7, height = 10, dpi = 200, bg = "white")
msg("Done: ", outdir)
