# =============================================================================
# migrant_net_roles.R — regular-equivalence roles on NET migrant flows
#
#   Rscript analysis/migrant_net_roles.R      (~5 min; ~25 min with run_valued)
#
# Gross Abel & Cohen flows are dense (almost every country sends to and receives
# from almost every type of country) and the estimation methods differ mainly in
# the size of counter-flows. Net flows (A->B minus B->A) are almost identical
# across the closed demographic-accounting methods (r = 0.998), so roles built
# on them should not depend on the method. Three versions of the net network:
#   net      one directed tie per pair, from net sender to net receiver (>= 1 person)
#   backbone only each country's main net destinations (ties carrying >= 10% of
#            its total net outflow): "who are you a principal destination for?"
#   valued   White's REGGE on log(1 + net flow) (slow; set run_valued <- TRUE)
# Each is run for the pseudo-Bayesian (da_pb_closed, pipeline default) and the
# minimisation (da_min_closed) estimates; agreement = adjusted Rand index.
#
# Writes to output/undesa_net/:
#   net_role_variants.csv     k, role sizes, silhouette, agreement across methods
#   net_roles_by_country_period.csv
#   N1_net_role_scatter.png, N2_method_agreement.png
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(readr); library(ggplot2)
})
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
cfg_file <- commandArgs(trailingOnly = TRUE)   # optional: Rscript <script> my_config.yml
cfg <- read_config(if (length(cfg_file)) cfg_file[1] else "config.yml")
run_valued <- FALSE
backbone_share <- 0.10
outdir <- dirname(out_path(cfg, "undesa_net", "x"))
rc <- cfg$roles; rc$k <- "auto"

read_net <- function(method) {
  c2 <- cfg; c2$flows$undesa$value <- method; c2$flows$undesa$net <- TRUE
  read_flows(c2, "undesa")
}
backbone <- function(fl) fl |> group_by(origin, period) |>
  filter(flow >= backbone_share * sum(flow)) |> ungroup()

fit <- function(fl, valued = FALSE) {
  sim <- if (valued) rege_valued(build_stacked_valued(fl, 1), rc$rege_iterations, "log1p")$similarity
         else rege(build_stacked_network(fl, 1), rc$rege_iterations, verbose = FALSE)$similarity
  cl <- cluster_roles(1 - sim, rc)
  list(role = cl$role, diag = cl$diagnostics, k = cl$k)
}

variants <- c("net", "backbone", if (run_valued) "valued")
methods <- c("da_pb_closed", "da_min_closed")
res <- list(); nets <- list()
for (m in methods) {
  fl <- read_net(m); nets[[m]] <- fl
  for (v in variants) {
    msg("  ", v, " / ", m)
    g <- switch(v, net = fl, backbone = backbone(fl), valued = fl)
    res[[paste(v, m)]] <- c(fit(g, v == "valued"), list(flows = g))
  }
}

# ---- table of variants ------------------------------------------------------
gross <- read_csv(out_path(cfg, "undesa", "roles_by_country_period.csv"), show_col_types = FALSE) |>
  mutate(node = node_id(iso3, period))
align <- function(r) unname(r$role[gross$node])
ari <- function(a, b) { ok <- !is.na(a) & !is.na(b); adjusted_rand(a[ok], b[ok]) }
tab <- map_dfr(variants, function(v) {
  a <- res[[paste(v, "da_pb_closed")]]; b <- res[[paste(v, "da_min_closed")]]
  tibble(network = v, k = a$k, sizes = paste(table(a$role), collapse = "/"),
         silhouette = a$diag$avg_silhouette[a$diag$k == a$k],
         k_min_closed = b$k, ari_pb_vs_min = ari(align(a), align(b)),
         ari_vs_gross_roles = ari(align(a), gross$role))
})
write_csv(tab, file.path(outdir, "net_role_variants.csv"))
print(as.data.frame(tab |> mutate(across(where(is.double), \(x) round(x, 2)))))

# ---- roles with profiles (pseudo-Bayesian) ----------------------------------
cl_list <- read_csv(in_root(cfg, cfg$country_list), show_col_types = FALSE) |>
  select(iso3 = iso3c, country = name)
nm <- node_metrics(nets$da_pb_closed, 1) |> mutate(node = node_id(iso3, period))
lab_net <- function(r, nm) {
  p <- nm |> mutate(role = r[node]) |> group_by(role) |>
    summarise(s = mean(pure_sender), rcv = mean(pure_receiver), .groups = "drop")
  setNames(case_when(p$s > 0.5 ~ "net sender (few or no net inflows)",
                     p$rcv > 0.5 ~ "net receiver (no or few net outflows)",
                     TRUE ~ "both (net sender to some, net receiver from others)"), p$role)
}
bb_nm <- node_metrics(res$`backbone da_pb_closed`$flows, 1) |> mutate(node = node_id(iso3, period))
lab_bb <- function(r) {
  p <- bb_nm |> mutate(role = r[node]) |> group_by(role) |>
    summarise(s = mean(pure_sender), rcv = mean(pure_receiver), .groups = "drop")
  setNames(case_when(p$s > 0.5 ~ "source: not a main destination for anyone",
                     p$rcv > 0.5 ~ "terminal destination: receives, no main outflow",
                     TRUE ~ "intermediate: main destination for some, sends onward"), p$role)
}
rn <- res$`net da_pb_closed`$role; rb <- res$`backbone da_pb_closed`$role
roles <- nm |>
  mutate(net_role = rn[node], net_label = lab_net(rn, nm)[as.character(net_role)],
         backbone_role = rb[node], backbone_label = lab_bb(rb)[as.character(backbone_role)]) |>
  left_join(cl_list, by = "iso3")
write_csv(roles |> select(-node), file.path(outdir, "net_roles_by_country_period.csv"))
print(roles |> count(net_label, backbone_label))

# ---- figures --------------------------------------------------------------------
theme_set(theme_minimal(11) + theme(plot.title.position = "plot", legend.position = "bottom"))
pd <- bind_rows(
  roles |> transmute(panel = "A. Net network (all net ties >= 1 person)", label = net_label,
                     x = out_strength, y = in_strength, iso3, period, country, vol = in_strength + out_strength),
  roles |> transmute(panel = "B. Backbone (each country's main net destinations)",
                     label = coalesce(backbone_label, "not in backbone"),
                     x = out_strength, y = in_strength, iso3, period, country, vol = in_strength + out_strength))
labs_pts <- pd |> filter(period == max(period)) |> group_by(panel, label) |>
  slice_max(vol, n = 3, with_ties = FALSE) |> ungroup()
pN1 <- ggplot(pd, aes(log10(x + 1), log10(y + 1), colour = label)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey75", linetype = 2) +
  geom_point(alpha = 0.5, size = 1.3) +
  ggrepel::geom_text_repel(data = labs_pts, aes(label = paste(coalesce(country, iso3), period)),
                           size = 2.4, show.legend = FALSE, max.overlaps = 30, min.segment.length = 0) +
  facet_wrap(~panel, ncol = 1) +
  scale_colour_manual(values = c(
    "both (net sender to some, net receiver from others)" = "#2c7bb6",
    "net receiver (no or few net outflows)" = "#3a9a5b",
    "net sender (few or no net inflows)" = "#b04a4a",
    "intermediate: main destination for some, sends onward" = "#e08a2c",
    "source: not a main destination for anyone" = "#7b5ea7",
    "terminal destination: receives, no main outflow" = "#333333",
    "not in backbone" = "#bbbbbb")) +
  labs(title = "Migrants, net-flow network: regular-equivalence roles",
       subtitle = "Each point is a country-period (pseudo-Bayesian estimates).\nAxes: total net outflow and net inflow, log10(people + 1).",
       x = "Net outflow (sum over countries it is a net sender to)",
       y = "Net inflow (sum over countries it is a net receiver from)", colour = NULL) +
  guides(colour = guide_legend(ncol = 2, override.aes = list(size = 3, alpha = 1)))
ggsave(file.path(outdir, "N1_net_role_scatter.png"), pN1, width = 8, height = 10, dpi = 200, bg = "white")

# agreement between the two Abel & Cohen methods, gross vs net
rob <- read_csv(out_path(cfg, "undesa", "robustness", "agreement_ari.csv"), show_col_types = FALSE)
gross_ari <- rob$`REGE base, da_min_closed`[rob$grouping == "REGE 1/1k/10k/100k (base), k auto"]
tier_ari  <- rob$`tiers 10k/100k/1M, da_min_closed`[rob$grouping == "fixed 10k/100k/1M (4 tiers)"]
ag <- tibble(grouping = c("Gross flows: REGE roles (1/1k/10k/100k)", "Gross flows: volume tiers (10k/100k/1M)",
                          paste0("Net flows: REGE roles (", tab$network, ")")),
             ari = c(gross_ari, tier_ari, tab$ari_pb_vs_min)) |>
  mutate(grouping = factor(grouping, rev(grouping)))
pN2 <- ggplot(ag, aes(ari, grouping)) + geom_col(fill = "#1f3a5f", width = 0.6) +
  geom_text(aes(label = sprintf("%.2f", ari)), hjust = -0.15, size = 3.5) +
  scale_x_continuous(limits = c(0, 1.05)) +
  labs(title = "Do the groups depend on the Abel & Cohen estimation method?",
       subtitle = "Agreement (adjusted Rand index) between groupings built from pseudo-Bayesian and minimisation estimates.\n1 = identical, 0 = no better than chance.",
       x = "Adjusted Rand index", y = NULL)
ggsave(file.path(outdir, "N2_method_agreement.png"), pN2, width = 8, height = 3.6, dpi = 200, bg = "white")
msg("Done: ", outdir)
