# =============================================================================
# validate_against_original.R
#
# Compares role assignments from this pipeline with the ones saved by the
# original FINAL_*_role_assessment.R scripts (results/*_full_roles_k3.csv).
#
#   Rscript validate_against_original.R "legacy/results"
#
# Findings when this was written (Sept 2026):
#   * UNHCR: identical partition (adjusted Rand index = 1.00).
#   * UN DESA: plain REGE on the >=1-migrant network puts every country in one
#     role (every country both sends and receives, so all are regularly
#     equivalent). The original 3 roles are almost entirely a function of
#     in/out degree (98% predicted by degree alone) and cannot be reproduced
#     with standard REGE.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
orig_dir <- if (length(args)) args[1] else "legacy/results"

for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
cfg <- read_config("config.yml")
clist <- read_csv(in_root(cfg, cfg$country_list), show_col_types = FALSE)
xw <- load_crosswalk(cfg)

read_original <- function(file) {
  x <- read_csv(file, show_col_types = FALSE)
  tibble(name = sub("\\.(19|20)\\d\\d$", "", x$id),
         period = as.integer(sub(".*\\.", "", x$id)),
         original_role = x$role) |>
    mutate(iso3 = coalesce(clist$iso3c[match(name, clist$name)],
                           xw$iso3[match(name, xw$alias)]))
}

for (ds in c("unhcr", "undesa")) {
  msg("==== ", ds)
  orig <- read_original(file.path(orig_dir, paste0(toupper(ds), "_full_roles_k3.csv")))
  flows <- read_flows(cfg, ds)

  # the original settings: one tie threshold of 1 person, k = 3
  sim <- rege(build_stacked_network(flows, 1), 3, verbose = FALSE)$similarity
  for (mode in c("rege_distance", "legacy_profiles")) {
    cl <- cluster_roles(1 - sim, list(cluster_on = mode, linkage = "average",
                                      k = 3, k_range = c(2, 6)))
    new <- tibble(node = names(cl$role), role = cl$role) |>
      separate(node, c("iso3", "period"), sep = "_", convert = TRUE) |>
      inner_join(orig, by = c("iso3", "period"))
    msg(sprintf("  %-16s adjusted Rand vs original = %.3f  (n = %d)", mode,
                adjusted_rand(new$role, new$original_role), nrow(new)))
    print(table(new = new$role, original = new$original_role))
  }

  # how much of the original role is explained by degree alone?
  nm <- node_metrics(flows, 1) |> inner_join(orig, by = c("iso3", "period"))
  fit <- nnet::multinom(factor(original_role) ~ log1p(in_degree) + log1p(out_degree),
                        data = nm, trace = FALSE)
  msg(sprintf("  original role predicted from in/out degree alone: %.1f%% correct",
              100 * mean(predict(fit) == nm$original_role)))
}
