# =============================================================================
# roles.R — from flows to roles
#
#   1. build_stacked_network(): one node per country-period ("KEN_2005").
#      Ties only link countries within the same period, but because all periods
#      sit in one matrix, REGE compares e.g. Kenya-1995 with Kenya-2015 and with
#      Uganda-2005 on the same scale. That is what lets us track role change.
#   2. rege() (rege.R) -> similarity; distance = 1 - similarity
#   3. cluster_roles(): hierarchical clustering + k diagnostics
#   4. node_metrics(), role_profile(), role_stability(): describe the roles
#
# run_roles() at the bottom runs steps 1-4 for one dataset and writes the
# outputs listed there.
# =============================================================================

# Node name for a country in a period, e.g. node_id("KEN", 2005) -> "KEN_2005".
# Split back into iso3 and period in run_roles() with separate(sep = "_").
node_id <- function(iso3, period) paste(iso3, period, sep = "_")

# Build one 0/1 adjacency matrix per tie threshold over the same node set.
#   flows       output of read_flows() (origin, destination, period, flow)
#   thresholds  vector of minimum flow sizes, e.g. c(1, 1000, 100000)
# The node set is every country-period with at least one flow at the LOWEST
# threshold, so all matrices have identical rows/columns (needed by rege()).
# A[i, j] = 1 if country-period i sent at least t people to country-period j.
build_stacked_network <- function(flows, thresholds = 1) {
  base <- flows |> filter(flow >= min(thresholds))
  nodes <- sort(unique(c(node_id(base$origin, base$period),
                         node_id(base$destination, base$period))))
  n <- length(nodes)
  A_list <- lapply(thresholds, function(t) {
    e <- flows |> filter(flow >= t)
    A <- matrix(0, n, n, dimnames = list(nodes, nodes))
    # set all ties at once using a two-column (row, column) index matrix
    A[cbind(match(node_id(e$origin, e$period), nodes),
            match(node_id(e$destination, e$period), nodes))] <- 1
    A
  })
  names(A_list) <- paste0(">=", thresholds)
  msg("  stacked network: ", n, " country-periods, ",
      sum(A_list[[1]]), " ties at threshold ", min(thresholds))
  A_list
}

# Within-cluster sum of squares computed from a distance matrix.
# For a cluster C with |C| members: WCSS_C = sum_{i,j in C} D[i,j]^2 / (2|C|),
# which equals the sum of squared distances to the cluster centroid when D is
# Euclidean. Reported only as an elbow diagnostic; k is chosen by silhouette.
wcss_from_dist <- function(D, cl) {
  sum(vapply(unique(cl), function(k) {
    idx <- which(cl == k)
    sum(D[idx, idx]^2) / (2 * length(idx))
  }, numeric(1)))
}

# Hierarchical clustering of the REGE distances and choice of the number of roles.
#   D          n x n distance matrix (1 - REGE similarity)
#   cfg_roles  the roles: block of the config (linkage, k, k_range, cluster_on)
# Returns the tree, the chosen k, each node's role, and a table of diagnostics
# for every k in k_range:
#   avg_silhouette   mean silhouette width s(i) = (b(i) - a(i)) / max(a(i), b(i)),
#                    a(i) = mean distance of i to its own cluster,
#                    b(i) = mean distance to the nearest other cluster.
#                    Ranges -1..1; higher = better separated roles.
#   wcss             within-cluster sum of squares (elbow plot)
#   smallest_cluster size of the smallest role (flags tiny splinter roles)
cluster_roles <- function(D, cfg_roles) {
  if (cfg_roles$cluster_on == "legacy_profiles") {
    # Reproduces the original scripts: factoextra::hcut(x = matrix) treats the
    # distance matrix as data, takes Euclidean distances between its rows, and
    # (because 'method =' is not an hcut argument) uses the default Ward.D2.
    d <- dist(D)
    hc <- hclust(d, method = "ward.D2")
  } else {
    # Main analysis: cluster the REGE distances directly with the configured
    # linkage (average linkage: distance between two clusters = mean distance
    # between their members).
    d <- as.dist(D)
    hc <- hclust(d, method = cfg_roles$linkage)
  }
  ks <- seq(cfg_roles$k_range[1], cfg_roles$k_range[2])
  Dm <- as.matrix(d)
  diag_tbl <- map_dfr(ks, function(k) {
    cl <- cutree(hc, k)          # cut the tree into k groups
    tibble(k = k,
           avg_silhouette = mean(cluster::silhouette(cl, d)[, "sil_width"]),
           wcss = wcss_from_dist(Dm, cl),
           smallest_cluster = min(table(cl)))
  })
  # k: "auto" = the k with the highest average silhouette; otherwise the fixed
  # number given in the config
  k <- if (identical(cfg_roles$k, "auto")) {
    diag_tbl$k[which.max(diag_tbl$avg_silhouette)]
  } else as.integer(cfg_roles$k)
  list(hclust = hc, k = k, role = cutree(hc, k), diagnostics = diag_tbl)
}

# Per country-period network statistics on the weighted graph (threshold 1).
# Computed separately for each period's network:
#   in_degree / out_degree      number of origins / destinations
#   in_strength / out_strength  total people received / sent
#   betweenness                 how often the country lies on shortest paths
#   net_flow                    in_strength - out_strength (> 0 = net receiver)
#   pure_sender / pure_receiver no inflows / no outflows at all
node_metrics <- function(flows, threshold = 1) {
  e <- flows |> filter(flow >= threshold)
  per <- split(e, e$period)
  imap_dfr(per, function(df, p) {
    # directed graph with the flow as edge weight
    g <- igraph::graph_from_data_frame(df |> select(origin, destination, weight = flow))
    tibble(
      iso3 = igraph::V(g)$name,
      period = as.integer(p),
      in_degree = igraph::degree(g, mode = "in"),
      out_degree = igraph::degree(g, mode = "out"),
      in_strength = igraph::strength(g, mode = "in"),
      out_strength = igraph::strength(g, mode = "out"),
      # betweenness on inverse flow: large flows = short paths
      betweenness = igraph::betweenness(g, weights = 1 / igraph::E(g)$weight)
    )
  }) |>
    mutate(net_flow = in_strength - out_strength,
           pure_sender = in_degree == 0, pure_receiver = out_degree == 0)
}

# Table 1 of the manuscript: mean statistics by role
role_profile <- function(roles_df) {
  roles_df |>
    group_by(role) |>
    summarise(n = n(),
              # medians first: across() below replaces the columns by their means
              median_net_flow = median(net_flow),
              share_net_receiving = mean(net_flow > 0),
              median_volume = median(in_strength + out_strength),
              across(c(in_degree, out_degree, in_strength, out_strength, betweenness),
                     \(x) mean(x, na.rm = TRUE)),
              share_pure_sender = mean(pure_sender),
              share_pure_receiver = mean(pure_receiver),
              .groups = "drop")
}

# Automatic, human-readable role labels from the profile (edit freely).
# Volume level = the role's median total flow (in + out) on fixed decade
# cut-offs; balance = share of the role's country-periods with inflow > outflow
# (a role is called net receiving/sending only if >= 2/3 of its members are).
# Roles where > 90% of members have no inflows (or no outflows) are labelled
# "sender only" ("receiver only"). The labels are descriptive only; the
# manuscript names (R1-R3, M1-M3) are assigned in the figure scripts.
label_roles <- function(profile) {
  lvl <- as.character(cut(profile$median_volume, c(-Inf, 1e4, 1e5, 1e6, Inf),
                          labels = c("low", "medium", "high", "very high")))
  balance <- case_when(profile$share_net_receiving >= 2/3 ~ "net receiving",
                       profile$share_net_receiving <= 1/3 ~ "net sending",
                       TRUE ~ sprintf("mixed (%.0f%% net receiving)", 100 * profile$share_net_receiving))
  profile |>
    mutate(label = case_when(
      share_pure_sender > 0.9 ~ "sender only",
      share_pure_receiver > 0.9 ~ "receiver only",
      TRUE ~ paste0("two-way, ", lvl, " volume, ", balance)
    )) |>
    select(role, label)
}

# How often does each country keep its role from one period to the next?
# Transitions are counted between a country's consecutive OBSERVED periods
# (lag() within country after sorting by period).
#   by_country           periods observed, number of distinct roles, share of
#                        periods in the most common role, number of switches
#   transition_matrix    counts and row probabilities P(role in t | role in t-1)
#   overall_persistence  share of all transitions where the role is unchanged
role_stability <- function(roles_df) {
  trans <- roles_df |>
    arrange(iso3, period) |>
    group_by(iso3) |>
    mutate(prev_role = lag(role), prev_period = lag(period)) |>
    ungroup() |>
    filter(!is.na(prev_role))          # first observed period has no predecessor
  list(
    by_country = roles_df |>
      group_by(iso3) |>
      summarise(periods = n(),
                n_roles_used = n_distinct(role),
                modal_role_share = max(table(role)) / n(),
                .groups = "drop") |>
      left_join(trans |> group_by(iso3) |> summarise(switches = sum(role != prev_role)),
                by = "iso3"),
    transition_matrix = trans |> count(from = prev_role, to = role) |>
      group_by(from) |> mutate(prob = n / sum(n)) |> ungroup(),
    overall_persistence = mean(trans$role == trans$prev_role)
  )
}

# Full role step for one dataset ("undesa" or "unhcr").
# Writes to <output_dir>/<dataset>/:
#   roles_by_country_period.csv   one row per country-period: role, label, metrics
#   role_profile.csv              Table 1 (statistics by role)
#   k_diagnostics.csv             silhouette / WCSS for each k
#   role_stability_by_country.csv, role_transitions.csv
#   rege_similarity.rds           REGE matrix and clustering tree (for figures)
run_roles <- function(cfg, dataset, flows) {
  rc <- cfg$roles
  # a dataset may override any roles setting (e.g. its own tie thresholds)
  rc <- modifyList(rc, cfg$flows[[dataset]]$roles %||% list())
  if (identical(rc$method, "valued")) {
    # robustness option: REGE on flow values (rege_valued in rege.R)
    W <- build_stacked_valued(flows, min(unlist(rc$thresholds)))
    msg("  valued REGE (", rc$rege_iterations, " iterations) on ", nrow(W), " nodes ...")
    sim <- rege_valued(W, rc$rege_iterations, rc$value_transform %||% "log1p",
                       rc$no_self_match %||% FALSE)$similarity
  } else {
    # main analysis: multi-threshold binary REGE
    A_list <- build_stacked_network(flows, unlist(rc$thresholds))
    msg("  REGE (", rc$rege_iterations, " iterations) on ", nrow(A_list[[1]]), " nodes ...")
    sim <- rege(A_list, rc$rege_iterations)$similarity
  }
  D <- 1 - sim                        # similarity in [0,1] -> distance in [0,1]
  cl <- cluster_roles(D, rc)
  msg("  k = ", cl$k, "; role sizes: ", paste(table(cl$role), collapse = " / "))

  # node names "KEN_2005" -> columns iso3 = "KEN", period = 2005 (integer)
  roles <- tibble(node = names(cl$role), role = unname(cl$role)) |>
    separate(node, c("iso3", "period"), sep = "_", convert = TRUE)
  # attach network statistics and country names / UN regions
  cl_list <- read_csv(in_root(cfg, cfg$country_list), show_col_types = FALSE)
  roles <- roles |>
    left_join(node_metrics(flows, min(unlist(rc$thresholds))), by = c("iso3", "period")) |>
    left_join(cl_list |> select(iso3 = iso3c, country = name, region = un_region, area = un_area),
              by = "iso3")

  profile <- role_profile(roles)
  labels <- label_roles(profile)
  roles <- roles |> left_join(labels, by = "role")
  stab <- role_stability(roles)

  write_csv(roles, out_path(cfg, dataset, "roles_by_country_period.csv"))
  write_csv(profile |> left_join(labels, by = "role"), out_path(cfg, dataset, "role_profile.csv"))
  write_csv(cl$diagnostics, out_path(cfg, dataset, "k_diagnostics.csv"))
  write_csv(stab$by_country, out_path(cfg, dataset, "role_stability_by_country.csv"))
  write_csv(stab$transition_matrix, out_path(cfg, dataset, "role_transitions.csv"))
  saveRDS(list(similarity = sim, hclust = cl$hclust), out_path(cfg, dataset, "rege_similarity.rds"))
  msg("  role persistence between consecutive periods: ",
      round(100 * stab$overall_persistence, 1), "%")

  list(roles = roles, profile = profile, labels = labels, diagnostics = cl$diagnostics,
       stability = stab, similarity = sim, hclust = cl$hclust)
}
