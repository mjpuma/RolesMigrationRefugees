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
# =============================================================================

node_id <- function(iso3, period) paste(iso3, period, sep = "_")

build_stacked_network <- function(flows, thresholds = 1) {
  base <- flows |> filter(flow >= min(thresholds))
  nodes <- sort(unique(c(node_id(base$origin, base$period),
                         node_id(base$destination, base$period))))
  n <- length(nodes)
  A_list <- lapply(thresholds, function(t) {
    e <- flows |> filter(flow >= t)
    A <- matrix(0, n, n, dimnames = list(nodes, nodes))
    A[cbind(match(node_id(e$origin, e$period), nodes),
            match(node_id(e$destination, e$period), nodes))] <- 1
    A
  })
  names(A_list) <- paste0(">=", thresholds)
  msg("  stacked network: ", n, " country-periods, ",
      sum(A_list[[1]]), " ties at threshold ", min(thresholds))
  A_list
}

# Within-cluster sum of squares computed from a distance matrix
wcss_from_dist <- function(D, cl) {
  sum(vapply(unique(cl), function(k) {
    idx <- which(cl == k)
    sum(D[idx, idx]^2) / (2 * length(idx))
  }, numeric(1)))
}

cluster_roles <- function(D, cfg_roles) {
  if (cfg_roles$cluster_on == "legacy_profiles") {
    # Reproduces the original scripts: factoextra::hcut(x = matrix) treats the
    # distance matrix as data, takes Euclidean distances between its rows, and
    # (because 'method =' is not an hcut argument) uses the default Ward.D2.
    d <- dist(D)
    hc <- hclust(d, method = "ward.D2")
  } else {
    d <- as.dist(D)
    hc <- hclust(d, method = cfg_roles$linkage)
  }
  ks <- seq(cfg_roles$k_range[1], cfg_roles$k_range[2])
  Dm <- as.matrix(d)
  diag_tbl <- map_dfr(ks, function(k) {
    cl <- cutree(hc, k)
    tibble(k = k,
           avg_silhouette = mean(cluster::silhouette(cl, d)[, "sil_width"]),
           wcss = wcss_from_dist(Dm, cl),
           smallest_cluster = min(table(cl)))
  })
  k <- if (identical(cfg_roles$k, "auto")) {
    diag_tbl$k[which.max(diag_tbl$avg_silhouette)]
  } else as.integer(cfg_roles$k)
  list(hclust = hc, k = k, role = cutree(hc, k), diagnostics = diag_tbl)
}

# Per country-period network statistics on the weighted graph (threshold 1)
node_metrics <- function(flows, threshold = 1) {
  e <- flows |> filter(flow >= threshold)
  per <- split(e, e$period)
  imap_dfr(per, function(df, p) {
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
              across(c(in_degree, out_degree, in_strength, out_strength, betweenness),
                     \(x) mean(x, na.rm = TRUE)),
              share_pure_sender = mean(pure_sender),
              share_pure_receiver = mean(pure_receiver),
              median_net_flow = median(net_flow),
              .groups = "drop")
}

# Automatic, human-readable role labels from the profile (edit freely)
label_roles <- function(profile) {
  total <- profile$in_strength + profile$out_strength
  two_way <- profile$share_pure_sender <= 0.9 & profile$share_pure_receiver <= 0.9
  lvl <- rep(NA_character_, nrow(profile))
  r <- rank(total[two_way], ties.method = "first")
  lvl[two_way] <- c("low", "medium", "high")[ceiling(3 * r / length(r))]
  balance <- ifelse(profile$median_net_flow > 0, "net receiving", "net sending")
  profile |>
    mutate(label = case_when(
      share_pure_sender > 0.9 ~ "sender only",
      share_pure_receiver > 0.9 ~ "receiver only",
      TRUE ~ paste0("two-way, ", lvl, " volume, ", balance)
    )) |>
    select(role, label)
}

# How often does each country keep its role from one period to the next?
role_stability <- function(roles_df) {
  trans <- roles_df |>
    arrange(iso3, period) |>
    group_by(iso3) |>
    mutate(prev_role = lag(role), prev_period = lag(period)) |>
    ungroup() |>
    filter(!is.na(prev_role))
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

run_roles <- function(cfg, dataset, flows) {
  rc <- cfg$roles
  # a dataset may override any roles setting (e.g. its own tie thresholds)
  rc <- modifyList(rc, cfg$flows[[dataset]]$roles %||% list())
  if (identical(rc$method, "valued")) {
    W <- build_stacked_valued(flows, min(unlist(rc$thresholds)))
    msg("  valued REGE (", rc$rege_iterations, " iterations) on ", nrow(W), " nodes ...")
    sim <- rege_valued(W, rc$rege_iterations, rc$value_transform %||% "log1p",
                       rc$no_self_match %||% FALSE)$similarity
  } else {
    A_list <- build_stacked_network(flows, unlist(rc$thresholds))
    msg("  REGE (", rc$rege_iterations, " iterations) on ", nrow(A_list[[1]]), " nodes ...")
    sim <- rege(A_list, rc$rege_iterations)$similarity
  }
  D <- 1 - sim
  cl <- cluster_roles(D, rc)
  msg("  k = ", cl$k, "; role sizes: ", paste(table(cl$role), collapse = " / "))

  roles <- tibble(node = names(cl$role), role = unname(cl$role)) |>
    separate(node, c("iso3", "period"), sep = "_", convert = TRUE)
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
