# =============================================================================
# refugee_rege_weighted.R — refugee roles from flow sizes (valued REGE), and how
# many REGE iterations are needed
#
#   Rscript analysis/refugee_rege_weighted.R [config.yml]     (after run_pipeline.R)
#
# The pipeline builds refugee roles from binary ties (at least one refugee). Here
# the same stacked network of country-years is analysed with the flow sizes,
# following Luczkovich, Borgatti, Johnson & Everett (2003):
#
#   ratio      REGE as described in that paper (src/rege_ratio.cpp): each tie is
#              matched by a tie of similar size, min/max of the two, to an
#              equivalent partner; every tie counts once
#   REGGE      White's REGGE as implemented in UCINET (src/rege_valued.cpp): ties
#              weighted by size, in- and out-ties to a partner matched together
#   REGGE log  the same on log(1 + flow), the version in the earlier code review
#   binary     the pipeline's REGE on 0/1 ties, for reference
#
# Iterations. REGE starts with every pair fully equivalent and lowers the
# coefficients at each iteration; the limit is the maximal regular equivalence.
# UCINET's default is 3 iterations (also used in the pipeline); the paper iterates
# until nothing changes. Each method is run up to 200 iterations, and at selected
# iterations the clustering is repeated, to see when the roles stop changing.
#
# Clustering: average linkage on 1 - REGE. Number of roles: the R-squared of the
# paper (pair coefficients regressed on a same-class dummy) and the silhouette.
#
# Writes to <output_dir>/unhcr_weighted/:
#   convergence.csv        largest change in any coefficient, every iteration
#   iterations.csv         at selected iterations: k by R-squared and silhouette,
#                          agreement with the roles after the last iteration and
#                          with the binary roles of the pipeline
#   scree.csv              R-squared and silhouette for k = 2..12 (last iteration)
#   roles_weighted.csv     role of every country-year (k = 3 and k by R-squared)
#   role_profiles.csv      size, flows and the binary roles inside each weighted role
#   W1_convergence.png, W2_iterations.png, W3_scree.png
# =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(purrr); library(readr); library(ggplot2) })
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
cfg_file <- commandArgs(trailingOnly = TRUE)
cfg <- read_config(if (length(cfg_file)) cfg_file[1] else "config.yml")
outdir <- dirname(out_path(cfg, "unhcr_weighted", "x"))
if (Sys.info()[["sysname"]] == "Linux") Sys.setenv(PKG_CXXFLAGS = "-fopenmp", PKG_LIBS = "-fopenmp")
Rcpp::sourceCpp("src/rege_ratio.cpp")
Rcpp::sourceCpp("src/rege_valued.cpp")
max_iter <- 200
keep <- c(1, 2, 3, 4, 5, 7, 10, 15, 20, 30, 50, 75, 100, 150, 200)
K <- 2:12

# the stacked country-year network of the pipeline (ties of at least one refugee)
fl <- read_flows(cfg, "unhcr") |> filter(origin != destination, flow >= 1)
W <- build_stacked_valued(fl, threshold = 1)
nodes <- rownames(W)
binary <- read_csv(out_path(cfg, "unhcr", "roles_by_country_period.csv"), show_col_types = FALSE) |>
  transmute(node = node_id(iso3, period), binary = role)
b_role <- binary$binary[match(nodes, binary$node)]
stopifnot(!anyNA(b_role))
msg("Refugee network: ", nrow(W), " country-years, ", sum(W > 0), " ties")

r2_partition <- function(S, cl) { u <- upper.tri(S); cor(S[u], as.numeric(outer(cl, cl, "==")[u]))^2 }
cluster_all <- function(S) {
  D <- as.dist(1 - S); hc <- hclust(D, method = "average")
  sc <- map_dfr(K, function(k) { cl <- cutree(hc, k)
    tibble(k = k, r2 = r2_partition(S, cl), silhouette = mean(cluster::silhouette(cl, D)[, "sil_width"]), smallest = min(table(cl))) })
  list(hc = hc, scree = sc)
}

# ---- run every method, keeping the similarity at selected iterations ------------
A <- (W > 0) * 1
runs <- list(
  ratio       = function() rege_ratio_cpp(W, max_iter, 0, 2, as.integer(keep)),
  REGGE       = function() rege_valued_cpp(W, max_iter, FALSE, as.integer(keep)),
  `REGGE log` = function() rege_valued_cpp(log1p(W), max_iter, FALSE, as.integer(keep)))
conv <- list(); its <- list(); final <- list()
for (nm in names(runs)) {
  t0 <- Sys.time(); r <- suppressMessages(capture.output(res <- runs[[nm]]()))
  conv[[nm]] <- tibble(method = nm, iteration = seq_along(res$change), max_change = res$change)
  msg("  ", nm, ": ", length(res$change), " iterations in ", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), " s; last change ", signif(tail(res$change, 1), 2))
  final[[nm]] <- res$snapshots[[length(res$snapshots)]]
  its[[nm]] <- res$snapshots
}
# binary REGE (R implementation of the pipeline) at the same iterations
its$binary <- setNames(lapply(keep[keep <= 50], function(t) rege(A, t, verbose = FALSE)$similarity), keep[keep <= 50])
final$binary <- its$binary[[length(its$binary)]]
b_conv <- rege(A, 50, verbose = FALSE)$change
conv$binary <- tibble(method = "binary", iteration = seq_along(b_conv), max_change = b_conv)
conv <- bind_rows(conv); write_csv(conv, file.path(outdir, "convergence.csv"))

# ---- does the clustering still change with more iterations? ----------------------
fin_cl <- map(final, function(S) cluster_all(S))
it_tab <- imap_dfr(its, function(snaps, nm) imap_dfr(snaps, function(S, t) {
  cc <- cluster_all(S)
  tibble(method = nm, iteration = as.integer(t),
         k_r2 = cc$scree$k[which.max(cc$scree$r2)], k_silhouette = cc$scree$k[which.max(cc$scree$silhouette)],
         r2_k3 = cc$scree$r2[cc$scree$k == 3], silhouette_k3 = cc$scree$silhouette[cc$scree$k == 3],
         mean_similarity = mean(S[upper.tri(S)]),
         ari_last_k3 = adjusted_rand(cutree(cc$hc, 3), cutree(fin_cl[[nm]]$hc, 3)),
         ari_binary_roles_k3 = adjusted_rand(cutree(cc$hc, 3), b_role))
}))
write_csv(it_tab, file.path(outdir, "iterations.csv"))
print(as.data.frame(it_tab |> mutate(across(where(is.double), \(x) round(x, 3)))))

# ---- roles after the last iteration --------------------------------------------------
sc <- imap_dfr(fin_cl, \(cc, nm) cc$scree |> mutate(method = nm)); write_csv(sc, file.path(outdir, "scree.csv"))
nf <- bind_rows(fl |> group_by(iso3 = origin, period) |> summarise(out_partners = n(), out_flow = sum(flow), .groups = "drop"),
                fl |> group_by(iso3 = destination, period) |> summarise(in_partners = n(), in_flow = sum(flow), .groups = "drop")) |>
  group_by(iso3, period) |> summarise(across(everything(), \(x) sum(x, na.rm = TRUE)), .groups = "drop")
roles <- imap_dfr(fin_cl, function(cc, nm) {
  kr <- cc$scree$k[which.max(cc$scree$r2)]
  tibble(node = nodes, method = nm, k3 = cutree(cc$hc, 3), k_r2 = kr, role_k_r2 = cutree(cc$hc, kr), binary = b_role)
}) |> tidyr::separate(node, c("iso3", "period"), sep = "_", convert = TRUE) |> left_join(nf, by = c("iso3", "period"))
write_csv(roles, file.path(outdir, "roles_weighted.csv"))
prof <- roles |> group_by(method, k3) |>
  summarise(n = n(), sends = mean(out_partners > 0), receives = mean(in_partners > 0),
            median_out = median(out_flow), median_in = median(in_flow),
            binary_sender_only = mean(binary == 1), binary_two_way = mean(binary == 2), binary_receiver_only = mean(binary == 3),
            examples = paste(head(unique(iso3[order(-(out_flow + in_flow))]), 8), collapse = " "), .groups = "drop")
write_csv(prof, file.path(outdir, "role_profiles.csv"))
print(as.data.frame(prof |> mutate(across(where(is.double), \(x) round(x, 2)))))

# ---- figures --------------------------------------------------------------------------
theme_set(theme_minimal(9) + theme(legend.position = "bottom", panel.grid.minor = element_blank()))
g1 <- ggplot(conv |> filter(max_change > 0), aes(iteration, max_change, colour = method)) + geom_line() +
  geom_vline(xintercept = 3, linetype = 2, colour = "grey50") + scale_y_log10() +
  labs(x = "Iteration", y = "Largest change in any REGE coefficient (log scale)", colour = NULL,
       caption = "Dashed line: 3 iterations, the UCINET default used in the pipeline")
ggsave(file.path(outdir, "W1_convergence.png"), g1, width = 6, height = 4, dpi = 200, bg = "white")
g2 <- it_tab |> select(method, iteration, `Agreement with roles after the last iteration` = ari_last_k3,
                       `Agreement with the binary roles` = ari_binary_roles_k3) |>
  pivot_longer(-c(method, iteration)) |>
  ggplot(aes(iteration, value, colour = method)) + geom_line() + geom_point(size = 1) + scale_x_log10(breaks = keep) +
  facet_wrap(~name) + labs(x = "Iterations (log scale)", y = "Adjusted Rand index, three roles", colour = NULL)
ggsave(file.path(outdir, "W2_iterations.png"), g2, width = 8, height = 4, dpi = 200, bg = "white")
g3 <- sc |> pivot_longer(c(r2, silhouette)) |> ggplot(aes(k, value, colour = method)) + geom_line() + geom_point(size = 1) +
  facet_wrap(~name, scales = "free_y") + labs(x = "Number of roles", y = NULL, colour = NULL)
ggsave(file.path(outdir, "W3_scree.png"), g3, width = 8, height = 4, dpi = 200, bg = "white")
msg("Done: ", outdir)
