# =============================================================================
# migrant_rege_luczkovich.R — migrant roles with the regular-equivalence method
# of Luczkovich, Borgatti, Johnson & Everett (2003)
#
#   Rscript analysis/migrant_rege_luczkovich.R [config.yml]
#   Environment: LUCZ_PARTS   any of: period stacked   (default: period)
#                LUCZ_METHODS Abel & Cohen estimates (default: da_pb_closed,da_min_closed)
#                LUCZ_CORES   threads for REGE (default 2)
#
# Luczkovich et al. assign trophic roles in food webs, applying REGE to valued
# carbon flows. Their steps, and how each is carried over to migration:
#
#   network      directed, valued: x_ik = carbon flowing from k to i
#                here: x_ik = gross migrant flow from country i to country k
#                (all flows, as estimated by Abel and Cohen, no threshold)
#   REGE         valued REGE, iterated until the coefficients stop changing
#                (src/rege_ratio.cpp; each tie is matched by a tie of similar
#                size to an equivalent partner, scored by min/max of the sizes)
#   clustering   hierarchical clustering of the REGE coefficients (average linkage)
#   number of    R-squared: the coefficients of all pairs regressed on a dummy
#   classes      that is 1 when the pair is in the same class. k is read from
#                the scree plot of R-squared against k (paper Figs. 4 and 8)
#   display      nonmetric MDS of the coefficients; image graph of the classes
#
# The paper analyses one food web. Two versions are run here:
#   part "period"   one network per five-year period, as in the paper
#   part "stacked"  all periods in one network of country-periods (as in the
#                   pipeline), so classes are comparable over time; slow
#
# For comparison, part "period" also runs White's REGGE as implemented in UCINET
# (src/rege_valued.cpp), in which ties are weighted by their size.
#
# Writes to <output_dir>/undesa_luczkovich/:
#   scree_<version>.csv           R-squared and silhouette for k = 2..15
#   classes_<version>.csv         class of every country-period, with its flows
#   class_profiles_<version>.csv  size, partners and flows of each class
#   agreement.csv                 adjusted Rand index with other role definitions
#   L1_scree.png, L2_mds_<version>.png, L3_image_<version>.png
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(readr); library(ggplot2)
})
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
cfg_file <- commandArgs(trailingOnly = TRUE)
cfg <- read_config(if (length(cfg_file)) cfg_file[1] else "config.yml")
outdir <- dirname(out_path(cfg, "undesa_luczkovich", "x"))
parts <- strsplit(Sys.getenv("LUCZ_PARTS", "period"), " +")[[1]]
methods <- strsplit(Sys.getenv("LUCZ_METHODS", "da_pb_closed,da_min_closed"), ",")[[1]]
cores <- as.integer(Sys.getenv("LUCZ_CORES", "2"))
K <- 2:15

# compile the two C++ REGE functions (OpenMP if the compiler has it)
if (Sys.info()[["sysname"]] == "Linux") Sys.setenv(PKG_CXXFLAGS = "-fopenmp", PKG_LIBS = "-fopenmp")
Rcpp::sourceCpp("src/rege_ratio.cpp")
Rcpp::sourceCpp("src/rege_valued.cpp")

# gross flows of one Abel & Cohen estimate (no netting, no threshold)
read_gross <- function(method) {
  c2 <- cfg; c2$flows$undesa$value <- method; c2$flows$undesa$net <- FALSE
  read_flows(c2, "undesa") |> filter(origin != destination, flow > 0)
}

# R-squared of the paper: regress the REGE coefficient of every pair on a dummy
# for "same class"; with one dummy this is the squared correlation of the two
r2_partition <- function(S, cl) {
  u <- upper.tri(S); same <- outer(cl, cl, "==")[u]
  if (all(same) || !any(same)) return(NA_real_)
  cor(S[u], as.numeric(same))^2
}
scree <- function(S, hc) {
  D <- as.dist(1 - S)
  map_dfr(K, function(k) {
    cl <- cutree(hc, k)
    tibble(k = k, r2 = r2_partition(S, cl),
           silhouette = mean(cluster::silhouette(cl, D)[, "sil_width"]),
           smallest_class = min(table(cl)))
  })
}
# k at the peak of R-squared (the paper reads k off the scree plot; the plot is saved too)
pick_k <- function(sc) sc$k[which.max(sc$r2)]

# per-node flows, to describe the classes
node_flows <- function(fl) bind_rows(
  fl |> group_by(iso3 = origin, period) |> summarise(out_partners = n(), out_flow = sum(flow), .groups = "drop"),
  fl |> group_by(iso3 = destination, period) |> summarise(in_partners = n(), in_flow = sum(flow), .groups = "drop")) |>
  group_by(iso3, period) |> summarise(across(everything(), \(x) sum(x, na.rm = TRUE)), .groups = "drop")

# class profiles, with classes ordered by total flow (1 = largest)
profile <- function(cls) cls |> group_by(version, method, period_set, class) |>
  summarise(n = n(), out_partners = median(out_partners), in_partners = median(in_partners),
            out_flow = median(out_flow), in_flow = median(in_flow),
            share_net_sender = mean(out_flow > in_flow),
            examples = paste(head(iso3[order(-(out_flow + in_flow))], 8), collapse = " "), .groups = "drop")

res_scree <- list(); res_cls <- list(); sims <- list()

# ---- part "period": one network per period, as in the paper --------------------
if ("period" %in% parts) for (meth in methods) {
  fl <- read_gross(meth); nf <- node_flows(fl)
  for (p in sort(unique(fl$period))) {
    W <- build_stacked_valued(fl |> filter(period == p), threshold = 0)
    for (alg in c("ratio", "REGGE")) {
      t0 <- Sys.time()
      r <- if (alg == "ratio") rege_ratio_cpp(W, 60, 1e-3, cores) else rege_valued_cpp(W, 10, FALSE)
      S <- r$similarity; dimnames(S) <- dimnames(W)
      hc <- hclust(as.dist(1 - S), method = "average")
      sc <- scree(S, hc) |> mutate(version = paste0("period, ", alg), method = meth, period_set = as.character(p),
                                   iterations = length(r$change), last_change = tail(r$change, 1))
      k <- pick_k(sc)
      msg("  ", meth, " ", p, " ", alg, ": ", length(r$change), " iterations (last change ", signif(tail(r$change, 1), 2),
          "), k = ", k, " in ", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), " s")
      res_scree[[paste(meth, p, alg)]] <- sc
      ids <- rownames(S)
      res_cls[[paste(meth, p, alg)]] <- tibble(node = ids, class = cutree(hc, k)) |>
        tidyr::separate(node, c("iso3", "period"), sep = "_", convert = TRUE) |>
        left_join(nf, by = c("iso3", "period")) |>
        mutate(version = paste0("period, ", alg), method = meth, period_set = as.character(p), k = k)
      if (meth == methods[1] && alg == "ratio") sims[[as.character(p)]] <- S
    }
  }
}

# ---- part "stacked": all periods in one network ----------------------------------
if ("stacked" %in% parts) for (meth in methods[1]) {
  fl <- read_gross(meth); nf <- node_flows(fl)
  W <- build_stacked_valued(fl, threshold = 0)
  t0 <- Sys.time()
  r <- rege_ratio_cpp(W, 60, 1e-3, cores)
  S <- r$similarity; dimnames(S) <- dimnames(W)
  saveRDS(list(S = S, change = r$change), file.path(outdir, "stacked_similarity.rds"))
  hc <- hclust(as.dist(1 - S), method = "average")
  sc <- scree(S, hc) |> mutate(version = "stacked, ratio", method = meth, period_set = "all",
                               iterations = length(r$change), last_change = tail(r$change, 1))
  k <- pick_k(sc)
  msg("  stacked ", meth, ": ", length(r$change), " iterations, k = ", k, " in ",
      round(as.numeric(difftime(Sys.time(), t0, units = "mins"))), " min")
  res_scree[["stacked"]] <- sc
  res_cls[["stacked"]] <- tibble(node = rownames(S), class = cutree(hc, k)) |>
    tidyr::separate(node, c("iso3", "period"), sep = "_", convert = TRUE) |>
    left_join(nf, by = c("iso3", "period")) |>
    mutate(version = "stacked, ratio", method = meth, period_set = "all", k = k)
  sims[["stacked"]] <- S
}

tag <- paste(parts, collapse = "_")
sc_all <- bind_rows(res_scree); cls <- bind_rows(res_cls)
# order classes by total flow so that class 1 is the largest movers in every run
cls <- cls |> group_by(version, method, period_set, class) |> mutate(tot = median(out_flow + in_flow)) |>
  group_by(version, method, period_set) |> mutate(class = match(class, unique(class[order(-tot)]))) |> ungroup() |> select(-tot)
write_csv(sc_all, file.path(outdir, paste0("scree_", tag, ".csv")))
write_csv(cls, file.path(outdir, paste0("classes_", tag, ".csv")))
write_csv(profile(cls), file.path(outdir, paste0("class_profiles_", tag, ".csv")))

# ---- agreement with other role definitions and with flow size -------------------
bb <- read_csv(out_path(cfg, "undesa_backbone", "backbone_roles.csv"), show_col_types = FALSE) |> select(iso3, period, backbone = role)
ag <- cls |> left_join(bb, by = c("iso3", "period")) |>
  group_by(version, method, period_set) |>
  summarise(k = first(k), n = n(),
            ari_backbone = adjusted_rand(class[!is.na(backbone)], backbone[!is.na(backbone)]),
            # how well do partners and volume alone predict the classes? (the original roles: 98% from partners)
            r2_partners = summary(lm(as.numeric(factor(class)) ~ log1p(out_partners + in_partners)))$r.squared,
            ari_volume_tiers = adjusted_rand(class, dplyr::ntile(out_flow + in_flow, first(k))), .groups = "drop")
if (length(methods) > 1 && "period" %in% parts) {
  w <- cls |> filter(grepl("^period", version)) |> select(version, method, period_set, iso3, class) |>
    pivot_wider(names_from = method, values_from = class)
  am <- w |> filter(!is.na(.data[[methods[1]]]), !is.na(.data[[methods[2]]])) |> group_by(version, period_set) |>
    summarise(ari_methods = adjusted_rand(.data[[methods[1]]], .data[[methods[2]]]), .groups = "drop")
  ag <- ag |> left_join(am, by = c("version", "period_set"))
}
write_csv(ag, file.path(outdir, paste0("agreement_", tag, ".csv")))
print(as.data.frame(ag |> mutate(across(where(is.double), \(x) round(x, 3)))))

# ---- figures ------------------------------------------------------------------------
theme_set(theme_minimal(9) + theme(legend.position = "bottom", panel.grid.minor = element_blank()))
p1 <- ggplot(sc_all |> filter(method == methods[1]), aes(k, r2, colour = period_set)) + geom_line() + geom_point(size = 0.8) +
  facet_wrap(~version) + labs(x = "Number of classes", y = "R-squared (pairs in the same class)", colour = "Period")
ggsave(file.path(outdir, paste0("L1_scree_", tag, ".png")), p1, width = 8, height = 4, dpi = 200, bg = "white")
for (nm in names(sims)) {
  S <- sims[[nm]]; ver <- if (nm == "stacked") "stacked, ratio" else "period, ratio"
  set.seed(cfg$seed)
  m <- MASS::isoMDS(as.dist(pmax(1 - S, 1e-6)), k = 2, trace = FALSE)
  d <- tibble(node = rownames(S), x = m$points[, 1], y = m$points[, 2]) |>
    tidyr::separate(node, c("iso3", "period"), sep = "_", convert = TRUE) |>
    left_join(cls |> filter(version == ver, method == methods[1], period_set == (if (nm == "stacked") "all" else nm)), by = c("iso3", "period"))
  pm <- ggplot(d, aes(x, y, colour = factor(class))) + geom_point(size = 1, alpha = 0.8) +
    ggrepel::geom_text_repel(data = d |> group_by(class) |> slice_max(out_flow + in_flow, n = 4), aes(label = iso3), size = 2, show.legend = FALSE) +
    labs(title = paste0("Nonmetric MDS of REGE coefficients, ", nm, " (stress ", round(m$stress, 1), ")"), x = NULL, y = NULL, colour = "Class")
  ggsave(file.path(outdir, paste0("L2_mds_", nm, ".png")), pm, width = 6, height = 5, dpi = 200, bg = "white")
}
# image graph: share of possible ties present, and total flow, from each class to each
# other class (paper Figs. 6B and 9B); one panel per period for the main method
img <- map_dfr(unique(cls$method), function(meth) {
  fl <- read_gross(meth)
  cls |> filter(method == meth, version %in% c("period, ratio", "stacked, ratio")) |>
    group_by(version, period_set) |> group_modify(function(cc, key) {
      e <- fl |> inner_join(cc |> select(origin = iso3, period, from = class), by = c("origin", "period")) |>
        inner_join(cc |> select(destination = iso3, period, to = class), by = c("destination", "period"))
      sz <- cc |> count(period, class)
      pairs <- sz |> inner_join(sz, by = "period", suffix = c("_f", "_t"), relationship = "many-to-many") |>
        mutate(possible = n_f * n_t - ifelse(class_f == class_t, n_f, 0)) |>
        group_by(from = class_f, to = class_t) |> summarise(possible = sum(possible), .groups = "drop")
      e |> group_by(from, to) |> summarise(ties = n(), flow = sum(flow), .groups = "drop") |>
        right_join(pairs, by = c("from", "to")) |> mutate(across(c(ties, flow), \(x) coalesce(x, 0)), density = ties / possible)
    }) |> ungroup() |> mutate(method = meth)
})
write_csv(img, file.path(outdir, paste0("image_graph_", tag, ".csv")))
pi_ <- ggplot(img |> filter(method == methods[1]), aes(factor(to), factor(from), fill = density)) + geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", density)), size = 2) + scale_fill_gradient(low = "white", high = "#2D6FA8", limits = c(0, 1)) +
  facet_wrap(~paste(version, period_set)) + labs(x = "To class", y = "From class", fill = "Share of possible ties")
ggsave(file.path(outdir, paste0("L3_image_", tag, ".png")), pi_, width = 9, height = 6, dpi = 200, bg = "white")
msg("Done: ", outdir)
