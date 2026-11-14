# =============================================================================
# displacement_weighted_roles.R — roles from flow sizes in the combined UNHCR
# data (refugee arrivals + asylum applications), one network per year
#
#   Rscript analysis/displacement_weighted_roles.R [config.yml] [dataset]
#   dataset: unhcr_displacement (default)
#
# With asylum applications, about 3,100 flows a year link about 186 countries,
# and three quarters of country-years both send and receive at least one person,
# so ties of "at least one person" no longer separate roles. Three ways of
# bringing flow size in are compared, each fitted to every year separately
# (as Luczkovich et al. 2003 fit one food web):
#
#   ratio REGE   src/rege_ratio.cpp, iterated to convergence, on
#                  raw     the flows themselves: 1,000 vs 10,000 people match 10%
#                  log     log(1 + flow): 1,000 vs 10,000 match 75%
#                  share   each flow as a share of the sender's total outflow
#                          (out-ties) and of the receiver's total inflow (in-ties)
#                          is not used here; "share" = flow / sender's outflow
#   backbone     keep a flow if it carries at least s of the origin's outflow or
#                at least s of the destination's inflow (s = 5, 10, 20%); binary
#                REGE on the kept ties
#   threshold    binary REGE on ties of at least 10, 100 or 1,000 people
#
# For every year and version: REGE at 1, 2, 3, 5, 10, 20 iterations and at
# convergence (largest change < 1e-4, at most 300 iterations); average-linkage
# clustering; number of roles by the R-squared scree of the paper and by
# silhouette. Reported: iterations needed, whether the roles at 3 iterations
# equal those at convergence, agreement between versions, and the profile of
# the three-role solution.
#
# Writes to <output_dir>/<dataset>_weighted/:
#   by_year.csv         one row per year and version
#   summary.csv         averages over years
#   agreement.csv       mean adjusted Rand index between versions (3 roles, converged)
#   roles_by_year.csv   3-role solution of every version (roles ordered by net flow)
#   profiles.csv        who is in each role
# =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(purrr); library(readr) })
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
a <- commandArgs(trailingOnly = TRUE)
cfg <- read_config(if (length(a)) a[1] else "config.yml")
ds <- if (length(a) > 1) a[2] else "unhcr_displacement"
outdir <- dirname(out_path(cfg, paste0(ds, "_weighted"), "x"))
if (Sys.info()[["sysname"]] == "Linux") Sys.setenv(PKG_CXXFLAGS = "-fopenmp", PKG_LIBS = "-fopenmp")
invisible(capture.output(Rcpp::sourceCpp("src/rege_ratio.cpp")))
KEEP <- c(1L, 2L, 3L, 5L, 10L, 20L); K <- 2:8

fl <- read_flows(cfg, ds) |> filter(origin != destination, flow >= 1)
years <- sort(unique(fl$period))

r2_partition <- function(S, cl) { u <- upper.tri(S); x <- as.numeric(outer(cl, cl, "==")[u]); if (sd(x) == 0 || sd(S[u]) == 0) NA_real_ else cor(S[u], x)^2 }
clus <- function(S) {
  D <- as.dist(1 - S); hc <- hclust(D, "average")
  sc <- map_dfr(K, function(k) { cl <- cutree(hc, k); tibble(k = k, r2 = r2_partition(S, cl), sil = mean(cluster::silhouette(cl, D)[, "sil_width"]), smallest = min(table(cl))) })
  list(hc = hc, sc = sc)
}
# binary REGE (R/rege.R) with snapshots, to convergence
rege_binary <- function(A, tol = 1e-4, max_iter = 300) {
  snaps <- list(); prev <- NULL
  for (t in c(KEEP, seq(30, max_iter, by = 10))) {
    S <- rege(A, t, verbose = FALSE)$similarity
    if (t %in% KEEP) snaps[[as.character(t)]] <- S
    if (!is.null(prev) && max(abs(S - prev)) < tol * 10) return(list(similarity = S, snapshots = snaps, iterations = t))
    prev <- S
  }
  list(similarity = S, snapshots = snaps, iterations = max_iter)
}
versions <- list(
  `ratio REGE, raw flows` = function(W) { invisible(capture.output(r <- rege_ratio_cpp(W, 300, 1e-4, 2, KEEP))); r$iterations <- length(r$change); r },
  `ratio REGE, log flows` = function(W) { invisible(capture.output(r <- rege_ratio_cpp(log1p(W), 300, 1e-4, 2, KEEP))); r$iterations <- length(r$change); r },
  `ratio REGE, outflow shares` = function(W) { X <- W / pmax(rowSums(W), 1e-12); invisible(capture.output(r <- rege_ratio_cpp(X, 300, 1e-4, 2, KEEP))); r$iterations <- length(r$change); r })
for (s in c(0.05, 0.10, 0.20)) local({ s0 <- s
  versions[[sprintf("backbone %d%%, binary REGE", round(100 * s0))]] <<- function(W) {
    keep <- (W / pmax(rowSums(W), 1e-12) >= s0) | (t(t(W) / pmax(colSums(W), 1e-12)) >= s0)
    rege_binary((W > 0 & keep) * 1) } })
for (th in c(10, 100, 1000)) local({ t0 <- th
  versions[[sprintf("ties of at least %s, binary REGE", format(t0, big.mark = ","))]] <<- function(W) rege_binary((W >= t0) * 1) })

by_year <- list(); roles <- list()
for (y in years) {
  W <- build_stacked_valued(fl |> filter(period == y), threshold = 1)
  iso <- sub("_.*", "", rownames(W)); out <- rowSums(W); inn <- colSums(W)
  for (v in names(versions)) {
    r <- versions[[v]](W); S <- r$similarity; cc <- clus(S)
    k_r2 <- cc$sc$k[which.max(cc$sc$r2)]; fin3 <- cutree(cc$hc, 3); fink <- cutree(cc$hc, k_r2)
    at <- function(t, k) { s <- r$snapshots[[as.character(t)]]; if (is.null(s)) NA_real_ else adjusted_rand(cutree(hclust(as.dist(1 - s), "average"), k), cutree(cc$hc, k)) }
    by_year[[paste(y, v)]] <- tibble(year = y, version = v, nodes = nrow(W), iterations = r$iterations,
      mean_similarity = mean(S[upper.tri(S)]), k_r2 = k_r2, max_r2 = max(cc$sc$r2, na.rm = TRUE), k_sil = cc$sc$k[which.max(cc$sc$sil)],
      r2_k3 = cc$sc$r2[cc$sc$k == 3], sil_k3 = cc$sc$sil[cc$sc$k == 3], smallest_k3 = cc$sc$smallest[cc$sc$k == 3],
      ari_iter3_vs_converged_k3 = at(3, 3), ari_iter10_vs_converged_k3 = at(10, 3), ari_iter3_vs_converged_kr2 = at(3, k_r2))
    # order the three roles by median net flow (1 = most net sending)
    net <- tapply(inn - out, fin3, median); fin3 <- match(fin3, as.integer(names(sort(net))))
    roles[[paste(y, v)]] <- tibble(year = y, version = v, iso3 = iso, role = fin3, out_flow = out, in_flow = inn,
                                   out_partners = rowSums(W > 0), in_partners = colSums(W > 0))
  }
  msg("  ", y, " done")
}
by_year <- bind_rows(by_year); roles <- bind_rows(roles)
write_csv(by_year, file.path(outdir, "by_year.csv")); write_csv(roles, file.path(outdir, "roles_by_year.csv"))
summ <- by_year |> group_by(version) |>
  summarise(iterations = median(iterations), max_iterations = max(iterations), mean_similarity = mean(mean_similarity),
            k_r2_range = paste(range(k_r2), collapse = "-"), r2_at_k_r2 = mean(max_r2), r2_k3 = mean(r2_k3), silhouette_k3 = mean(sil_k3),
            years_with_singleton_role_k3 = sum(smallest_k3 <= 2),
            ari_iter3_vs_converged = mean(ari_iter3_vs_converged_k3, na.rm = TRUE), ari_iter10_vs_converged = mean(ari_iter10_vs_converged_k3, na.rm = TRUE), .groups = "drop")
write_csv(summ, file.path(outdir, "summary.csv"))
ag <- roles |> select(year, version, iso3, role) |> inner_join(roles |> select(year, version2 = version, iso3, role2 = role), by = c("year", "iso3"), relationship = "many-to-many") |>
  filter(version < version2) |> group_by(version, version2, year) |> summarise(ari = adjusted_rand(role, role2), .groups = "drop") |>
  group_by(version, version2) |> summarise(mean_ari = mean(ari), min_ari = min(ari), .groups = "drop")
write_csv(ag, file.path(outdir, "agreement.csv"))
prof <- roles |> group_by(version, role) |>
  summarise(share_of_country_years = n() / n_distinct(paste(year)) , n = n(), sends = mean(out_partners > 0), receives = mean(in_partners > 0),
            median_out = median(out_flow), median_in = median(in_flow), median_out_partners = median(out_partners), median_in_partners = median(in_partners),
            examples = paste(names(sort(table(iso3), decreasing = TRUE))[1:8], collapse = " "), .groups = "drop")
# year-to-year persistence of the role of a country
pers <- roles |> arrange(version, iso3, year) |> group_by(version, iso3) |> mutate(same = role == lag(role) & year == lag(year) + 1) |>
  group_by(version) |> summarise(persistence = mean(same, na.rm = TRUE), .groups = "drop")
summ <- summ |> left_join(pers, by = "version"); write_csv(summ, file.path(outdir, "summary.csv"))
write_csv(prof, file.path(outdir, "profiles.csv"))
options(width = 250)
print(as.data.frame(summ |> mutate(across(where(is.double), \(x) round(x, 2)))))
print(as.data.frame(prof |> mutate(across(where(is.double), \(x) round(x, 2)))))
print(as.data.frame(ag |> mutate(across(where(is.double), \(x) round(x, 2)))))
msg("Done: ", outdir)
