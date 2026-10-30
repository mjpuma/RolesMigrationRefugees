# =============================================================================
# appendix_figures.R — network diagrams and role-separation plots
#
#   Rscript appendix_figures.R      (after run_pipeline.R)
#
# For each dataset (undesa, unhcr) writes to output/<dataset>/figures/appendix/:
#   A_separation_scatter.png  in- vs out-flow per country-period, coloured by role
#   B_rege_space.png          2-D map (classical MDS) of REGE distances, by role
#   C_role_image.png          role-level network: total flows between roles
#   D_country_network_<yr>.png country network for one year (largest flows)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2)
  library(igraph); library(ggrepel)
})
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
cfg <- read_config("config.yml")

pal <- function(k) role_palette(k)
theme_set(theme_minimal(11) + theme(legend.position = "bottom", plot.title.position = "plot"))

network_year <- c(undesa = 2020, unhcr = 2015)   # year shown in diagram D
top_edges    <- c(undesa = 100,  unhcr = 150)    # strongest flows drawn in D

for (ds in cfg$run) {
  msg("Appendix figures: ", ds)
  outdir <- out_path(cfg, ds, "figures", "appendix", "x"); outdir <- dirname(outdir)
  roles <- read_csv(out_path(cfg, ds, "roles_by_country_period.csv"), show_col_types = FALSE) |>
    mutate(role_lab = factor(paste0(role, ": ", label)), name = coalesce(country, iso3))
  k <- n_distinct(roles$role)
  cols <- setNames(pal(k), levels(roles$role_lab))
  flows <- read_flows(cfg, ds) |> filter(flow >= 1)
  title_ds <- if (ds == "unhcr") "Refugees (UNHCR)" else "Migrants (UN DESA / Abel & Cohen)"

  # ---- A. separation scatter ----------------------------------------------------
  last <- max(roles$period)
  lab <- roles |> filter(period == last) |>
    group_by(role) |> slice_max(in_strength + out_strength, n = 3, with_ties = FALSE) |> ungroup()
  pA <- ggplot(roles, aes(log10(out_strength + 1), log10(in_strength + 1), colour = role_lab)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey75", linetype = 2) +
    geom_point(alpha = 0.55, size = 1.6) +
    geom_text_repel(data = lab, aes(label = paste(name, period)), size = 2.8,
                    show.legend = FALSE, max.overlaps = 30, min.segment.length = 0) +
    scale_colour_manual(values = cols) +
    labs(title = paste0(title_ds, ": how the roles separate"),
         subtitle = "Each point is one country in one period. Axes: log10(people + 1). Dashed line = equal in- and outflow.",
         x = "Outflow (people leaving)", y = "Inflow (people arriving)", colour = "Role") +
    guides(colour = guide_legend(ncol = 2, override.aes = list(size = 3, alpha = 1)))
  ggsave(file.path(outdir, "A_separation_scatter.png"), pA, width = 8, height = 7, dpi = 200, bg = "white")

  # ---- B. REGE space (classical MDS of 1 - REGE) ---------------------------------
  sim <- readRDS(out_path(cfg, ds, "rege_similarity.rds"))$similarity
  mds <- cmdscale(as.dist(1 - sim), k = 2, eig = TRUE)
  share <- round(100 * mds$eig[1:2] / sum(pmax(mds$eig, 0)), 1)
  B <- tibble(node = rownames(sim), d1 = mds$points[, 1], d2 = mds$points[, 2]) |>
    separate(node, c("iso3", "period"), sep = "_", convert = TRUE) |>
    left_join(roles |> select(iso3, period, role_lab, name), by = c("iso3", "period"))
  pB <- ggplot(B, aes(d1, d2, colour = role_lab)) +
    geom_point(alpha = 0.6, size = 1.6, position = position_jitter(width = 0.002, height = 0.002, seed = 1)) +
    scale_colour_manual(values = cols) +
    labs(title = paste0(title_ds, ": countries in regular-equivalence space"),
         subtitle = sprintf("Classical MDS of 1 - REGE similarity (dimension 1: %.1f%%, dimension 2: %.1f%% of variation).\nClose points = similar roles.", share[1], share[2]),
         x = "MDS dimension 1", y = "MDS dimension 2", colour = "Role") +
    guides(colour = guide_legend(ncol = 2, override.aes = list(size = 3, alpha = 1)))
  ggsave(file.path(outdir, "B_rege_space.png"), pB, width = 8, height = 7, dpi = 200, bg = "white")

  # ---- C. role-level network (blockmodel image) ---------------------------------
  rl <- roles |> select(iso3, period, role, role_lab)
  E <- flows |>
    inner_join(rl, by = c("origin" = "iso3", "period")) |>
    inner_join(rl, by = c("destination" = "iso3", "period"), suffix = c("_from", "_to")) |>
    group_by(role_lab_from, role_lab_to) |>
    summarise(people = sum(flow), links = n(), .groups = "drop") |>
    mutate(share = people / sum(people))
  write_csv(E, file.path(outdir, "C_role_image_flows.csv"))
  sizes <- roles |> count(role_lab, name = "n")
  g <- graph_from_data_frame(E |> transmute(from = as.character(role_lab_from),
                                            to = as.character(role_lab_to), share),
                             vertices = sizes |> transmute(name = as.character(role_lab), n))
  png(file.path(outdir, "C_role_image.png"), width = 1800, height = 1500, res = 200)
  par(mar = c(1, 1, 3, 1))
  set.seed(2)
  lay <- layout_in_circle(g)
  w <- E(g)$share
  plot(g, layout = lay,
       vertex.size = 12 + 28 * sqrt(V(g)$n / max(V(g)$n)),
       vertex.color = cols[V(g)$name], vertex.frame.color = "white",
       vertex.label = sub(":.*", "", V(g)$name), vertex.label.color = "white",
       vertex.label.font = 2, vertex.label.cex = 1.1,
       edge.width = 0.5 + 18 * w, edge.arrow.size = 0.4 + 0.8 * w,
       edge.color = sapply(0.35 + 0.6 * w / max(w), function(a) adjustcolor("grey30", a)),
       edge.curved = 0.25, edge.loop.angle = pi / 4,
       edge.label = ifelse(w >= 0.02, sprintf("%.0f%%", 100 * w), NA),
       edge.label.cex = 0.8, edge.label.color = "black")
  title(paste0(title_ds, ": flows between roles (all periods)"), cex.main = 1)
  legend("topright", legend = names(cols), fill = cols, border = NA, bty = "n", cex = 0.7)
  mtext("Node size = number of country-periods; arrow width and label = share of all people moving (labels shown for >= 2%).",
        side = 1, line = -0.5, cex = 0.65)
  dev.off()

  # ---- C2. role-to-role destination shares (heatmap) ---------------------------
  H <- E |> group_by(role_lab_from) |> mutate(row_share = people / sum(people)) |> ungroup() |>
    mutate(from = factor(sub(":.*", "", role_lab_from)), to = factor(sub(":.*", "", role_lab_to)))
  pC2 <- ggplot(H, aes(to, forcats::fct_rev(from), fill = row_share)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = sprintf("%.0f%%", 100 * row_share),
                  colour = row_share > 0.5), size = 3.5, show.legend = FALSE) +
    scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = "white")) +
    scale_fill_gradient(low = "#f2f5f9", high = "#1f3a5f", labels = scales::percent, limits = c(0, 1)) +
    labs(title = paste0(title_ds, ": where each role's outflows go"),
         subtitle = "Rows: role of the origin country. Columns: role of the destination.\nCells: share of the origin role's total outflow (all periods).",
         x = "Destination role", y = "Origin role", fill = "Share") +
    theme(panel.grid = element_blank())
  ggsave(file.path(outdir, "C2_role_destination_shares.png"), pC2, width = 7.5, height = 5.5, dpi = 200, bg = "white")

  # ---- D. country network for one year -------------------------------------------
  yr <- network_year[[ds]]
  # each country's 2 largest outflows (so every role appears), plus the overall largest flows
  fy <- bind_rows(
    flows |> filter(period == yr) |> group_by(origin) |> slice_max(flow, n = 2, with_ties = FALSE) |> ungroup(),
    flows |> filter(period == yr) |> slice_max(flow, n = top_edges[[ds]])) |> distinct()
  ry <- roles |> filter(period == yr)
  gy <- graph_from_data_frame(fy |> select(origin, destination, flow), directed = TRUE)
  vr <- ry[match(V(gy)$name, ry$iso3), ]
  strength_tot <- vr$in_strength + vr$out_strength
  top_lab <- rank(-strength_tot, ties.method = "first") <= 20 |
    ave(-strength_tot, vr$role, FUN = function(x) rank(x, ties.method = "first")) <= 3
  png(file.path(outdir, paste0("D_country_network_", yr, ".png")), width = 2200, height = 2000, res = 200)
  par(mar = c(1, 1, 3, 1))
  set.seed(3)
  lay <- layout_with_fr(gy, weights = log10(E(gy)$flow))
  plot(gy, layout = lay,
       vertex.size = 1.5 + 7 * sqrt(strength_tot / max(strength_tot, na.rm = TRUE)),
       vertex.color = cols[as.character(vr$role_lab)], vertex.frame.color = NA,
       vertex.label = ifelse(top_lab, vr$name, NA), vertex.label.cex = 0.6,
       vertex.label.color = "black", vertex.label.dist = 0.6,
       edge.width = 0.3 + 3 * (log10(E(gy)$flow) - min(log10(E(gy)$flow))) /
                    diff(range(log10(E(gy)$flow))),
       edge.arrow.size = 0.15, edge.color = adjustcolor("grey40", 0.35), edge.curved = 0.15)
  n_all <- sum(flows$period == yr)
  shows_all <- nrow(fy) == n_all
  title(if (shows_all) sprintf("%s: all %d flows recorded in %d", title_ds, n_all, yr)
        else sprintf("%s: main flows in %d", title_ds, yr), cex.main = 1)
  present <- names(cols)[names(cols) %in% as.character(vr$role_lab)]
  legend("bottomleft", legend = present, fill = cols[present], border = NA, bty = "n", cex = 0.7)
  mtext(paste(if (shows_all) sprintf("Edges: all %d flows recorded in %d.", n_all, yr)
              else sprintf("Edges: each country's 2 largest outflows plus the %d largest flows overall (%d of %d flows).", top_edges[[ds]], nrow(fy), n_all),
              "Node colour = role; size = total in + out flow; labels = 20 largest countries + 3 largest per role."),
        side = 1, line = -0.5, cex = 0.65)
  dev.off()
}
msg("Done: output/<dataset>/figures/appendix/")
