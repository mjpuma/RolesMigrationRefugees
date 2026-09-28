# =============================================================================
# figures.R — standard figures for any flow dataset
# =============================================================================

role_palette <- function(k) {
  base <- c("#2a78b3", "#d9812b", "#3c9a5f", "#b34b4b", "#7b5ea7", "#8c7a5b",
            "#c265a6", "#6b6b6b", "#a8a232", "#2fa3b5")
  rep(base, length.out = k)
}

plot_k_diagnostics <- function(diag, dataset) {
  diag |>
    pivot_longer(c(avg_silhouette, wcss)) |>
    ggplot(aes(k, value)) + geom_line() + geom_point() +
    facet_wrap(~ name, scales = "free_y", ncol = 1,
               labeller = as_labeller(c(avg_silhouette = "Average silhouette (higher = better)",
                                        wcss = "Within-cluster sum of squares (elbow)"))) +
    labs(title = paste0(dataset, ": choosing the number of roles"), x = "k", y = NULL) +
    theme_minimal(11)
}

plot_strength_by_role <- function(roles, dataset) {
  roles |>
    mutate(role = paste0(role, ": ", label)) |>
    group_by(period, role) |>
    summarise(`In-strength` = mean(in_strength), `Out-strength` = mean(out_strength),
              .groups = "drop") |>
    pivot_longer(c(`In-strength`, `Out-strength`)) |>
    ggplot(aes(period, value + 1, colour = role)) +
    geom_line(linewidth = 0.8) + geom_point() +
    scale_y_log10(labels = scales::label_comma()) +
    scale_colour_manual(values = role_palette(n_distinct(roles$role))) +
    facet_wrap(~ name) +
    labs(title = paste0(dataset, ": mean flow volume by role"), x = NULL,
         y = "People (log scale)", colour = "Role") +
    guides(colour = guide_legend(ncol = 2)) +
    theme_minimal(11) + theme(legend.position = "bottom")
}

plot_role_timeline <- function(roles, dataset, area = NULL) {
  d <- roles
  if (!is.null(area)) d <- d |> filter(.data$area == !!area)
  d |>
    mutate(role = factor(paste0(role, ": ", label)),
           country = forcats::fct_reorder(coalesce(country, iso3), as.numeric(role), .fun = mean)) |>
    ggplot(aes(period, country, fill = role)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    scale_fill_manual(values = role_palette(n_distinct(roles$role))) +
    labs(title = paste0(dataset, ": role of each country over time",
                        if (!is.null(area)) paste0(" (", area, ")")),
         x = NULL, y = NULL, fill = "Role") +
    guides(fill = guide_legend(ncol = 2)) +
    theme_minimal(9) + theme(legend.position = "bottom", panel.grid = element_blank())
}

plot_odds_ratios <- function(bin, dataset) {
  bin |>
    filter(!grepl("Intercept|period", term), !unstable) |>
    mutate(sig = p.value < 0.05) |>
    ggplot(aes(odds_ratio, term, colour = sig)) +
    geom_vline(xintercept = 1, colour = "grey60") +
    geom_errorbarh(aes(xmin = or_low, xmax = or_high), height = 0) +
    geom_point() +
    scale_x_log10() +
    scale_colour_manual(values = c(`FALSE` = "grey55", `TRUE` = "#b34b4b"),
                        labels = c("p >= 0.05", "p < 0.05")) +
    facet_wrap(~ paste("Role", role)) +
    labs(title = paste0(dataset, ": drivers of role membership"),
         x = "Odds ratio per 1 SD (95% CI, log scale)", y = NULL, colour = NULL) +
    theme_minimal(10) + theme(legend.position = "bottom")
}

save_figures <- function(cfg, dataset, role_res, model_res = NULL) {
  save <- function(p, name, w = 8, h = 6) {
    ggsave(out_path(cfg, dataset, "figures", name), p, width = w, height = h, dpi = 200, bg = "white")
  }
  save(plot_k_diagnostics(role_res$diagnostics, dataset), "k_diagnostics.png", 7, 6)
  save(plot_strength_by_role(role_res$roles, dataset), "strength_by_role.png", 9, 5)
  save(plot_role_timeline(role_res$roles, dataset, "Africa"), "role_timeline_africa.png", 8, 11)
  if (!is.null(model_res)) save(plot_odds_ratios(model_res$binary, dataset), "odds_ratios.png", 10, 7)
}
