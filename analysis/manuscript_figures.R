# =============================================================================
# manuscript_figures.R — main-text and SI figures for the manuscript draft
#
#   Rscript analysis/manuscript_figures.R
#   (after run_pipeline.R, analysis/migrant_backbone.R and
#    analysis/refugee_role_changes.R)
#
# Role codes used in the manuscript (same position = same colour):
#   refugees (UNHCR):  R1 sender only    R2 two-way        R3 receiver only
#   migrants (backbone of net flows):
#                      M1 source         M2 intermediate   M3 terminal destination
#
# Writes to output/manuscript/:
#   Fig1_role_images.png      role-level "image" graphs + in/out separation
#   Fig2_role_maps.png        modal role per country, and how often it holds
#   Fig3_driver_odds.png      drivers of role membership (both networks)
#   Fig4_role_change_odds.png drivers of refugee role changes
#   Fig5_drought_by_income.png drought odds ratios by baseline-income tercile
#   FigS_*.png                supplementary figures
#   table1_role_summary.csv, image_flows_*.csv, manuscript_numbers.txt
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(readr); library(ggplot2)
  library(patchwork); library(ggrepel); library(igraph)
})
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
cfg_file <- commandArgs(trailingOnly = TRUE)   # optional: Rscript <script> my_config.yml
cfg <- read_config(if (length(cfg_file)) cfg_file[1] else "config.yml")
outdir <- dirname(out_path(cfg, "manuscript", "x"))
FONT <- "Liberation Sans"            # Arial-metric
# all figures: 300 dpi PNG, white background, widths in inches (7.2 = PNAS Nexus two-column width)
save_png <- function(p, name, w, h) ggsave(file.path(outdir, name), p, width = w, height = h,
                                           dpi = 300, bg = "white", device = ragg::agg_png)
theme_set(theme_minimal(9, base_family = FONT) +
            theme(plot.title.position = "plot", panel.grid.minor = element_blank(),
                  plot.tag = element_text(face = "bold", size = 11),
                  plot.title = element_text(face = "bold", size = 9.5),
                  legend.position = "bottom"))

# position colours: sending / intermediate / receiving (validated, CVD-safe)
POS_COL <- c(send = "#C8553D", mid = "#2D6FA8", recv = "#E3A33B")
REF_LAB <- c(`1` = "R1 sender only", `2` = "R2 two-way", `3` = "R3 receiver only")
MIG_LAB <- c(`1` = "M1 source", `2` = "M2 intermediate", `3` = "M3 terminal destination")
pos_of  <- c(`1` = "send", `2` = "mid", `3` = "recv")
ref_cols <- setNames(POS_COL[pos_of], REF_LAB)
mig_cols <- setNames(POS_COL[pos_of], MIG_LAB)
# Every number quoted in the manuscript is written to manuscript_numbers.txt
# by note(), so the text can be checked against the run that made the figures.
nums <- character()                  # numbers quoted in the text
note <- function(...) { s <- paste0(...); nums <<- c(nums, s); message(s) }

# ---- data -----------------------------------------------------------------------
ref_roles <- read_csv(out_path(cfg, "unhcr", "roles_by_country_period.csv"), show_col_types = FALSE) |>
  mutate(name = coalesce(country, iso3))
# safety check: role 1 must be the sender-only role and role 3 the receiver-only
# role, otherwise the R1-R3 codes would be attached to the wrong roles
lab_chk <- ref_roles |> distinct(role, label) |> arrange(role)
stopifnot(grepl("sender only", lab_chk$label[1]), grepl("receiver only", lab_chk$label[3]))
ref_roles <- ref_roles |> mutate(code = REF_LAB[as.character(role)])

mig_roles <- read_csv(out_path(cfg, "undesa_backbone", "backbone_roles.csv"), show_col_types = FALSE) |>
  mutate(name = coalesce(country, iso3))
# same check for the backbone roles (migrant_backbone.R orders them source, intermediate, terminal)
stopifnot(grepl("^source", mig_roles$label[mig_roles$role == 1][1]),
          grepl("^intermediate", mig_roles$label[mig_roles$role == 2][1]))
mig_roles <- mig_roles |> mutate(code = MIG_LAB[as.character(role)])

# flows used for the image graphs and networks: refugee arrivals (>= 1 person);
# net migrant flows; and the backbone (>= 10% of the origin's net outflow,
# same rule as migrant_backbone.R)
ref_flows <- read_flows(cfg, "unhcr") |> filter(flow >= 1)
c2 <- cfg; c2$flows$undesa$net <- TRUE
net_all <- read_flows(c2, "undesa")
net_flows <- net_all |> filter(flow >= 1)
bb_flows <- net_all |> group_by(origin, period) |> filter(flow >= 0.10 * sum(flow)) |> ungroup()

# ---- Fig 1 A/B: role-level image graphs -------------------------------------------
# Role-level "image" of a network: total people moving from each role to each
# role (both ends looked up in the same period), and the share of all people.
image_edges <- function(flows, roles) {
  rl <- roles |> select(iso3, period, role)
  flows |>
    inner_join(rl, by = c("origin" = "iso3", "period")) |>
    inner_join(rl, by = c("destination" = "iso3", "period"), suffix = c("_from", "_to")) |>
    group_by(from = role_from, to = role_to) |>
    summarise(people = sum(flow), links = n(), .groups = "drop") |>
    mutate(share = people / sum(people))
}
E_ref <- image_edges(ref_flows, ref_roles)
E_mig <- image_edges(bb_flows, mig_roles)
write_csv(E_ref, file.path(outdir, "image_flows_refugees.csv"))
write_csv(E_mig, file.path(outdir, "image_flows_migrants_backbone.csv"))

# Draw the image graph: one circle per role (area ~ number of country-periods),
# arrows between roles with width and label = share of all people moving,
# self-loops for flows within a role. Pure geometry; no results are computed here.
plot_image <- function(E, roles, labs, title, caption) {
  E <- E |> filter(share >= 0.001)   # hide links carrying < 0.1% of the flow
  P <- tibble(role = 1:3, x = c(0, 1, 2), y = c(0, 0.85, 0),
              n = as.integer(table(factor(roles$role, 1:3)))) |>
    mutate(r = 0.13 + 0.17 * sqrt(n / max(n)), col = POS_COL[pos_of[as.character(role)]],
           lab = sub(" .*", "", labs[as.character(role)]),
           desc = sub("^\\S+ ", "", labs[as.character(role)]))
  cx <- mean(P$x); cy <- mean(P$y)
  wid <- function(s) 0.4 + 6.5 * s
  # links between different roles
  L <- E |> filter(from != to) |> left_join(P, by = c("from" = "role")) |>
    left_join(P, by = c("to" = "role"), suffix = c("", "_t")) |>
    mutate(dx = x_t - x, dy = y_t - y, len = sqrt(dx^2 + dy^2), ux = dx / len, uy = dy / len,
           # shift each direction sideways so A->B and B->A do not overlap
           nx = uy, ny = -ux, off = 0.045,
           x0 = x + (r + 0.03) * ux + off * nx, y0 = y + (r + 0.03) * uy + off * ny,
           x1 = x_t - (r_t + 0.05) * ux + off * nx, y1 = y_t - (r_t + 0.05) * uy + off * ny,
           lx = (x0 + x1) / 2 + 0.085 * nx, ly = (y0 + y1) / 2 + 0.085 * ny,
           lw = wid(share))
  # self-loops, drawn on the side away from the triangle's centre
  S <- E |> filter(from == to) |> left_join(P, by = c("from" = "role")) |>
    mutate(ang = atan2(y - cy, x - cx), lr = 0.13, lw = wid(share))
  loop_path <- S |> rowwise() |> reframe({
    ccx <- x + (r + lr * 0.7) * cos(ang); ccy <- y + (r + lr * 0.7) * sin(ang)
    t <- seq(ang + pi + 1.1, ang + pi + 2 * pi - 1.1, length.out = 60)
    tibble(from = from, px = ccx + lr * cos(t), py = ccy + lr * sin(t), lw = lw,
           lx = ccx + (lr + 0.17) * cos(ang - pi / 2), ly = ccy + (lr + 0.17) * sin(ang - pi / 2), share = share)
  })
  lbl <- function(s) ifelse(s >= 0.005, sprintf("%.0f%%", 100 * s), "<1%")
  ggplot() +
    geom_segment(data = L, aes(x0, y0, xend = x1, yend = y1, linewidth = lw),
                 colour = "grey35", lineend = "butt", linejoin = "mitre",
                 arrow = arrow(length = unit(2.2, "mm"), type = "closed")) +
    geom_path(data = loop_path, aes(px, py, group = from, linewidth = lw), colour = "grey35",
              arrow = arrow(length = unit(2.2, "mm"), type = "closed")) +
    ggforce_circle(P) +
    geom_text(data = P, aes(x, y, label = lab), colour = "white", fontface = "bold",
              family = FONT, size = 3.6) +
    geom_text(data = P |> filter(role != 2), aes(x, y - r - 0.07, label = paste0(desc, "\n(n = ", scales::comma(n), ")")),
              family = FONT, size = 2.6, vjust = 1, lineheight = 0.9, colour = "grey20") +
    geom_text(data = P |> filter(role == 2), aes(x - r - 0.06, y, label = paste0(desc, "\n(n = ", scales::comma(n), ")")),
              family = FONT, size = 2.6, hjust = 1, lineheight = 0.9, colour = "grey20") +
    geom_label(data = L, aes(lx, ly, label = lbl(share)), family = FONT, size = 2.5,
               label.size = 0, label.padding = unit(0.6, "mm"), fill = "white", colour = "grey10") +
    geom_label(data = distinct(loop_path, from, lx, ly, share), aes(lx, ly, label = lbl(share)),
               family = FONT, size = 2.5, label.size = 0, label.padding = unit(0.6, "mm"),
               fill = "white", colour = "grey10") +
    scale_linewidth_identity() +
    coord_equal(xlim = c(-0.55, 2.55), ylim = c(-0.62, 1.62), expand = FALSE, clip = "off") +
    labs(title = title, subtitle = caption) +
    theme_void(9, base_family = FONT) +
    theme(plot.title = element_text(face = "bold", size = 9.5, hjust = 0.5),
          plot.subtitle = element_text(size = 7.5, colour = "grey30", hjust = 0.5),
          plot.tag = element_text(face = "bold", size = 11))
}
# filled circles without ggforce
ggforce_circle <- function(P) {
  C <- P |> rowwise() |> reframe(tibble(role = role, t = seq(0, 2 * pi, length.out = 90)) |>
                                   mutate(cx = x + r * cos(t), cy = y + r * sin(t), col = col))
  geom_polygon(data = C, aes(cx, cy, group = role, fill = I(col)), colour = "white", linewidth = 0.6)
}

pA <- plot_image(E_ref, ref_roles, REF_LAB, "Refugees (UNHCR arrivals, 1990–2020)",
                 "Arrow width and label: share of all refugees moving")
pB <- plot_image(E_mig, mig_roles, MIG_LAB, "Migrants (backbone of net flows, 1995–2020)",
                 "Arrow width and label: share of all net migration on the backbone")

# ---- Fig 1 C/D: separation scatter -------------------------------------------------
# label the two largest country-periods of each role in the last period
lab_pts <- function(d, xv, yv) d |> filter(period == max(period)) |> group_by(role) |>
  slice_max(.data[[xv]] + .data[[yv]], n = 2, with_ties = FALSE) |> ungroup()
pC_d <- ref_roles |> mutate(xo = log10(out_strength + 1), yi = log10(in_strength + 1))
pC <- ggplot(pC_d, aes(xo, yi, colour = code)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey75", linetype = 2) +
  geom_point(alpha = 0.55, size = 1, position = position_jitter(0.04, 0.04, seed = 1)) +
  geom_text_repel(data = lab_pts(pC_d, "xo", "yi"), aes(label = paste(name, period)), size = 2.3,
                  family = FONT, show.legend = FALSE, min.segment.length = 0, seed = 2, colour = "grey15") +
  scale_colour_manual(values = ref_cols, name = NULL) +
  theme(legend.text = element_text(size = 7)) +
  labs(x = "Refugees leaving, log10(people + 1)", y = "Refugees arriving, log10(people + 1)") +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 2.5, alpha = 1)))

# net in- and outflow of every country-period on the full net network
nm <- node_metrics(net_flows, 1) |> select(iso3, period, net_in = in_strength, net_out = out_strength)
pD_d <- mig_roles |> select(-any_of(c("net_in", "net_out"))) |> left_join(nm, by = c("iso3", "period")) |>
  mutate(xo = log10(coalesce(net_out, 0) + 1), yi = log10(coalesce(net_in, 0) + 1))
pD <- ggplot(pD_d, aes(xo, yi, colour = code)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey75", linetype = 2) +
  geom_point(alpha = 0.55, size = 1) +
  geom_text_repel(data = pD_d |> filter(period == 2020, name %in% c("United States", "Germany", "India",
                    "Pakistan", "Papua New Guinea", "South Africa", "China", "Russia")),
                  aes(label = paste(name, period)), size = 2.3, family = FONT, show.legend = FALSE,
                  min.segment.length = 0, seed = 2, colour = "grey15") +
  scale_colour_manual(values = mig_cols, name = NULL) +
  theme(legend.text = element_text(size = 7)) +
  labs(x = "Net migrant outflow, log10(people + 1)", y = "Net migrant inflow, log10(people + 1)") +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 2.5, alpha = 1)))

fig1 <- (pA | pB) / (pC | pD) + plot_layout(heights = c(1, 1.05)) +
  plot_annotation(tag_levels = "A")
save_png(fig1, "Fig1_role_images.png", 7.2, 7.4)

note("Refugee image shares: ", paste(sprintf("R%d->R%d %.1f%%", E_ref$from, E_ref$to, 100 * E_ref$share), collapse = "; "))
note("Migrant backbone image shares: ", paste(sprintf("M%d->M%d %.1f%%", E_mig$from, E_mig$to, 100 * E_mig$share), collapse = "; "))

# ---- Fig 2: maps of modal role and its stability -----------------------------------
world <- sf::read_sf("data_ref/ne_110m_admin_0_countries.geojson") |>   # Natural Earth 1:110m (public domain)
  mutate(iso3 = ifelse(ISO_A3 == "-99", ADM0_A3, ISO_A3)) |>
  mutate(iso3 = recode(iso3, SDS = "SSD", KOS = "XKX")) |>
  filter(ADM0_A3 != "ATA") |> select(iso3, geometry)
# per country: periods observed, most frequent (modal) role, share of periods
# in that role, and number of role switches
modal <- function(roles, labs) roles |> group_by(iso3) |>
  summarise(n = n(), mode = as.integer(names(which.max(table(role)))),
            share = max(table(role)) / n(), switches = sum(role != lag(role), na.rm = TRUE), .groups = "drop") |>
  mutate(code = labs[as.character(mode)],
         stab = cut(share, c(0, 0.6, 0.85, 1.01), labels = c("< 60%", "60–85%", "> 85%")))
# World map of the modal role; grey where no role holds >= 60% of periods or
# (refugees) where fewer than min_n years are recorded. Robinson projection.
map_panel <- function(m, cols, title, min_n = 1) {
  extra <- c(`No dominant role (< 60% of periods)` = "grey55", `Fewer than 5 years recorded` = "grey80")
  m <- m |> mutate(cat = case_when(n < min_n ~ "Fewer than 5 years recorded",
                                   share < 0.6 ~ "No dominant role (< 60% of periods)",
                                   TRUE ~ code))
  lv <- c(names(cols), names(extra)); if (min_n <= 1) lv <- setdiff(lv, "Fewer than 5 years recorded")
  d <- world |> left_join(m, by = "iso3") |> mutate(cat = factor(cat, lv)) |> sf::st_transform("+proj=robin")
  ggplot(d) +
    geom_sf(aes(fill = cat), colour = "white", linewidth = 0.1) +
    scale_fill_manual(values = c(cols, extra), na.value = "grey93", name = NULL, drop = FALSE,
                      na.translate = FALSE) +
    labs(title = title) +
    theme_void(9, base_family = FONT) +
    theme(legend.position = "right", legend.text = element_text(size = 7),
          legend.key.size = unit(3.5, "mm"), plot.title = element_text(face = "bold", size = 9.5),
          plot.tag = element_text(face = "bold", size = 11))
}
m_ref <- modal(ref_roles, REF_LAB); m_mig <- modal(mig_roles, MIG_LAB)
fig2 <- map_panel(m_ref, ref_cols, "Refugees: dominant role, 1990–2020", min_n = 5) /
  map_panel(m_mig, mig_cols, "Migrants (backbone): dominant role, 1995–2020") +
  plot_annotation(tag_levels = "A")
save_png(fig2, "Fig2_role_maps.png", 7.2, 5.6)
note("Refugee countries with < 5 years recorded: ", sum(m_ref$n < 5), " of ", nrow(m_ref), "; no dominant role (<60%): refugees ", sum(m_ref$n >= 5 & m_ref$share < 0.6), ", migrants ", sum(m_mig$share < 0.6), " of ", nrow(m_mig))
miss_iso <- setdiff(unique(c(m_ref$iso3, m_mig$iso3)), world$iso3)
note("Map: country codes without a 110m polygon (small states/territories, not drawn): ", length(miss_iso))

stab_ref <- role_stability(ref_roles); stab_mig <- role_stability(mig_roles)
note("Persistence refugees (consecutive years) = ", round(100 * stab_ref$overall_persistence, 1), "%")
note("Persistence migrants backbone (consecutive periods) = ", round(100 * stab_mig$overall_persistence, 1), "%")
note("Refugee transition matrix: ", paste(sprintf("%d->%d %.2f", stab_ref$transition_matrix$from,
     stab_ref$transition_matrix$to, stab_ref$transition_matrix$prob), collapse = "; "))
note("Backbone transition matrix: ", paste(sprintf("%d->%d %.2f", stab_mig$transition_matrix$from,
     stab_mig$transition_matrix$to, stab_mig$transition_matrix$prob), collapse = "; "))
note("Refugee countries always in same role: ", sum(m_ref$share == 1 & m_ref$n > 1), " of ", sum(m_ref$n > 1),
     "; migrants backbone: ", sum(m_mig$share == 1 & m_mig$n > 1), " of ", sum(m_mig$n > 1))

# ---- Fig 3: drivers of role membership ----------------------------------------------
# Display names for model terms (only terms listed here are plotted).
# The tercile terms come from the income_groups / slopes_by_group options in
# config.yml: <var>_x_g<k> is the slope of <var> among countries in baseline-
# income tercile k (1 = lowest); income_g2 / income_g3 are tercile main effects.
TERC <- c(g1 = "low income", g2 = "middle income", g3 = "high income")
DRV <- c(gdp_pc = "GDP per capita",
         income_g2 = "Middle income tercile (vs low)", income_g3 = "High income tercile (vs low)",
         agri_share_base = "Agriculture share of GDP (baseline)",
         conflict_events = "Conflict events", conflict_deaths = "Conflict deaths",
         drought12_ext = "Drought extent (SPEI-12)", drought12_ext_lag1 = "Drought extent, previous year",
         setNames(paste0("Drought extent: ", TERC), paste0("drought12_ext_x_", names(TERC))),
         setNames(paste0("Drought extent, prev. year: ", TERC), paste0("drought12_ext_lag1_x_", names(TERC))),
         spei48_dec = "SPEI-48 at period start",
         setNames(paste0("SPEI-48 at period start: ", TERC), paste0("spei48_dec_x_", names(TERC))),
         drought12_ext_x_gdp_pc = "Drought extent x GDP", drought12_ext_lag1_x_gdp_pc = "Drought (prev. year) x GDP",
         drought12_ext_x_agri_share_base = "Drought extent x agriculture share",
         drought12_ext_lag1_x_agri_share_base = "Drought (prev. year) x agriculture share",
         spei48_dec_x_agri_share_base = "SPEI-48 x agriculture share",
         spei_min = "SPEI minimum", spei_mean = "SPEI mean", spei_mean_cv = "SPEI variability (CV)",
         drought_affected = "Drought: affected", drought_deaths = "Drought: deaths",
         flood_affected = "Flood: affected", flood_deaths = "Flood: deaths",
         storm_affected = "Storm: affected", storm_deaths = "Storm: deaths",
         landslide_affected = "Landslide: affected", landslide_deaths = "Landslide: deaths")
# panel group of each term: economy (GDP, income terciles, agriculture share),
# conflict, climate (drought indicators and their interactions), disasters (EM-DAT)
GRP <- setNames(ifelse(grepl("^gdp|^income_g|^agri", names(DRV)), "Economy",
                ifelse(grepl("^conflict", names(DRV)), "Conflict",
                ifelse(grepl("spei|drought12", names(DRV)), "Climate", "Disasters"))), names(DRV))
# model results from run_pipeline.R (refugees) and migrant_backbone.R (migrants)
b_ref <- read_csv(out_path(cfg, "unhcr", "binary_role_models.csv"), show_col_types = FALSE) |>
  mutate(code = REF_LAB[as.character(role)], net = "Refugees")
b_mig <- read_csv(out_path(cfg, "undesa_backbone", "binary_role_models.csv"), show_col_types = FALSE) |>
  mutate(code = MIG_LAB[as.character(role)], net = "Migrants")
# Fig 3 shows every driver except the tercile-specific drought slopes, which
# have their own figure (Fig 5). Roles with < 20 members are not shown.
B <- bind_rows(b_ref, b_mig) |> filter(term %in% names(DRV), !grepl("_x_g[0-9]$", term), !unstable, n_events >= 20) |>
  mutate(driver = factor(DRV[term], rev(DRV)), grp = factor(GRP[term], c("Economy", "Conflict", "Climate", "Disasters")),
         panel = factor(paste0(code, "\n(", n_events, " of ", n_obs, ")"),
                        unique(paste0(code, "\n(", n_events, " of ", n_obs, ")"))),
         pos = pos_of[as.character(role)], sig = p.value < 0.05)
# filled, role-coloured points = p < 0.05; open grey points = p >= 0.05
fig3 <- ggplot(B, aes(odds_ratio, driver)) +
  geom_vline(xintercept = 1, colour = "grey60", linewidth = 0.3) +
  geom_errorbarh(aes(xmin = or_low, xmax = or_high, colour = ifelse(sig, pos, "ns")), height = 0, linewidth = 0.5) +
  geom_point(aes(colour = ifelse(sig, pos, "ns"), shape = sig), size = 1.6, fill = "white") +
  scale_colour_manual(values = c(POS_COL, ns = "grey60"), guide = "none") +
  scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 16), labels = c("p ≥ 0.05", "p < 0.05"), name = NULL) +
  scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4), labels = c("0.25", "0.5", "1", "2", "4")) +
  coord_cartesian(xlim = c(0.2, 5)) +
  facet_grid(grp ~ panel, scales = "free_y", space = "free_y") +
  labs(x = "Odds ratio per 1 SD (95% CI, log scale)", y = NULL) +
  theme(strip.text.x = element_text(size = 7.5, face = "bold"), strip.text.y = element_text(size = 7.5, angle = 0, hjust = 0),
        panel.spacing.x = unit(3, "mm"), axis.text.y = element_text(size = 7.5), axis.text.x = element_text(size = 7))
save_png(fig3, "Fig3_driver_odds.png", 7.2, 6.0)
note("Significant role drivers (p<0.05):\n", paste(capture.output(print(as.data.frame(
  B |> filter(sig) |> transmute(code, term, OR = round(odds_ratio, 2), lo = round(or_low, 2), hi = round(or_high, 2), p = signif(p.value, 2))))), collapse = "\n"))
note("Model sizes: ", paste(unique(sprintf("%s n_obs=%d events=%d", B$code, B$n_obs, B$n_events)), collapse = "; "))

# ---- Fig 4: refugee role changes ------------------------------------------------------
# results from refugee_role_changes.R: rows = specification (level / change), columns = type of change
T <- read_csv(out_path(cfg, "unhcr", "role_changes", "transition_models.csv"), show_col_types = FALSE) |>
  # tercile-specific drought slopes are shown in Fig 5, not here
  filter(!unstable, driver %in% names(DRV), !grepl("_x_g[0-9]$", driver)) |>
  mutate(driver = factor(DRV[driver], rev(DRV)), sig = p.value < 0.05,
         change = factor(paste0(change, "\n(", events, " of ", n, ")"), unique(paste0(change, "\n(", events, " of ", n, ")"))),
         spec = factor(ifelse(spec == "level in year t", "Driver level in year of change", "Change in driver from previous year"),
                       c("Driver level in year of change", "Change in driver from previous year")))
fig4 <- ggplot(T, aes(odds_ratio, driver)) +
  geom_vline(xintercept = 1, colour = "grey60", linewidth = 0.3) +
  geom_errorbarh(aes(xmin = or_low, xmax = or_high, colour = sig), height = 0, linewidth = 0.5) +
  geom_point(aes(colour = sig, shape = sig), size = 1.6, fill = "white") +
  scale_colour_manual(values = c(`FALSE` = "grey60", `TRUE` = "#1f3a5f"), guide = "none") +
  scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 16), labels = c("p ≥ 0.05", "p < 0.05"), name = NULL) +
  scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4), labels = c("0.25", "0.5", "1", "2", "4")) +
  coord_cartesian(xlim = c(0.2, 5)) + facet_grid(spec ~ change) +
  labs(x = "Odds ratio per 1 SD (95% CI, log scale)", y = NULL) +
  theme(strip.text = element_text(size = 7.5, face = "bold"), axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 7))
save_png(fig4, "Fig4_role_change_odds.png", 7.2, 7)
note("Role-change significant (p<0.05):\n", paste(capture.output(print(as.data.frame(T |> filter(sig) |>
  transmute(change = gsub("\n.*", "", change), spec, driver, OR = round(odds_ratio, 2), lo = round(or_low, 2),
            hi = round(or_high, 2), p = signif(p.value, 2), engine)))), collapse = "\n"))

# ---- Fig 5: drought slopes by baseline-income tercile ---------------------------------
# One odds ratio per drought indicator and income tercile, for every role model
# and every role-change model (level spec), with the Wald test of equal slopes
# across terciles printed in each panel. Shown only when the models were run
# with slopes_by_group (main config).
tw <- bind_rows(
  bind_rows(b_ref, b_mig) |> filter(n_events >= 20) |>
    transmute(model = paste0(code, "\n(", n_events, " of ", n_obs, ")"), term, odds_ratio, or_low, or_high,
              p.value, std.error, estimate, unstable),
  read_csv(out_path(cfg, "unhcr", "role_changes", "transition_models.csv"), show_col_types = FALSE) |>
    filter(spec == "level in year t") |>
    transmute(model = paste0(change, "\n(", events, " of ", n, ")"), term, odds_ratio, or_low, or_high,
              p.value, std.error, estimate, unstable))
if (any(grepl("_x_g[0-9]$", tw$term))) {
  model_order <- unique(tw$model)
  slopes <- tw |> filter(grepl("_x_g[0-9]$", term), !unstable) |>
    mutate(var = sub("_x_g[0-9]$", "", term), tercile = TERC[sub(".*_x_", "", term)],
           var_lab = c(drought12_ext = "Drought extent (year t / 5-yr mean)",
                       drought12_ext_lag1 = "Drought extent, previous year",
                       spei48_dec = "SPEI-48 at period start (higher = wetter)")[var],
           sig = p.value < 0.05, model = factor(model, model_order),
           tercile = factor(tercile, TERC))
  # Wald test of equal slopes, one per model and drought indicator (printed at the right of each row)
  wald <- tw |> filter(grepl("^wald_equal_slopes:", term)) |>
    mutate(var = sub("^wald_equal_slopes:", "", term), model = factor(model, model_order),
           var_lab = c(drought12_ext = "Drought extent (year t / 5-yr mean)",
                       drought12_ext_lag1 = "Drought extent, previous year",
                       spei48_dec = "SPEI-48 at period start (higher = wetter)")[var],
           txt = sprintf("equal slopes: p = %.2f", p.value))
  note("Equal-slope Wald tests:\n", paste(capture.output(print(as.data.frame(
    tw |> filter(grepl("^wald_", term)) |> transmute(model = gsub("\n.*", "", model), term,
                                                     chi2 = round(estimate, 2), df = std.error, p = signif(p.value, 2))))),
    collapse = "\n"))
  note("Tercile slopes:\n", paste(capture.output(print(as.data.frame(slopes |>
    transmute(model = gsub("\n.*", "", model), var, tercile, OR = round(odds_ratio, 2), lo = round(or_low, 2),
              hi = round(or_high, 2), p = signif(p.value, 2))))), collapse = "\n"))
  TERC_COL <- setNames(c("#C8553D", "#8C8C8C", "#2D6FA8"), TERC)   # low / middle / high income
  # one row per model, one panel per drought indicator; colour = income tercile
  slopes <- slopes |> mutate(model = factor(model, rev(model_order)))
  wald <- wald |> mutate(model = factor(model, rev(model_order)))
  fig5 <- ggplot(slopes, aes(odds_ratio, model, colour = tercile, group = tercile)) +
    geom_vline(xintercept = 1, colour = "grey60", linewidth = 0.3) +
    geom_errorbarh(aes(xmin = or_low, xmax = or_high), height = 0, linewidth = 0.5,
                   position = position_dodge(width = 0.65)) +
    geom_point(aes(shape = sig), size = 1.6, fill = "white", position = position_dodge(width = 0.65)) +
    geom_text(data = wald, aes(x = 9, y = model, label = txt), inherit.aes = FALSE, hjust = 1,
              size = 2.2, family = FONT, colour = "grey30") +
    scale_colour_manual(values = TERC_COL, name = "Baseline income tercile", breaks = TERC) +
    scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 16), labels = c("p ≥ 0.05", "p < 0.05"), name = NULL) +
    scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4), labels = c("0.25", "0.5", "1", "2", "4")) +
    coord_cartesian(xlim = c(0.25, 9)) +
    facet_grid(var_lab ~ ., scales = "free_y", space = "free_y", labeller = label_wrap_gen(18)) +
    labs(x = "Odds ratio per 1 SD of the drought indicator (95% CI, log scale)", y = NULL) +
    theme(strip.text.y = element_text(size = 7.5, face = "bold", angle = 0, hjust = 0),
          axis.text.y = element_text(size = 7.5), axis.text.x = element_text(size = 7),
          legend.text = element_text(size = 7), legend.box = "vertical",
          panel.grid.major.y = element_blank()) +
    guides(colour = guide_legend(order = 1), shape = guide_legend(order = 2))
  save_png(fig5, "Fig5_drought_by_income.png", 7.2, 6.2)
}

# ---- Table 1: role summary -------------------------------------------------------------
# Table 1: size, degrees, median strengths and share net receiving by role,
# plus the four largest countries of each role in the last period
summ <- function(roles, net) roles |> group_by(code) |>
  summarise(network = net, n = n(), countries = n_distinct(iso3),
            in_degree = mean(in_degree, na.rm = TRUE), out_degree = mean(out_degree, na.rm = TRUE),
            in_strength_median = median(in_strength, na.rm = TRUE), out_strength_median = median(out_strength, na.rm = TRUE),
            share_net_receiving = mean(net_flow > 0, na.rm = TRUE), .groups = "drop")
ref_bb <- ref_roles   # degrees on the refugee network
mig_bb <- mig_roles   # in/out degree here = backbone degree; strengths = all net flows
T1 <- bind_rows(summ(ref_bb, "Refugees"), summ(mig_bb, "Migrants (backbone)"))
examples <- function(roles) roles |> filter(period == max(period)) |> group_by(code) |>
  slice_max(in_strength + out_strength, n = 4, with_ties = FALSE) |>
  summarise(examples = paste(name, collapse = ", "), .groups = "drop")
T1 <- T1 |> left_join(bind_rows(examples(ref_roles), examples(mig_roles)), by = "code")
write_csv(T1, file.path(outdir, "table1_role_summary.csv"))
note("Table 1:\n", paste(capture.output(print(as.data.frame(T1 |> mutate(across(where(is.double), \(x) signif(x, 3)))))), collapse = "\n"))

# ---- SI: Africa timelines, side by side -------------------------------------------------
# role of every African country over time (tiles), refugees and migrants side by side
afr <- function(roles, cols, title, yrs) roles |> filter(area == "Africa") |>
  ggplot(aes(period, forcats::fct_rev(name), fill = code)) +
  geom_tile(colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = cols, name = NULL) + scale_x_continuous(breaks = yrs, expand = c(0, 0)) +
  labs(title = title, x = NULL, y = NULL) +
  theme(panel.grid = element_blank(), axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 6.5),
        legend.position = "bottom") + guides(fill = guide_legend(ncol = 1))
figS_afr <- afr(ref_roles, ref_cols, "Refugees (annual)", seq(1990, 2020, 5)) |
  afr(mig_roles, mig_cols, "Migrants, backbone (5-year periods)", seq(1995, 2020, 5))
save_png(figS_afr + plot_annotation(tag_levels = "A"), "FigS_africa_timelines.png", 7.2, 8.5)

# ---- SI: backbone country network (same style as the gross-flow network figure) ---------
add_shape("ellipse", clip = shape_noclip,
  plot = function(coords, v = NULL, params) {
    col <- params("vertex", "color"); size <- params("vertex", "size") / 200
    asp <- params("vertex", "size2")
    if (!is.null(v)) { col <- col[v]; if (length(size) > 1) size <- size[v]; asp <- asp[v] }
    n <- nrow(coords); col <- rep_len(col, n); size <- rep_len(size, n); asp <- rep_len(asp, n)
    t <- seq(0, 2 * pi, length.out = 49)
    for (i in seq_len(n)) polygon(coords[i, 1] + size[i] * sqrt(asp[i]) * cos(t),
                                  coords[i, 2] + size[i] / sqrt(asp[i]) * sin(t), col = col[i], border = NA)
  })
# Country network for one year: Fruchterman-Reingold layout weighted by log
# flow; node colour = role, size = volume, ellipse shape = in/out degree ratio
# (capped at 4), edge width = log flow; the 25 largest countries are labelled.
network_plot <- function(fy, ry, cols, file, title, foot, top_n = 25) {
  gy <- graph_from_data_frame(fy |> select(origin, destination, flow), directed = TRUE)
  vr <- ry[match(V(gy)$name, ry$iso3), ]
  vol <- coalesce(vr$vol, 0)
  asp <- pmin(pmax((vr$in_deg_all + 1e-9) / (vr$out_deg_all + 1e-9), 1 / 4), 4)
  lab <- rank(-vol, ties.method = "first") <= top_n
  ragg::agg_png(file.path(outdir, file), width = 2200, height = 2000, res = 300)
  par(mar = c(2.6, 0.5, 1.8, 0.5), family = FONT)
  set.seed(3)
  lay <- layout_with_fr(gy, weights = log10(E(gy)$flow + 1))
  lw <- log10(E(gy)$flow); lw <- 0.3 + 2.5 * (lw - min(lw)) / diff(range(lw))
  plot(gy, layout = lay, vertex.shape = "ellipse",
       vertex.size = 1.5 + 7 * sqrt(vol / max(vol)), vertex.size2 = asp,
       vertex.color = cols[as.character(vr$code)], vertex.frame.color = NA,
       vertex.label = ifelse(lab, vr$name, NA), vertex.label.cex = 0.45, vertex.label.family = FONT,
       vertex.label.color = "black", vertex.label.dist = 0.5,
       edge.width = lw, edge.arrow.size = 0.12, edge.color = adjustcolor("grey40", 0.4), edge.curved = 0.15)
  title(title, cex.main = 0.8, family = FONT)
  legend("bottomleft", legend = names(cols), fill = cols, border = NA, bty = "n", cex = 0.6)
  mtext(paste(strwrap(foot, 150), collapse = "\n"), side = 1, line = 1.2, cex = 0.45)
  dev.off()
}
yr <- 2020
fy <- bb_flows |> filter(period == yr)
deg_all <- node_metrics(net_flows |> filter(period == yr), 1) |>
  transmute(iso3, in_deg_all = in_degree, out_deg_all = out_degree, vol = in_strength + out_strength)
ry <- mig_roles |> filter(period == yr) |> select(iso3, name, code) |> left_join(deg_all, by = "iso3")
network_plot(fy, ry, mig_cols, "FigS_backbone_network_2020.png",
             "Migrants: backbone of main net corridors, 2015–2020",
             sprintf("Edges: net flows carrying ≥ 10%% of the origin's total net outflow (%d of %d net flows). Colour = role; size = total net in + out flow; shape = net in-degree / out-degree (wide: receives from more countries than it sends to). Labels: 25 largest.",
                     nrow(fy), sum(net_flows$period == yr)))
yr_r <- 2015
fyr <- ref_flows |> filter(period == yr_r)
deg_r <- node_metrics(fyr, 1) |> transmute(iso3, in_deg_all = in_degree, out_deg_all = out_degree, vol = in_strength + out_strength)
ryr <- ref_roles |> filter(period == yr_r) |> select(iso3, name, code) |> left_join(deg_r, by = "iso3")
network_plot(fyr, ryr, ref_cols, "FigS_refugee_network_2015.png",
             "Refugees: all flows recorded in 2015",
             sprintf("Edges: all %d refugee-arrival flows in the UNHCR file for 2015. Colour = role; size = total in + out flow; shape = in-degree / out-degree. Labels: 25 largest.", nrow(fyr)))
note("Backbone 2020: ", nrow(fy), " of ", sum(net_flows$period == yr), " net flows; refugees 2015: ", nrow(fyr), " flows")

# ---- SI: k diagnostics and MDS ---------------------------------------------------------
kd <- bind_rows(read_csv(out_path(cfg, "unhcr", "k_diagnostics.csv"), show_col_types = FALSE) |> mutate(net = "Refugees"))
sens <- read_csv(out_path(cfg, "undesa_backbone", "backbone_sensitivity.csv"), show_col_types = FALSE)
note("Backbone sensitivity:\n", paste(capture.output(print(as.data.frame(sens |> mutate(across(where(is.double), \(x) round(x, 3)))))), collapse = "\n"))
note("Refugee k diagnostics: ", paste(sprintf("k=%d sil=%.3f", kd$k, kd$avg_silhouette), collapse = "; "))

# headline counts
note("Refugee network: ", nrow(ref_flows), " flows; ", nrow(ref_roles), " country-years; ", n_distinct(ref_roles$iso3), " countries; role sizes ",
     paste(table(ref_roles$code), collapse = "/"))
note("Migrant backbone: ", nrow(bb_flows), " backbone ties of ", nrow(net_flows), " net flows; ", nrow(mig_roles), " country-periods; ",
     n_distinct(mig_roles$iso3), " countries; role sizes ", paste(table(mig_roles$code), collapse = "/"))
writeLines(nums, file.path(outdir, "manuscript_numbers.txt"))
msg("Done: ", outdir)
