# =============================================================================
# spei_size_check.R — do country drought indices depend on country size?
#
#   Rscript analysis/spei_size_check.R [config.yml]
#
# Compares the long-run (1990-2020) country mean of three drought indices with
# country land area: the original SPEI-01 annual minimum, Jordan's SPEI-48
# annual average of monthly minima (fromJordan/), and the area-weighted SPEI-12
# drought extent from data_prep/build_spei_country.R. An index built from the
# minimum over grid cells scales with the number of cells, i.e. with size.
# Writes output/manuscript/FigS_spei_size_check.png (in the config's output_dir).
# =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr); library(ggplot2); library(sf) })
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
cfg_file <- commandArgs(trailingOnly = TRUE)
cfg <- read_config(if (length(cfg_file)) cfg_file[1] else "config.yml")
outdir <- dirname(out_path(cfg, "manuscript", "x"))
xw <- load_crosswalk(cfg)

# country land area (km2) from the Natural Earth 10m polygons, in an
# equal-area projection (Mollweide, ESRI:54009)
sf_use_s2(FALSE)
shp <- in_root(cfg, "data/raw/ne_10m_admin_0_countries/ne_10m_admin_0_countries.shp")
w <- st_read(shp, quiet = TRUE)
w$iso3 <- ifelse(w$ISO_A3 == "-99", w$ADM0_A3, w$ISO_A3)
area <- w |> st_transform("ESRI:54009") |> mutate(km2 = as.numeric(st_area(geometry)) / 1e6) |>
  st_drop_geometry() |> group_by(iso3) |> summarise(km2 = sum(km2))

# long-run (all years) country mean of a legacy wide SPEI file
# (rows = Natural Earth names, columns X1990 ... X2020)
wide_mean <- function(f, enc = "UTF-8") {
  read_csv(in_root(cfg, f), show_col_types = FALSE, locale = locale(encoding = enc), col_select = -1) |>
    mutate(iso3 = name_to_iso3(countriesADMIN, xw, f)) |> filter(!is.na(iso3)) |>
    pivot_longer(starts_with("X")) |> group_by(iso3) |> summarise(value = mean(value, na.rm = TRUE))
}
d <- bind_rows(
  wide_mean("data/raw/spei01-annual-min.csv", "latin1") |> mutate(index = "A. SPEI-01 annual minimum (original)"),
  wide_mean("fromJordan/spei48-annual-avg-min-expanded.csv") |> mutate(index = "B. SPEI-48 mean of monthly minima (6 Nov update)"),
  read_csv(in_root(cfg, "data/raw/spei_country_year.csv"), show_col_types = FALSE) |>
    group_by(iso3) |> summarise(value = mean(drought12_ext, na.rm = TRUE)) |>
    mutate(index = "C. SPEI-12 drought extent, area-weighted (this study)")) |>
  inner_join(area, by = "iso3")
# Spearman rank correlation of each index with log10 area. A value near 0
# means the index does not depend on country size.
rho <- d |> group_by(index) |>
  summarise(rho = cor(value, log10(km2), method = "spearman", use = "complete.obs"), n = n(), .groups = "drop")
print(rho)
p <- ggplot(d, aes(log10(km2), value)) +
  geom_point(alpha = 0.5, size = 1, colour = "#2D6FA8") +
  geom_smooth(method = "loess", se = FALSE, colour = "grey30", linewidth = 0.5, formula = y ~ x) +
  geom_text(data = rho, aes(x = -Inf, y = Inf, label = sprintf("Spearman rho = %.2f (n = %d)", rho, n)),
            hjust = -0.05, vjust = 1.3, size = 2.8, family = "Liberation Sans") +
  facet_wrap(~index, scales = "free_y", ncol = 3, labeller = label_wrap_gen(30)) +
  labs(x = "Country land area, log10(km²)", y = "Country mean, 1990–2020") +
  theme_minimal(9, base_family = "Liberation Sans") +
  theme(strip.text = element_text(face = "bold", size = 8, hjust = 0), panel.grid.minor = element_blank())
ggsave(file.path(outdir, "FigS_spei_size_check.png"), p, width = 7.2, height = 2.9, dpi = 300, bg = "white",
       device = ragg::agg_png)
write_csv(rho, file.path(outdir, "spei_size_check.csv"))
