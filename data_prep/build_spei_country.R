# =============================================================================
# build_spei_country.R — country-year drought indicators from gridded SPEIbase
#
#   Rscript data_prep/build_spei_country.R
#
# Inputs (put in data/raw/, names set below or in config.yml under
# drivers: spei_grid:):
#   spei12.nc, spei48.nc   SPEIbase v2.9 monthly grids, 0.5 deg
#                          (https://spei.csic.es/database.html)
#   countries              any polygon file with an ISO3 column (default:
#                          Natural Earth admin-0, 10m or 50m scale so small
#                          islands are included)
#   population (optional)  a population raster (e.g. GPWv4 or GHS-POP),
#                          used for population-weighted versions
#
# Why not the minimum over grid cells: SPEI is already standardized in every
# grid cell, so a country's area-weighted mean is comparable across countries.
# The minimum over cells is not: a country with more cells has a lower expected
# minimum, so that index measures country size (rho = -0.8 to -0.9 with log area
# for the earlier SPEI files). Every indicator below is a weighted mean or a
# weighted share and does not depend on the number of cells.
#
# Output: data/raw/spei_country_year.csv with one row per iso3 x year:
#   spei12_dec        area-weighted mean SPEI-12 in December
#                     (water balance of the calendar year relative to normal)
#   drought12_ext     share of the country's area in drought (SPEI-12 <= -1),
#                     averaged over the 12 months of the year (0..1)
#   drought12_severe  as above, SPEI-12 <= -1.5
#   spei48_dec        area-weighted mean SPEI-48 in December
#                     (four-year cumulative water balance)
#   *_pop versions    same, weighted by population instead of area (if given)
#
# Formulas. For country c, month m, grid cells g with SPEI value S_gm:
#   weight     w_gc = f_gc * A_g   (f_gc = fraction of cell g inside c,
#                                   A_g = cell area in km2; x population for _pop)
#   mean       Sbar_cm = sum_g w_gc S_gm / sum_g w_gc
#   extent     X_cm    = sum_g w_gc 1[S_gm <= -1] / sum_g w_gc     (0..1)
#   annual     drought12_ext_cy = (1/12) sum_{m in year y} X_cm   (SPEI-12)
#              spei12_dec_cy    = Sbar_c,Dec(y)                    (SPEI-12)
#              spei48_dec_cy    = Sbar_c,Dec(y)                    (SPEI-48)
# Cells with missing SPEI (e.g. ice, deserts with no precipitation) are left
# out of both numerator and denominator in that month.
# =============================================================================

suppressPackageStartupMessages({ library(terra); library(sf); library(dplyr); library(tidyr); library(readr); library(exactextractr) })

# Options (defaults below); any can be overridden on the command line as
# name=value, e.g. Rscript data_prep/build_spei_country.R years=1980:2020
args <- commandArgs(trailingOnly = TRUE)
opt <- list(
  spei12     = "data/raw/spei12.nc",
  spei48     = "data/raw/spei48.nc",
  countries  = "data/raw/ne_10m_admin_0_countries/ne_10m_admin_0_countries.shp",
  iso_col    = "ISO_A3",
  population = NA,                       # e.g. "data/raw/gpw_v4_population_count_2015_30_min.tif"
  years      = 1990:2020,
  threshold  = -1, severe = -1.5,       # SPEI cut-offs for drought and severe drought
  out        = "data/raw/spei_country_year.csv")
for (a in args) {                      # e.g. years=1990:2020 population=data/raw/pop.tif
  kv <- strsplit(a, "=")[[1]]
  # "1990:2020" is evaluated to a vector of years; other values are converted
  # to numbers / logicals where possible
  opt[[kv[1]]] <- if (grepl("^\\d{4}:\\d{4}$", kv[2])) eval(parse(text = kv[2])) else type.convert(kv[2], as.is = TRUE)
}

# Read one SPEIbase NetCDF file as a multi-layer raster (one layer per month)
# and keep the layers of the requested years, named "YYYY-MM".
read_months <- function(f, years) {
  r <- rast(f)
  tm <- time(r)
  if (is.null(tm) || all(is.na(tm))) stop("No time axis in ", f)
  keep <- which(as.integer(format(tm, "%Y")) %in% years)
  r <- r[[keep]]; names(r) <- format(tm[keep], "%Y-%m"); r
}

sf_use_s2(FALSE)                       # planar union; Natural Earth has some invalid rings
poly <- st_make_valid(st_read(opt$countries, quiet = TRUE))
iso <- poly[[opt$iso_col]]
# Natural Earth codes some countries "-99" in ISO_A3; fall back to ADM0_A3
if ("ADM0_A3" %in% names(poly)) iso <- ifelse(is.na(iso) | iso == "-99", poly$ADM0_A3, iso)
poly$iso3 <- iso
# dissolve all polygons with the same ISO3 into one (multi)polygon per country
poly <- poly |> group_by(iso3) |> summarise(.groups = "drop")

# Weights: exact share of each cell inside the country (exactextractr) x cell
# area (km2), optionally x population. Cells that only touch a small island
# still count, which is what fixes the island coverage problem. The weights are
# computed once per country and reused for every month.
r12 <- read_months(opt$spei12, opt$years)
r48 <- read_months(opt$spei48, opt$years)
# area of every 0.5-degree cell in km2 (shrinks towards the poles)
area <- values(cellSize(r12[[1]], unit = "km"), mat = FALSE)
# optional population per SPEI cell (population raster summed onto the SPEI grid)
pop <- if (!is.na(opt$population)) values(resample(rast(opt$population), r12[[1]], method = "sum"), mat = FALSE) else NULL
# for every country: the cell numbers it overlaps and the fraction of each
# cell inside it (coverage_fraction). Both SPEI grids share the same cells.
cov <- exactextractr::exact_extract(r12[[1]], poly, include_cell = TRUE, progress = FALSE)
names(cov) <- poly$iso3

# All annual indicators for one weighting.
#   wvec    weight per grid cell (area, or population)
#   suffix  "" for area weights, "_pop" for population weights
annual <- function(wvec, suffix) {
  one_index <- function(r, label) {
    V <- values(r)                                   # cells x months
    ym <- names(r)
    res <- lapply(names(cov), function(iso) {
      cc <- cov[[iso]]
      w <- cc$coverage_fraction * wvec[cc$cell]      # w_gc for the country's cells
      X <- V[cc$cell, , drop = FALSE]                 # SPEI of those cells, all months
      # cell-months with missing SPEI or zero weight get weight 0 (left out)
      ok <- !is.na(X) & !is.na(w) & w > 0
      W <- ifelse(ok, w, 0); X[!ok] <- 0
      sw <- colSums(W)                                 # sum of weights per month
      tibble(iso3 = iso, ym = ym,
             mean = ifelse(sw > 0, colSums(W * X) / sw, NA_real_),
             ext  = ifelse(sw > 0, colSums(W * (X <= opt$threshold & ok)) / sw, NA_real_),
             sev  = ifelse(sw > 0, colSums(W * (X <= opt$severe & ok)) / sw, NA_real_))
    })
    bind_rows(res) |> mutate(year = as.integer(substr(ym, 1, 4)), month = as.integer(substr(ym, 6, 7)))
  }
  m12 <- one_index(r12, "spei12"); m48 <- one_index(r48, "spei48")
  # monthly -> annual: December value of the mean; 12-month average of the extents
  a12 <- m12 |> group_by(iso3, year) |>
    summarise(spei12_dec = mean[month == 12][1], drought12_ext = mean(ext),
              drought12_severe = mean(sev), .groups = "drop")
  a48 <- m48 |> filter(month == 12) |> transmute(iso3, year, spei48_dec = mean)
  out <- full_join(a12, a48, by = c("iso3", "year")) |>
    filter(!(is.na(spei12_dec) & is.na(drought12_ext) & is.na(spei48_dec)))
  if (suffix != "") names(out)[-(1:2)] <- paste0(names(out)[-(1:2)], suffix)
  out
}
res <- annual(area, "")               # area-weighted (used in the analysis)
if (!is.null(pop)) res <- full_join(res, annual(pop, "_pop"), by = c("iso3", "year"))
write_csv(res |> arrange(iso3, year), opt$out)
message("Wrote ", opt$out, ": ", n_distinct(res$iso3), " countries x ", n_distinct(res$year), " years")
