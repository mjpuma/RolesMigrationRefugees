# =============================================================================
# drivers.R — covariates as one tidy country-year panel (iso3, year, ...)
#
# Each loader reads one raw source and returns long data keyed on ISO3 + year.
# build_driver_panel() joins them; aggregate_drivers() collapses the annual
# panel into the periods of a given flow dataset using the rules in config.yml.
#
# Sources and the variables they produce (names as used in config.yml):
#   load_population   UN WPP           pop (people, mid-year)
#   load_gdp          World Bank WDI   gdp_pc (current US$ per person)
#   load_spei         legacy SPEI CSVs spei01_min etc. (only in the
#                                      config_original_spei01 / Jordan configs)
#   load_conflict     UCDP-GED         conflict_events, conflict_deaths
#   load_disasters    EM-DAT           <type>_deaths, <type>_affected
#   load_climate      our SPEI build   drought12_ext, spei48_dec, ... (from
#                                      data_prep/build_spei_country.R)
#   load_long_files   any iso3/year    e.g. agri_share (World Bank NV.AGR.TOTL.ZS)
#
# Missing sources are skipped: a loader whose config entry is absent returns
# NULL, and compact() in build_driver_panel() drops it.
# =============================================================================

# UN World Population Prospects: total population on 1 July, given in
# thousands in the file, converted to persons. Used (a) to turn deaths and
# affected counts into per-capita rates and (b) as a descriptive variable.
load_population <- function(cfg) {
  read_csv(in_root(cfg, cfg$drivers$population$file), show_col_types = FALSE,
           col_select = c("ISO3_code", "LocTypeName", "Time", "TPopulation1July")) |>
    # keep countries/areas only (drop regions, income groups, the world total)
    filter(LocTypeName == "Country/Area", !is.na(ISO3_code)) |>
    transmute(iso3 = ISO3_code, year = as.integer(Time), pop = TPopulation1July * 1000)
}

# World Bank GDP per capita (current US$, series NY.GDP.PCAP.CD).
# The WDI export is wide (one column per year, named "1990 [YR1990]"); ".."
# marks missing values.
load_gdp <- function(cfg, xw = NULL) {
  gc <- cfg$drivers$gdp
  gdp <- read_csv(in_root(cfg, gc$file), show_col_types = FALSE,
                  na = c("..", "", "NA")) |>
    filter(!is.na(`Country Code`), `Series Code` == "NY.GDP.PCAP.CD") |>
    select(iso3 = `Country Code`, matches("^\\d{4} \\[YR")) |>
    pivot_longer(-iso3, names_to = "year", values_to = "gdp_pc") |>
    # "1990 [YR1990]" -> 1990
    mutate(year = as.integer(substr(year, 1, 4)), gdp_pc = as.numeric(gdp_pc))
  # Optional supplement for territories missing from the World Bank file
  # (wide: first column = country name, one column per year). It only fills
  # country-years that the World Bank file leaves empty.
  # In the main config this is data/raw/new_gdp_percapita_rows.csv: the rows of
  # Jordan's new_gdp.csv that are GDP per capita in US$ (the French overseas
  # departments, which are total GDP in million EUR, are excluded).
  if (!is.null(gc$supplement) && !is.null(xw)) {
    sup <- read_csv(in_root(cfg, gc$supplement), show_col_types = FALSE,
                    locale = locale(encoding = gc$supplement_encoding %||% "UTF-8")) |>
      rename(name = 1) |>
      pivot_longer(-name, names_to = "year", values_to = "gdp_sup") |>
      mutate(iso3 = name_to_iso3(name, xw, "GDP supplement"), year = as.integer(year),
             gdp_sup = as.numeric(gdp_sup)) |>
      filter(!is.na(iso3), !is.na(gdp_sup)) |>
      select(iso3, year, gdp_sup)
    gdp <- full_join(gdp, sup, by = c("iso3", "year"))
    # count how many World Bank gaps the supplement fills (reported in the log)
    n_fill <- sum(is.na(gdp$gdp_pc) & !is.na(gdp$gdp_sup))
    # coalesce(): World Bank value if present, otherwise the supplement
    gdp <- gdp |> mutate(gdp_pc = coalesce(gdp_pc, gdp_sup)) |> select(-gdp_sup)
    msg("  GDP supplement filled ", n_fill, " country-years")
  }
  gdp
}

# Legacy SPEI files (one CSV per index, rows = Natural Earth country names,
# columns X1990 ... X2020). Kept so the original SPEI-01 minimum and Jordan's
# SPEI-48 average-minimum can be rerun for comparison (Table S1). The main
# analysis uses load_climate() instead.
#   cfg$drivers$spei$files  named list: variable name -> file
load_spei <- function(cfg, xw) {
  rules <- cfg$drivers$aggregation
  imap(cfg$drivers$spei$files, function(f, var) {
    # rule used to combine several polygons with the same ISO3 (default mean)
    fun <- match.fun(rules[[var]] %||% "mean")
    read_csv(in_root(cfg, f), show_col_types = FALSE, col_select = -1) |>   # drop the row-number column
      select(name = countriesADMIN, starts_with("X")) |>
      pivot_longer(-name, names_to = "year", values_to = "value") |>
      mutate(iso3 = name_to_iso3(name, xw, "SPEI"),
             year = as.integer(sub("X", "", year))) |>
      filter(!is.na(iso3)) |>
      # several polygons can map to one ISO3 (e.g. Gaza + West Bank)
      group_by(iso3, year) |>
      summarise(!!var := fun(value, na.rm = TRUE), .groups = "drop") |>
      # optional: standardize within each country (removes level differences,
      # e.g. minimum-over-grid-cell indices that scale with country area)
      group_by(iso3) |>
      mutate(!!var := if (var %in% unlist(cfg$drivers$spei$within_country))
        (.data[[var]] - mean(.data[[var]], na.rm = TRUE)) / sd(.data[[var]], na.rm = TRUE)
        else .data[[var]]) |>
      ungroup()
  }) |>
    # one column per SPEI file, joined by country and year
    reduce(full_join, by = c("iso3", "year"))
}

# UCDP Georeferenced Event Dataset: one row per event (a day of organized
# violence at a location with at least one death). Summed to country-years.
#   conflict_events  number of GED events
#   conflict_deaths  sum of the "best" fatality estimate
load_conflict <- function(cfg, xw) {
  read_csv(in_root(cfg, cfg$drivers$conflict$file), show_col_types = FALSE,
           col_select = c("country", "year", "best")) |>
    mutate(iso3 = name_to_iso3(country, xw, "UCDP-GED")) |>
    filter(!is.na(iso3)) |>
    group_by(iso3, year = as.integer(year)) |>
    summarise(conflict_events = n(),          # number of GED events (not conflicts)
              conflict_deaths = sum(best, na.rm = TRUE),
              .groups = "drop")
}

# EM-DAT disasters: one row per disaster event. Kept types are listed in
# cfg$drivers$disasters$types (drought, flood, storm, landslide in the main
# config). Output: one column per type and measure, e.g. flood_deaths,
# flood_affected, summed over the events of a country-year. Country-years with
# no event are 0 (values_fill = 0 here, and zero_fill in aggregate_drivers()).
load_disasters <- function(cfg) {
  dc <- cfg$drivers$disasters
  f <- in_root(cfg, dc$file)
  csv <- sub("\\.xlsx?$", ".csv", f)
  raw <- if (file.exists(csv)) {
    # a CSV export of the same sheet (header row first) avoids needing readxl
    read_csv(csv, show_col_types = FALSE, guess_max = 1e5)
  } else {
    readxl::read_xlsx(f, skip = dc$skip_rows %||% 0, .name_repair = "minimal")
  }
  raw |>
    mutate(`Disaster Type` = str_trim(`Disaster Type`)) |>   # EM-DAT has "Extreme temperature " with a trailing space
    filter(`Disaster Type` %in% dc$types) |>
    transmute(iso3 = ISO, year = as.integer(Year),
              type = str_replace_all(tolower(`Disaster Type`), " ", "_"),   # "Flood" -> "flood"
              deaths = as.numeric(`Total Deaths`),
              affected = as.numeric(`Total Affected`)) |>
    group_by(iso3, year, type) |>
    summarise(across(c(deaths, affected), \(x) sum(x, na.rm = TRUE)), .groups = "drop") |>
    # long (one row per type) -> wide (flood_deaths, flood_affected, ...)
    pivot_wider(names_from = type, values_from = c(deaths, affected),
                names_glue = "{type}_{.value}", values_fill = 0)
}

# Country-year climate indicators already keyed on ISO3 (long format: iso3, year, ...),
# e.g. data/raw/spei_country_year.csv from data_prep/build_spei_country.R.
# That file holds the area-weighted drought indicators used in the main
# analysis: drought12_ext (share of land area with SPEI-12 <= -1, averaged
# over the 12 months), spei48_dec (area-weighted mean December SPEI-48), and
# the other columns documented in data_prep/SPEI_README.md.
load_climate <- function(cfg) {
  f <- cfg$drivers$climate$file
  if (is.null(f)) return(NULL)
  read_csv(in_root(cfg, f), show_col_types = FALSE) |> mutate(year = as.integer(year))
}

# Any further long-format country-year files (iso3, year, <variables>), e.g.
# data/raw/wdi_agri_share.csv. Each file's variables join the panel by iso3 + year.
# distinct() guards against a file that repeats a country-year (keeps the first).
load_long_files <- function(cfg) {
  fs <- unlist(cfg$drivers$long_files)
  if (!length(fs)) return(NULL)
  map(fs, \(f) read_csv(in_root(cfg, f), show_col_types = FALSE) |>
        mutate(year = as.integer(year)) |> distinct(iso3, year, .keep_all = TRUE)) |>
    reduce(full_join, by = c("iso3", "year"))
}

# Annual country panel with every driver.
# full_join keeps every country-year that appears in ANY source; missing
# values are handled later (zero-fill for event counts, otherwise the row
# drops out of the models). add_successor_rows() then copies modern-state
# rows onto historical flow codes (SUD, SCG; see utils.R).
build_driver_panel <- function(cfg) {
  xw <- load_crosswalk(cfg)
  msg("Loading driver data ...")
  panel <- list(load_population(cfg), load_gdp(cfg, xw), load_spei(cfg, xw),
                load_conflict(cfg, xw), load_disasters(cfg), load_climate(cfg),
                load_long_files(cfg)) |>
    compact() |>                                   # drop sources that returned NULL
    reduce(full_join, by = c("iso3", "year")) |>
    add_successor_rows(cfg)
  msg("  driver panel: ", n_distinct(panel$iso3), " countries x ",
      n_distinct(panel$year), " years, ", ncol(panel) - 2, " variables")
  panel
}

# Collapse annual drivers into the periods used by one flow dataset.
#   panel  annual panel from build_driver_panel()
#   flows  output of read_flows(); supplies each period's window
# For each country and period, every driver is summarised over the years
# window_start..window_end by its rule in cfg$drivers$aggregation:
#   mean  average over the window (default; e.g. 5-year mean drought extent)
#   sum   total over the window
#   max / min
#   first value in the first year with data (main config: spei48_dec, i.e.
#         December SPEI-48 of the window's first year, which summarises the
#         water balance of the 48 months up to that December: the conditions
#         at the start of the period)
#   last  value in the last year with data
# For annual refugee data the window is one year, so every rule returns that
# year's value. Per-capita variables are divided by the window's mean
# population AFTER aggregation.
aggregate_drivers <- function(cfg, panel, flows) {
  dcfg <- cfg$drivers
  windows <- flows |> distinct(period, window_start, window_end)
  vars <- setdiff(names(panel), c("iso3", "year", "pop"))

  # zero-fill: no record in UCDP / EM-DAT = no event, but ONLY within the years
  # the source covers. Beyond its last year (e.g. 2021+ when the file ends in
  # 2020) values stay missing, so they are not mistaken for "no events".
  zf <- vars[matches_any(vars, unlist(dcfg$zero_fill))]
  for (v in zf) {
    yrs <- panel$year[!is.na(panel[[v]])]
    if (!length(yrs)) next
    inside <- panel$year >= min(yrs) & panel$year <= max(yrs)
    panel[[v]][inside & is.na(panel[[v]])] <- 0
  }
  # report driver coverage so gaps with the flow years are visible
  cover <- vapply(vars, function(v) {
    y <- panel$year[!is.na(panel[[v]])]
    if (length(y)) paste0(min(y), "-", max(y)) else "none"
  }, character(1))
  flow_years <- range(windows$window_start, windows$window_end)
  late <- names(cover)[vapply(vars, function(v) {
    y <- panel$year[!is.na(panel[[v]])]; length(y) && max(y) < flow_years[2] }, logical(1))]
  if (length(late)) msg("  note: flows run to ", flow_years[2], " but these drivers end earlier: ",
                        paste(sprintf("%s (%s)", late, cover[late]), collapse = ", "),
                        " -> later periods have missing drivers and drop out of the models")

  # Expand each period into its calendar years (period 1995 -> 1991..1995)
  # and attach every country's driver values for those years. Sorting by year
  # within country-period matters for the "first" and "last" rules.
  long <- windows |>
    rowwise() |>
    mutate(year = list(seq(window_start, window_end))) |>
    unnest(year) |>
    inner_join(panel, by = "year", relationship = "many-to-many") |>
    arrange(iso3, period, year)

  # Apply one aggregation rule to the values of one country-period.
  # All-missing input gives NA (not 0, -Inf or NaN).
  agg_one <- function(x, rule) {
    if (all(is.na(x))) return(NA_real_)
    switch(rule, mean = mean(x, na.rm = TRUE), max = max(x, na.rm = TRUE),
           min = min(x, na.rm = TRUE), sum = sum(x, na.rm = TRUE),
           last = dplyr::last(x[!is.na(x)]),     # value in the final year of the period
           first = dplyr::first(x[!is.na(x)]))
  }
  # rule for a variable: its own entry in the config, else the default rule
  rule_for <- function(v) dcfg$aggregation[[v]] %||% dcfg$aggregation$default

  out <- long |>
    group_by(iso3, period) |>
    summarise(pop = mean(pop, na.rm = TRUE),
              across(all_of(vars), \(x) agg_one(x, rule_for(cur_column()))),
              .groups = "drop") |>
    mutate(pop = ifelse(is.nan(pop), NA, pop))     # mean of nothing is NaN -> NA

  # per-capita conversion for variables matching cfg$drivers$per_capita
  # (deaths and affected counts): count / mean population of the window
  pc <- vars[matches_any(vars, unlist(dcfg$per_capita))]
  out |> mutate(across(all_of(pc), \(x) x / pop))
}
