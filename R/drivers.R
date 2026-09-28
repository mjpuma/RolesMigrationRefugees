# =============================================================================
# drivers.R — covariates as one tidy country-year panel (iso3, year, ...)
#
# Each loader reads one raw source and returns long data keyed on ISO3 + year.
# build_driver_panel() joins them; aggregate_drivers() collapses the annual
# panel into the periods of a given flow dataset using the rules in config.yml.
# =============================================================================

load_population <- function(cfg) {
  read_csv(in_root(cfg, cfg$drivers$population$file), show_col_types = FALSE,
           col_select = c("ISO3_code", "LocTypeName", "Time", "TPopulation1July")) |>
    filter(LocTypeName == "Country/Area", !is.na(ISO3_code)) |>
    transmute(iso3 = ISO3_code, year = as.integer(Time), pop = TPopulation1July * 1000)
}

load_gdp <- function(cfg) {
  read_csv(in_root(cfg, cfg$drivers$gdp$file), show_col_types = FALSE,
           na = c("..", "", "NA")) |>
    filter(!is.na(`Country Code`), `Series Code` == "NY.GDP.PCAP.CD") |>
    select(iso3 = `Country Code`, matches("^\\d{4} \\[YR")) |>
    pivot_longer(-iso3, names_to = "year", values_to = "gdp_pc") |>
    mutate(year = as.integer(substr(year, 1, 4)), gdp_pc = as.numeric(gdp_pc))
}

load_spei <- function(cfg, xw) {
  rules <- cfg$drivers$aggregation
  imap(cfg$drivers$spei$files, function(f, var) {
    fun <- match.fun(rules[[var]] %||% "mean")
    read_csv(in_root(cfg, f), show_col_types = FALSE, col_select = -1) |>
      select(name = countriesADMIN, starts_with("X")) |>
      pivot_longer(-name, names_to = "year", values_to = "value") |>
      mutate(iso3 = name_to_iso3(name, xw, "SPEI"),
             year = as.integer(sub("X", "", year))) |>
      filter(!is.na(iso3)) |>
      # several polygons can map to one ISO3 (e.g. Gaza + West Bank)
      group_by(iso3, year) |>
      summarise(!!var := fun(value, na.rm = TRUE), .groups = "drop")
  }) |>
    reduce(full_join, by = c("iso3", "year"))
}

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

load_disasters <- function(cfg) {
  dc <- cfg$drivers$disasters
  readxl::read_xlsx(in_root(cfg, dc$file), skip = dc$skip_rows %||% 0,
                    .name_repair = "minimal") |>
    filter(`Disaster Type` %in% dc$types) |>
    transmute(iso3 = ISO, year = as.integer(Year),
              type = str_replace_all(tolower(`Disaster Type`), " ", "_"),
              deaths = as.numeric(`Total Deaths`),
              affected = as.numeric(`Total Affected`)) |>
    group_by(iso3, year, type) |>
    summarise(across(c(deaths, affected), \(x) sum(x, na.rm = TRUE)), .groups = "drop") |>
    pivot_wider(names_from = type, values_from = c(deaths, affected),
                names_glue = "{type}_{.value}", values_fill = 0)
}

# Annual country panel with every driver
build_driver_panel <- function(cfg) {
  xw <- load_crosswalk(cfg)
  msg("Loading driver data ...")
  panel <- list(load_population(cfg), load_gdp(cfg), load_spei(cfg, xw),
                load_conflict(cfg, xw), load_disasters(cfg)) |>
    reduce(full_join, by = c("iso3", "year")) |>
    add_successor_rows(cfg)
  msg("  driver panel: ", n_distinct(panel$iso3), " countries x ",
      n_distinct(panel$year), " years, ", ncol(panel) - 2, " variables")
  panel
}

# Collapse annual drivers into the periods used by one flow dataset
aggregate_drivers <- function(cfg, panel, flows) {
  dcfg <- cfg$drivers
  windows <- flows |> distinct(period, window_start, window_end)
  vars <- setdiff(names(panel), c("iso3", "year", "pop"))

  # zero-fill: no record in UCDP / EM-DAT = no event
  zf <- vars[matches_any(vars, unlist(dcfg$zero_fill))]
  panel <- panel |> mutate(across(all_of(zf), \(x) coalesce(x, 0)))

  long <- windows |>
    rowwise() |>
    mutate(year = list(seq(window_start, window_end))) |>
    unnest(year) |>
    inner_join(panel, by = "year", relationship = "many-to-many")

  agg_one <- function(x, rule) {
    if (all(is.na(x))) return(NA_real_)
    switch(rule, mean = mean(x, na.rm = TRUE), max = max(x, na.rm = TRUE),
           min = min(x, na.rm = TRUE), sum = sum(x, na.rm = TRUE))
  }
  rule_for <- function(v) dcfg$aggregation[[v]] %||% dcfg$aggregation$default

  out <- long |>
    group_by(iso3, period) |>
    summarise(pop = mean(pop, na.rm = TRUE),
              across(all_of(vars), \(x) agg_one(x, rule_for(cur_column()))),
              .groups = "drop") |>
    mutate(pop = ifelse(is.nan(pop), NA, pop))

  pc <- vars[matches_any(vars, unlist(dcfg$per_capita))]
  out |> mutate(across(all_of(pc), \(x) x / pop))
}
