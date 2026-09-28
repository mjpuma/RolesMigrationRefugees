# =============================================================================
# flows.R — read any origin-destination flow table into one standard format
#
# Output columns (one row per origin, destination, period):
#   origin, destination  ISO3 codes
#   period               label year = END year of the interval the flow covers
#   window_start/_end    first and last calendar year used to aggregate drivers
#   flow                 number of people
# =============================================================================

read_flows <- function(cfg, dataset) {
  spec <- cfg$flows[[dataset]]
  if (is.null(spec)) stop("Dataset '", dataset, "' is not defined in config.yml")
  raw <- read_csv(in_root(cfg, spec$file), show_col_types = FALSE, guess_max = 1e5)

  if (spec$format == "long") {
    df <- raw |>
      transmute(origin      = .data[[spec$origin]],
                destination = .data[[spec$destination]],
                year        = as.integer(.data[[spec$period]]),
                flow        = as.numeric(.data[[spec$value]]))
  } else if (spec$format == "wide") {
    value_cols <- names(raw)[startsWith(names(raw), spec$value_prefix)]
    df <- raw |>
      select(origin = all_of(spec$origin), destination = all_of(spec$destination),
             all_of(value_cols)) |>
      pivot_longer(all_of(value_cols), names_to = "year", values_to = "flow") |>
      mutate(year = as.integer(str_remove(year, fixed(spec$value_prefix))),
             flow = as.numeric(flow))
  } else {
    stop("format must be 'long' or 'wide'")
  }

  L <- spec$period_length
  df <- df |>
    mutate(period = if (spec$period_anchor == "start") year + L else year,
           window_start = period - L + 1L,
           window_end   = period) |>
    filter(!is.na(origin), !is.na(destination), origin != destination,
           !is.na(flow)) |>
    select(origin, destination, period, window_start, window_end, flow)

  msg("  ", dataset, ": ", nrow(df), " dyad-periods, periods ",
      paste(sort(unique(df$period)), collapse = ", "))
  df
}
