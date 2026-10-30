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
    # `value` can name one column or several; several are summed
    # (e.g. refugee arrivals + asylum applications). A row counts as missing
    # only if all of the named columns are missing.
    vcols <- unlist(spec$value)
    vals  <- as.data.frame(lapply(raw[vcols], as.numeric))
    total <- rowSums(vals, na.rm = TRUE)
    total[rowSums(!is.na(vals)) == 0] <- NA
    df <- raw |>
      transmute(origin      = .data[[spec$origin]],
                destination = .data[[spec$destination]],
                year        = as.integer(.data[[spec$period]]),
                flow        = total)
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

  # optional year range, e.g. years: [2000, 2024]
  if (!is.null(spec$years)) {
    yr <- unlist(spec$years)
    df <- df |> filter(year >= yr[1], year <= yr[2])
  }

  L <- spec$period_length
  df <- df |>
    mutate(period = if (spec$period_anchor == "start") year + L else year,
           window_start = period - L + 1L,
           window_end   = period) |>
    filter(!is.na(origin), !is.na(destination), origin != destination,
           !is.na(flow)) |>
    select(origin, destination, period, window_start, window_end, flow)

  # optional net-flow network (net: true): for each pair of countries keep one
  # directed link, from the net sender to the net receiver, with the net flow
  # (A->B minus B->A). Gross estimation methods that differ mainly in the
  # size of counter-flows (e.g. Abel & Cohen closed demographic accounting,
  # pseudo-Bayesian vs minimisation) give almost identical net flows.
  if (isTRUE(spec$net)) {
    df <- df |>
      mutate(a = pmin(origin, destination), b = pmax(origin, destination),
             signed = if_else(origin == a, flow, -flow)) |>
      group_by(a, b, period, window_start, window_end) |>
      summarise(net = sum(signed), .groups = "drop") |>
      filter(net != 0) |>
      transmute(origin = if_else(net > 0, a, b), destination = if_else(net > 0, b, a),
                period, window_start, window_end, flow = abs(net))
  }

  msg("  ", dataset, ": ", nrow(df), " dyad-periods, periods ",
      paste(sort(unique(df$period)), collapse = ", "))
  df
}
