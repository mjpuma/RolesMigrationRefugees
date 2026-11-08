# =============================================================================
# flows.R — read any origin-destination flow table into one standard format
#
# Every flow dataset (Abel & Cohen migrant flows, UNHCR refugee arrivals, or a
# new one) is described by a block under `flows:` in config.yml. read_flows()
# turns the raw file into the same columns, so nothing downstream needs to know
# where the data came from.
#
# Output columns (one row per origin, destination, period):
#   origin, destination  ISO3 codes
#   period               label year = END year of the interval the flow covers
#   window_start/_end    first and last calendar year used to aggregate drivers
#   flow                 number of people
#
# Period convention. Abel & Cohen report flows for five-year intervals labelled
# by their FIRST year (year0 = 1990 means mid-1990 to mid-1995). We label every
# period by its END year (1995) and aggregate the drivers over 1991-1995
# (window_start = period - length + 1). UNHCR flows are annual and labelled by
# the year itself. The config fields period_length and period_anchor encode
# this for each dataset.
# =============================================================================

read_flows <- function(cfg, dataset) {
  # the config block for this dataset, e.g. cfg$flows$undesa
  spec <- cfg$flows[[dataset]]
  if (is.null(spec)) stop("Dataset '", dataset, "' is not defined in config.yml")
  # guess_max = 1e5: let readr look at many rows before deciding column types,
  # so a column with blanks at the top is not mistaken for text
  raw <- read_csv(in_root(cfg, spec$file), show_col_types = FALSE, guess_max = 1e5)

  if (spec$format == "long") {
    # LONG layout: one row per origin, destination and period.
    # `value` can name one column or several; several are summed
    # (e.g. refugee arrivals + asylum applications). A row counts as missing
    # only if all of the named columns are missing.
    vcols <- unlist(spec$value)
    vals  <- as.data.frame(lapply(raw[vcols], as.numeric))
    total <- rowSums(vals, na.rm = TRUE)
    total[rowSums(!is.na(vals)) == 0] <- NA       # all named columns empty -> NA, not 0
    df <- raw |>
      transmute(origin      = .data[[spec$origin]],
                destination = .data[[spec$destination]],
                year        = as.integer(.data[[spec$period]]),
                flow        = total)
  } else if (spec$format == "wide") {
    # WIDE layout: one row per origin-destination pair and one column per year,
    # e.g. refflows.1990 ... refflows.2020 in the UNHCR extract. The columns
    # are found by their common prefix and reshaped to long format.
    value_cols <- names(raw)[startsWith(names(raw), spec$value_prefix)]
    df <- raw |>
      select(origin = all_of(spec$origin), destination = all_of(spec$destination),
             all_of(value_cols)) |>
      pivot_longer(all_of(value_cols), names_to = "year", values_to = "flow") |>
      mutate(year = as.integer(str_remove(year, fixed(spec$value_prefix))),   # "refflows.1990" -> 1990
             flow = as.numeric(flow))
  } else {
    stop("format must be 'long' or 'wide'")
  }

  # optional year range, e.g. years: [2000, 2024] (applied to the raw year,
  # before the period label is computed)
  if (!is.null(spec$years)) {
    yr <- unlist(spec$years)
    df <- df |> filter(year >= yr[1], year <= yr[2])
  }

  # Period label and driver window (see the header). For Abel & Cohen
  # (anchor "start", length 5): year0 1990 -> period 1995, window 1991-1995.
  # For UNHCR (anchor "end", length 1): year 2005 -> period 2005, window 2005-2005.
  L <- spec$period_length
  df <- df |>
    mutate(period = if (spec$period_anchor == "start") year + L else year,
           window_start = period - L + 1L,
           window_end   = period) |>
    # drop rows that cannot be used: unknown origin or destination, flows of a
    # country to itself, and missing flow values (0 flows are kept here; ties
    # are defined later by a threshold)
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
      # put each pair in a fixed order (a < b alphabetically) and give the flow
      # a sign: positive if it goes a -> b, negative if it goes b -> a
      mutate(a = pmin(origin, destination), b = pmax(origin, destination),
             signed = if_else(origin == a, flow, -flow)) |>
      group_by(a, b, period, window_start, window_end) |>
      summarise(net = sum(signed), .groups = "drop") |>     # net = flow(a->b) - flow(b->a)
      filter(net != 0) |>                                    # exactly balanced pairs have no net link
      # direction of the net link follows the sign; the flow is its absolute size
      transmute(origin = if_else(net > 0, a, b), destination = if_else(net > 0, b, a),
                period, window_start, window_end, flow = abs(net))
  }

  msg("  ", dataset, ": ", nrow(df), " dyad-periods, periods ",
      paste(sort(unique(df$period)), collapse = ", "))
  df
}
