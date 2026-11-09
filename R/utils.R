# =============================================================================
# utils.R — small helpers shared by every step of the pipeline
#
# Every other file in R/ is sourced after this one (run_pipeline.R sources all
# files in R/ in alphabetical order; the helpers below have no dependencies on
# the other files, so the order does not matter).
#
# Contents
#   msg()                 time-stamped progress messages
#   read_config()         read config.yml and remember the project folder
#   in_root(), out_path() build input and output paths from the config
#   load_crosswalk()      country-name -> ISO3 lookup table
#   name_to_iso3()        apply the lookup and REPORT names it cannot match
#   add_successor_rows()  copy drivers of a modern state onto a historical code
#   matches_any()         does a variable name match any of several patterns?
#   adjusted_rand()       agreement between two groupings (adjusted Rand index)
# =============================================================================

# Packages used throughout. suppressPackageStartupMessages() only hides the
# start-up banners; any real loading error still stops the script.
suppressPackageStartupMessages({
  library(dplyr)     # data manipulation verbs (filter, mutate, summarise, joins)
  library(tidyr)     # reshaping (pivot_longer / pivot_wider)
  library(readr)     # fast, explicit CSV reading and writing
  library(purrr)     # map() / imap() style iteration
  library(stringr)   # string helpers
  library(ggplot2)   # figures
})

# Progress message with a clock time, e.g. "14:03:12    REGE iteration 1 ...".
# Written with message(), so it goes to stderr and does not mix with printed
# result tables.
msg <- function(...) message(format(Sys.time(), "%H:%M:%S"), "  ", ...)

# Read config.yml (or any alternative config file) into a nested R list.
# cfg$root records the folder that contains the config file; every relative
# path in the config (data files, output folder) is resolved against it, so the
# pipeline gives the same result whatever the R working directory is.
read_config <- function(path = "config.yml") {
  cfg <- yaml::read_yaml(path)
  # Project root: the nearest folder, starting from the config file's own, that
  # holds run_pipeline.R. The main config.yml sits in the root; variants are in configs/.
  root <- normalizePath(dirname(path))
  while (!file.exists(file.path(root, "run_pipeline.R"))) {
    if (dirname(root) == root) stop("Cannot find the project root (run_pipeline.R) above ", path)
    root <- dirname(root)
  }
  cfg$root <- root
  cfg
}

# Absolute path of an input file named relative to the project folder.
in_root <- function(cfg, p) file.path(cfg$root, p)

# Absolute path of an output file inside the config's output folder, e.g.
# out_path(cfg, "unhcr", "roles_by_country_period.csv")
#   -> <project>/output/unhcr/roles_by_country_period.csv
# The containing folder is created if it does not exist yet, so callers can
# write to the returned path straight away.
out_path <- function(cfg, ...) {
  p <- file.path(cfg$root, cfg$output_dir, ...)
  dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
  p
}

# ---- Country names -> ISO3 ---------------------------------------------------
# One lookup table (data_ref/country_name_crosswalk.csv, columns alias, iso3)
# replaces the hundreds of hand-written renaming lines in the original scripts.
# Sources that identify countries only by name (SPEI files, UCDP-GED, the GDP
# supplement) are matched through it; sources that already carry ISO3 codes
# (Abel & Cohen flows, UNHCR flows, World Bank GDP, EM-DAT) are not.
load_crosswalk <- function(cfg) {
  read_csv(in_root(cfg, cfg$country_crosswalk), show_col_types = FALSE,
           # read every column as text so codes such as "NA" (Namibia's ISO2)
           # or leading zeros are never converted to missing values or numbers
           col_types = cols(.default = "c")) |>
    # if an alias appears twice, the first row wins (keeps the lookup unique)
    distinct(alias, .keep_all = TRUE)
}

# Translate a vector of country names into ISO3 codes.
#   x          names as they appear in the source file
#   crosswalk  table from load_crosswalk()
#   source     label used in the message (e.g. "SPEI")
# Names that are not in the crosswalk at all are listed in a message, so a
# spelling that silently drops a country is visible in the run log. Names that
# ARE in the crosswalk with an empty iso3 (e.g. "Aland", "Siachen Glacier")
# are deliberately dropped without a message: they are territories that do not
# exist in the flow data.
name_to_iso3 <- function(x, crosswalk, source = "") {
  iso <- crosswalk$iso3[match(x, crosswalk$alias)]
  miss <- sort(unique(x[is.na(iso) & !(x %in% crosswalk$alias)]))
  if (length(miss)) {
    msg("  [", source, "] names not in crosswalk (dropped; add to ",
        "country_name_crosswalk.csv if needed): ", paste(miss, collapse = "; "))
  }
  iso
}

# Copy driver rows from a modern code onto a historical code used by a flow
# dataset, per cfg$iso3_successors (e.g. SUD: SDN). Abel & Cohen code Sudan
# before 2011 as "SUD" (incl. South Sudan) and Serbia & Montenegro as "SCG";
# the driver sources only have the modern codes. Without this step those
# country-periods would have no drivers and drop out of the models.
#   panel  country-year driver table (iso3, year, ...)
# The modern country keeps its own rows; the historical code receives a copy.
add_successor_rows <- function(panel, cfg) {
  succ <- cfg$iso3_successors
  if (is.null(succ) || !length(succ)) return(panel)
  # succ is a named list: names = historical codes, values = modern codes
  extra <- imap_dfr(succ, function(modern, historic) {
    panel |> filter(iso3 == modern) |> mutate(iso3 = historic)
  })
  # drop any rows the panel already had under the historical code, then append
  bind_rows(panel |> filter(!iso3 %in% names(succ)), extra)
}

# TRUE for every element of x that matches at least one of the regular
# expressions in `patterns`. Used to apply config rules such as
# per_capita: ["_deaths$", "_affected$"] or log_transform: ["^gdp_pc$", ...].
matches_any <- function(x, patterns) {
  if (!length(patterns)) return(rep(FALSE, length(x)))
  Reduce(`|`, lapply(patterns, function(p) grepl(p, x)))
}

# Adjusted Rand index (Hubert & Arabie 1985) between two groupings of the same
# items, e.g. roles from two methods. 1 = identical groupings, 0 = agreement
# no better than chance, negative = worse than chance.
#   ARI = (sum_ij C(n_ij,2) - E) / (0.5 [sum_i C(a_i,2) + sum_j C(b_j,2)] - E)
#   E   = sum_i C(a_i,2) * sum_j C(b_j,2) / C(n,2)
# where n_ij is the cross-tabulation of the two groupings, a_i and b_j its row
# and column totals, n the number of items and C(m,2) = m(m-1)/2 the number of
# pairs. Group labels do not need to correspond (role 1 in one solution can be
# role 3 in the other).
adjusted_rand <- function(a, b) {
  tab <- table(a, b)
  comb2 <- function(n) n * (n - 1) / 2
  sum_ij <- sum(comb2(tab))
  sum_a  <- sum(comb2(rowSums(tab)))
  sum_b  <- sum(comb2(colSums(tab)))
  expected <- sum_a * sum_b / comb2(sum(tab))
  (sum_ij - expected) / (0.5 * (sum_a + sum_b) - expected)
}
