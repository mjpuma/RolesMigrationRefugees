# =============================================================================
# utils.R — small helpers shared by every step
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(stringr)
  library(ggplot2)
})

msg <- function(...) message(format(Sys.time(), "%H:%M:%S"), "  ", ...)

# Read config.yml and resolve relative paths against the project folder
read_config <- function(path = "config.yml") {
  cfg <- yaml::read_yaml(path)
  cfg$root <- normalizePath(dirname(path))
  cfg
}

in_root <- function(cfg, p) file.path(cfg$root, p)

out_path <- function(cfg, ...) {
  p <- file.path(cfg$root, cfg$output_dir, ...)
  dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
  p
}

# ---- Country names -> ISO3 ---------------------------------------------------
# One lookup table (data_ref/country_name_crosswalk.csv, columns alias, iso3)
# replaces the hundreds of hand-written renaming lines in the original scripts.
load_crosswalk <- function(cfg) {
  read_csv(in_root(cfg, cfg$country_crosswalk), show_col_types = FALSE,
           col_types = cols(.default = "c")) |>
    distinct(alias, .keep_all = TRUE)
}

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
# dataset (e.g. SDN -> SUD), per cfg$iso3_successors.
add_successor_rows <- function(panel, cfg) {
  succ <- cfg$iso3_successors
  if (is.null(succ) || !length(succ)) return(panel)
  extra <- imap_dfr(succ, function(modern, historic) {
    panel |> filter(iso3 == modern) |> mutate(iso3 = historic)
  })
  bind_rows(panel |> filter(!iso3 %in% names(succ)), extra)
}

# Match variable names against a vector of regular expressions
matches_any <- function(x, patterns) {
  if (!length(patterns)) return(rep(FALSE, length(x)))
  Reduce(`|`, lapply(patterns, function(p) grepl(p, x)))
}

# Adjusted Rand index between two partitions (for comparing role solutions)
adjusted_rand <- function(a, b) {
  tab <- table(a, b)
  comb2 <- function(n) n * (n - 1) / 2
  sum_ij <- sum(comb2(tab))
  sum_a  <- sum(comb2(rowSums(tab)))
  sum_b  <- sum(comb2(colSums(tab)))
  expected <- sum_a * sum_b / comb2(sum(tab))
  (sum_ij - expected) / (0.5 * (sum_a + sum_b) - expected)
}
