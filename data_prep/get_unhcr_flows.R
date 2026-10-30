# =============================================================================
# get_unhcr_flows.R — download UNHCR's current forced-displacement flow data
#
#   Rscript data_prep/get_unhcr_flows.R
#
# Source: UNHCR "Forcibly displaced flow data" (Refugee Data Finder), distributed
# in UNHCR's `refugees` R package (dataset `flows`). One row per country of
# origin x country of asylum x year, with counts of NEW arrivals/applications:
#   refugees          refugee arrivals: group-recognised (prima facie) refugees and
#                     people under temporary protection; refugee-like situations
#                     from 2007. (This is what the original analysis used.)
#   asylum_seekers    new individual asylum applications (limited to some
#                     high-income countries before 2000; broad coverage from 2000)
#   oip               other people in need of international protection (from 2018)
#   returned_refugees refugees returning to their country of origin
# Explainer: https://www.unhcr.org/refugee-statistics/insights/explainers/forcibly-displaced-flow-data.html
#
# Writes:
#   data/raw/unhcr_flows_latest.csv        the flow table used by config.yml
#   data/raw/unhcr_flows_latest_README.txt package version, download date, years
# =============================================================================

# 1. Install or update UNHCR's package (UNHCR's R-universe has the newest release;
#    CRAN can lag behind by a year or more)
repos <- c("https://unhcrverse.r-universe.dev", "https://cloud.r-project.org")
install.packages("refugees", repos = repos)

flows <- refugees::flows
version <- as.character(packageVersion("refugees"))
message("refugees package version ", version, ": ", nrow(flows), " rows, years ",
        min(flows$year), "-", max(flows$year))

# 2. Keep the identifying columns and the flow measures
keep <- c("year", "coo_name", "coo_iso", "coa_name", "coa_iso",
          "refugees", "asylum_seekers", "oip", "returned_refugees")
missing <- setdiff(keep, names(flows))
if (length(missing)) stop("Columns not found in refugees::flows: ", paste(missing, collapse = ", "),
                          "\nThe package structure may have changed; check names(refugees::flows).")
out <- as.data.frame(flows)[, keep]

# 3. Drop rows without a valid country code on either side (unknown origin,
#    stateless, "various") — they cannot be placed in a country network
valid <- function(x) !is.na(x) & grepl("^[A-Z]{3}$", x) & !x %in% c("UNK", "STA", "XXA", "XXB", "XXC")
dropped <- out[!(valid(out$coo_iso) & valid(out$coa_iso)), ]
out <- out[valid(out$coo_iso) & valid(out$coa_iso), ]
message("dropped ", nrow(dropped), " rows without a valid origin/asylum ISO3 code (",
        round(sum(dropped$refugees, na.rm = TRUE)), " refugee arrivals, ",
        round(sum(dropped$asylum_seekers, na.rm = TRUE)), " asylum applications)")

# 4. Save
dir.create("data/raw", showWarnings = FALSE, recursive = TRUE)
write.csv(out, "data/raw/unhcr_flows_latest.csv", row.names = FALSE, na = "")

cover <- function(v) {
  y <- out$year[!is.na(out[[v]]) & out[[v]] > 0]
  if (!length(y)) return("none")
  sprintf("%d-%d, %d dyad-years with flow > 0", min(y), max(y), length(y))
}
writeLines(c(
  "UNHCR forced-displacement flow data",
  paste("Source: refugees R package version", version, "(dataset `flows`)"),
  paste("Downloaded:", format(Sys.Date())),
  paste("Rows kept:", nrow(out), "; rows dropped (no valid ISO3):", nrow(dropped)),
  paste("refugees:         ", cover("refugees")),
  paste("asylum_seekers:   ", cover("asylum_seekers")),
  paste("oip:              ", cover("oip")),
  paste("returned_refugees:", cover("returned_refugees"))
), "data/raw/unhcr_flows_latest_README.txt")

message("Wrote data/raw/unhcr_flows_latest.csv")
message(paste(readLines("data/raw/unhcr_flows_latest_README.txt"), collapse = "\n"))
