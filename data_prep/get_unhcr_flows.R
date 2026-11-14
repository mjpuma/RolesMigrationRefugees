# =============================================================================
# get_unhcr_flows.R — download UNHCR's current forced-displacement flow data
#
#   Rscript data_prep/get_unhcr_flows.R           installs UNHCR's package, then reads it
#   Rscript data_prep/get_unhcr_flows.R github    reads the data file from UNHCR's GitHub
#                                                 repository (no package installation)
#
# Source: UNHCR "Forcibly displaced flow data" (Refugee Data Finder), distributed
# in UNHCR's `refugees` R package (dataset `flows`). One row per country of
# origin x country of asylum x year, with counts of NEW arrivals/applications:
#   refugees          refugee arrivals: group-recognised (prima facie) refugees and
#                     people under temporary protection; refugee-like situations
#                     from 2007. (This is what the original analysis used.)
#   asylum_seekers    new individual asylum applications (limited to some
#                     high-income countries before 2000; broad coverage from 2000)
#   refugee_like      people in refugee-like situations (not included in `refugees`)
#   oip               other people in need of international protection (from 2018;
#                     Venezuelans)
# From the package's `solutions` table, joined by origin, asylum country and year:
#   resettlement      refugees resettled to a third country (coa = resettlement country)
#   returned_refugees refugees returning to their country of origin
# These two are kept for reference and are not used as flows by default: a
# resettlement is counted by origin and resettlement country, and a return runs
# in the opposite direction.
#
# No double counting between `refugees` and `asylum_seekers`: the first counts
# group (prima facie) recognitions and temporary protection, the second counts
# new individual applications.
# Explainer: https://www.unhcr.org/refugee-statistics/insights/explainers/forcibly-displaced-flow-data.html
#
# Writes:
#   data/raw/unhcr_flows_latest.csv        the flow table used by config.yml
#   data/raw/unhcr_flows_latest_README.txt package version, download date, years
# =============================================================================

# 1. Get the package data: install or update UNHCR's package (UNHCR's R-universe
#    has the newest release; CRAN can lag behind), or read the .rda files from GitHub
from_github <- identical(commandArgs(trailingOnly = TRUE)[1], "github")
if (from_github) {
  base <- Sys.getenv("UNHCR_REFUGEES_BASE", "https://raw.githubusercontent.com/PopulationStatistics/refugees/main/")
  get_rda <- function(name) {
    tmp <- tempfile(fileext = ".rda"); download.file(paste0(base, "data/", name, ".rda"), tmp, mode = "wb", quiet = TRUE)
    e <- new.env(); load(tmp, envir = e); e[[name]]
  }
  flows <- get_rda("flows"); solutions <- get_rda("solutions")
  version <- sub("Version: ", "", grep("^Version:", readLines(paste0(base, "DESCRIPTION")), value = TRUE))
  version <- paste(version, "(GitHub main branch)")
} else {
  repos <- c("https://unhcrverse.r-universe.dev", "https://cloud.r-project.org")
  install.packages("refugees", repos = repos)
  flows <- refugees::flows; solutions <- refugees::solutions
  version <- as.character(packageVersion("refugees"))
}
message("refugees package version ", version, ": ", nrow(flows), " rows, years ",
        min(flows$year), "-", max(flows$year))

# 2. Keep the identifying columns and the flow measures; add resettlement and returns
keep <- c("year", "coo_name", "coo_iso", "coa_name", "coa_iso",
          "refugees", "asylum_seekers", "refugee_like", "oip")
missing <- setdiff(keep, names(flows))
if (length(missing)) stop("Columns not found in the flows table: ", paste(missing, collapse = ", "),
                          "\nThe package structure may have changed; check names(refugees::flows).")
out <- as.data.frame(flows)[, keep]
sol <- as.data.frame(solutions)[, c("year", "coo_name", "coo_iso", "coa_name", "coa_iso", "resettlement", "returned_refugees")]
out <- merge(out, sol, by = c("year", "coo_name", "coo_iso", "coa_name", "coa_iso"), all = TRUE)

# 3. Drop rows without a valid country code on either side (unknown origin,
#    stateless, "various") — they cannot be placed in a country network
#    (UKN = unknown, STA = stateless, TIB = Tibetan: UNHCR codes, not countries)
valid <- function(x) !is.na(x) & grepl("^[A-Z]{3}$", x) & !x %in% c("UKN", "UNK", "STA", "TIB", "XXA", "XXB", "XXC")
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
  paste("refugee_like:     ", cover("refugee_like")),
  paste("oip:              ", cover("oip")),
  paste("resettlement:     ", cover("resettlement")),
  paste("returned_refugees:", cover("returned_refugees"))
), "data/raw/unhcr_flows_latest_README.txt")

message("Wrote data/raw/unhcr_flows_latest.csv")
message(paste(readLines("data/raw/unhcr_flows_latest_README.txt"), collapse = "\n"))
