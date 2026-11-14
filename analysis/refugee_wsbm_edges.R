# =============================================================================
# refugee_wsbm_edges.R — write the refugee edge list for analysis/refugee_wsbm.py
#
#   Rscript analysis/refugee_wsbm_edges.R [config.yml] [dataset]
#   dataset: unhcr (default) or unhcr_displacement
#
# Writes <output_dir>/<dataset>_wsbm/edges.csv: origin, destination, period, flow
# (flows of at least one person, no self-ties).
# =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(readr) })
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
a <- commandArgs(trailingOnly = TRUE)
cfg <- read_config(if (length(a)) a[1] else "config.yml")
ds <- if (length(a) > 1) a[2] else "unhcr"
fl <- read_flows(cfg, ds) |> filter(origin != destination, flow >= 1) |> select(origin, destination, period, flow)
write_csv(fl, out_path(cfg, paste0(ds, "_wsbm"), "edges.csv"))
msg("Wrote ", nrow(fl), " flows for ", ds)
