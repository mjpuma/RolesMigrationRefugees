# =============================================================================
# run_pipeline.R — roles of countries in international mobility networks
#
#   Rscript run_pipeline.R                 # uses config.yml
#   Rscript run_pipeline.R my_config.yml   # alternative settings
#
# Steps for every dataset listed under `run:` in the config:
#   1. read flows          (R/flows.R)
#   2. REGE + clustering   (R/rege.R, R/roles.R)   -> output/<dataset>/roles_*.csv
#   3. aggregate drivers   (R/drivers.R)
#   4. role models         (R/models.R)            -> output/<dataset>/*_models.csv
#   5. figures             (R/figures.R)           -> output/<dataset>/figures/
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
config_file <- if (length(args)) args[1] else "config.yml"

for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
cfg <- read_config(config_file)
set.seed(cfg$seed)

driver_panel <- build_driver_panel(cfg)
results <- list()

for (dataset in cfg$run) {
  msg("==== ", dataset, ": ", cfg$flows[[dataset]]$label, " ====")
  flows <- read_flows(cfg, dataset)
  role_res <- run_roles(cfg, dataset, flows)

  drivers <- aggregate_drivers(cfg, driver_panel, flows)
  write_csv(drivers, out_path(cfg, dataset, "drivers_by_country_period.csv"))

  model_res <- run_models(cfg, dataset, role_res$roles, drivers)
  save_figures(cfg, dataset, role_res, model_res)
  results[[dataset]] <- list(roles = role_res, models = model_res)
}

saveRDS(results, out_path(cfg, "all_results.rds"))
writeLines(capture.output(sessionInfo()), out_path(cfg, "sessionInfo.txt"))
msg("Done. Results in ", file.path(cfg$root, cfg$output_dir))
