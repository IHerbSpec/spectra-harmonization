# From the repository root, in RStudio:
# source("analyses/run_all.R")
# Terminal: Rscript analyses/run_all.R [analyses/config/local.R]
# To select a config when sourcing: options(iherbspec.config="config/local.R")
local({
  args <- commandArgs(trailingOnly = TRUE)
  root <- if (file.exists("analyses/auxiliary/io.R")) "analyses" else 
    if (file.exists("auxiliary/io.R")) "." else
      stop("Run from the repository root or analyses/ directory.")
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  config_file <- getOption("iherbspec.config", "config/analysis.R")
  if (length(args) && any(grepl("^--file=.*run_all[.]R$", commandArgs()))) config_file <- args[1]
  # Also accepts analyses/config/local.R from the repository root.
  if (!file.exists(config_file)) config_file <- file.path(root, config_file)
  config_file <- normalizePath(config_file, winslash = "/", mustWork = TRUE)
  source(config_file, local = TRUE)
  for (f in c("auxiliary/io.R", "auxiliary/preprocessing.R",
              "auxiliary/characterization.R", "plotting/plotting_functions_spectra.R"))
    source(file.path(root, f), local = TRUE)
  config <- validate_config(config)
  ctx <- list(root = root, config = config, config_file = config_file)
  for (f in c("R/01_import_preprocess.R", "R/02_characterize.R", "R/03_plot.R")) {
    message("Running ", f)
    source(file.path(root, f), local = TRUE)
  }
  message("Complete: ", ctx$out)
})
