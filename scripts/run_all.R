#!/usr/bin/env Rscript

# Rebuild the publication figure collections from the repository root:
#
#   Rscript --vanilla scripts/run_all.R

arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- if (length(arg)) {
  normalizePath(
    gsub("~[+]~", " ", sub("^--file=", "", arg[[1]])),
    winslash = "/",
    mustWork = TRUE
  )
} else {
  normalizePath("scripts/run_all.R", winslash = "/", mustWork = TRUE)
}
project_root <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/", mustWork = TRUE)
setwd(project_root)

renv_lib <- Sys.glob(file.path(project_root, "renv", "library", "*", "R-*", "*"))
if (length(renv_lib)) {
  .libPaths(c(normalizePath(renv_lib[[1]], winslash = "/", mustWork = TRUE), .libPaths()))
}

run_script <- function(path, args = character(), env = character()) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  message("Running ", path)
  status <- system2(
    file.path(R.home("bin"), "Rscript"),
    args = c("--vanilla", shQuote(path), args),
    env = env
  )
  if (!identical(status, 0L)) {
    stop("Publication script failed: ", path, call. = FALSE)
  }
}

# Cleaned interim data are tracked. Uncomment only when rebuilding them from
# the raw inputs.
# run_script(file.path(project_root, "scripts", "auxiliary", "1-data-cleaning.R"))
# run_script(file.path(project_root, "scripts", "auxiliary", "3-cleaning-sensor-data.R"))
# run_script(file.path(project_root, "scripts", "auxiliary", "4-impute-swc-gam.R"))

bootstrap_cores <- as.integer(Sys.getenv("ALINV_BOOT_CORES", unset = "8"))
if (!is.finite(bootstrap_cores) || bootstrap_cores < 1L) bootstrap_cores <- 1L

run_script(file.path(project_root, "scripts", "main_figures", "make_all_figures.R"))
run_script(
  file.path(project_root, "scripts", "supplementary_figures", "make-v1-figures.R"),
  env = c(
    "ALINV_SUPP_BOOT_B=1000",
    paste0("ALINV_SUPP_BOOT_CORES=", bootstrap_cores)
  )
)
run_script(file.path(
  project_root, "scripts", "supplementary_figures",
  "fig_s16.R"
))
run_script(file.path(
  project_root, "scripts", "supplementary_figures",
  "fig_s17.R"
))

message("Publication figures rebuilt under output/main_figures and output/supplementary.")
