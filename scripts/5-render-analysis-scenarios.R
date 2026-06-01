#!/usr/bin/env Rscript

renv_lib <- Sys.glob(file.path("renv", "library", "*", "R-*", "*"))
if (length(renv_lib)) {
  .libPaths(c(normalizePath(renv_lib[[1]], winslash = "/", mustWork = TRUE), .libPaths()))
}

suppressPackageStartupMessages({
  library(rmarkdown)
  library(dplyr)
})

source("functions/_source.R")

args <- commandArgs(trailingOnly = TRUE)
analysis_date <- if (length(args) >= 1) args[[1]] else as.character(Sys.Date())
proj_root <- .alinv_project_root()
output_root <- file.path(proj_root, "output")

notebooks <- c(
  "notebooks/1-treatment-effects.Rmd",
  "notebooks/2-swc-interpolation.Rmd",
  "notebooks/3-sem-aggregation.Rmd",
  "notebooks/4-data-qc.Rmd",
  "notebooks/5-size-trajectories.Rmd"
)

global_notebooks <- c(
  "notebooks/6-swp-provenance.Rmd"
)

scenarios <- alinv_scenario_grid()

for (i in seq_len(nrow(scenarios))) {
  scenario_i <- scenarios[i, ]

  ctx <- alinv_set_analysis_context(
    scenario_label = scenario_i$scenario_label,
    soil_filter = scenario_i$soil_filter,
    include_soil_treatment = scenario_i$include_soil_treatment,
    analysis_date = analysis_date,
    output_root = output_root
  )

  message("\n=== Rendering scenario: ", ctx$scenario_label, " ===")

  for (nb in notebooks) {
    input_file <- file.path(proj_root, nb)
    output_file <- file.path(
      ctx$notebooks_root,
      paste0(tools::file_path_sans_ext(basename(nb)), ".html")
    )

    message("Rendering ", basename(nb), " -> ", output_file)

    rmarkdown::render(
      input = input_file,
      output_file = output_file,
      knit_root_dir = proj_root,
      params = list(
        scenario_label = ctx$scenario_label,
        soil_filter = ctx$soil_filter,
        include_soil_treatment = ctx$include_soil_treatment,
        analysis_date = ctx$analysis_date,
        output_root = output_root
      ),
      envir = new.env(parent = globalenv()),
      quiet = FALSE
    )
  }
}

global_ctx <- alinv_set_analysis_context(
  scenario_label = "global",
  soil_filter = "both",
  include_soil_treatment = FALSE,
  analysis_date = analysis_date,
  output_root = output_root
)

message("\n=== Rendering global notebooks ===")

for (nb in global_notebooks) {
  input_file <- file.path(proj_root, nb)
  output_file <- file.path(
    global_ctx$notebooks_root,
    paste0(tools::file_path_sans_ext(basename(nb)), ".html")
  )

  message("Rendering ", basename(nb), " -> ", output_file)

  rmarkdown::render(
    input = input_file,
    output_file = output_file,
    knit_root_dir = proj_root,
    params = list(
      scenario_label = global_ctx$scenario_label,
      soil_filter = global_ctx$soil_filter,
      include_soil_treatment = global_ctx$include_soil_treatment,
      analysis_date = global_ctx$analysis_date,
      output_root = output_root
    ),
    envir = new.env(parent = globalenv()),
    quiet = FALSE
  )
}

message("\nScenario rendering completed for analysis date ", analysis_date, ".")
