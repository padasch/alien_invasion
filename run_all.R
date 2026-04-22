#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(rmarkdown)
})

source("functions/_source.R")

# Optional full rebuild from raw inputs.
# The repository already contains the cleaned interim files, so most users
# do not need to run these lines. Uncomment them only when you want to
# recreate the tracked cleaned data from the raw source files.
# source("scripts/1-data-cleaning.R")
# source("scripts/4-impute-swc-gam.R")

analysis_date <- as.character(Sys.Date())
proj_root <- .alinv_project_root()
output_root <- file.path(proj_root, "output")

notebooks <- c(
  "notebooks/1-treatment-effects.Rmd",
  "notebooks/3-sem-aggregation.Rmd"
)

# Slice (4): Run both soils without soil as treatment.
scenarios <- alinv_scenario_grid() |> slice(4)

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

message("\nAll notebook renders completed for analysis date ", analysis_date, ".")
