#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
})

source("functions/_source.R")
source("functions/3-effect-size-factorial.R")

args <- commandArgs(trailingOnly = TRUE)
proj_root <- .alinv_project_root()
output_root <- file.path(proj_root, "output")

latest_date_dir <- alinv_resolve_output_date_dir(
  output_root = output_root,
  requested_date = args[[1]] %||% NULL
)
scenario_dirs <- list.dirs(latest_date_dir, recursive = FALSE, full.names = TRUE)

if (!length(scenario_dirs)) {
  stop("No scenario folders found under ", latest_date_dir, call. = FALSE)
}

shared_y_limits <- alinv_temporal_effect_y_limits()

for (scenario_dir in scenario_dirs) {
  model_dir <- file.path(scenario_dir, "data", "model-factorial-effect")
  if (!dir.exists(model_dir)) next

  effect_files <- list.files(
    model_dir,
    pattern = "-effects\\.csv$",
    full.names = TRUE
  )
  if (!length(effect_files)) next

  plot_dir <- file.path(model_dir, "plots_from_csv")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  message(
    "\nScenario: ", basename(scenario_dir),
    " | fixed y-limits = [", shared_y_limits[[1]], ", ", shared_y_limits[[2]], "]"
  )

  for (effect_file in effect_files) {
    meta_file <- sub("-effects\\.csv$", "-meta.csv", effect_file)

    effects_df <- readr::read_csv(effect_file, show_col_types = FALSE)
    meta_df <- if (file.exists(meta_file)) {
      readr::read_csv(meta_file, show_col_types = FALSE)
    } else {
      tibble::tibble()
    }

    plot_obj <- plot_temporal_effects_from_csv(
      effects_df = effects_df,
      meta_df = meta_df,
      y_limits = shared_y_limits
    )

    out_png <- file.path(
      plot_dir,
      paste0(tools::file_path_sans_ext(basename(effect_file)), ".png")
    )
    ggsave(out_png, plot_obj, width = 9, height = 5.6, dpi = 300)
    message("Saved: ", out_png)
  }
}
