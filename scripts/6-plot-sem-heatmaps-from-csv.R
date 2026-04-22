#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
})

source("functions/_source.R")
source("functions/4-structural-equation-model.R")

args <- commandArgs(trailingOnly = TRUE)
proj_root <- .alinv_project_root()
output_root <- file.path(proj_root, "output")

resolve_latest_date_dir <- function(output_root, requested_date = NULL) {
  if (!is.null(requested_date) && nzchar(requested_date)) {
    date_dir <- file.path(output_root, requested_date)
    if (!dir.exists(date_dir)) {
      stop("Requested output date folder does not exist: ", date_dir, call. = FALSE)
    }
    return(date_dir)
  }

  date_dirs <- list.dirs(output_root, full.names = TRUE, recursive = FALSE)
  date_tbl <- tibble(
    path = date_dirs,
    label = basename(date_dirs),
    date = suppressWarnings(as.Date(basename(date_dirs)))
  ) %>%
    filter(!is.na(.data$date)) %>%
    arrange(desc(.data$date))

  if (!nrow(date_tbl)) {
    stop("No dated output folders found under ", output_root, call. = FALSE)
  }

  date_tbl$path[[1]]
}

read_heatmap_csv <- function(path) {
  df <- readr::read_csv(path, show_col_types = FALSE)
  row_col <- names(df)[[1]]
  col_order <- names(df)[-1]
  row_order <- df[[row_col]]

  df %>%
    tidyr::pivot_longer(-dplyr::all_of(row_col), names_to = "col_label", values_to = "value") %>%
    dplyr::rename(row_label = !!row_col) %>%
    dplyr::mutate(
      row_label = factor(.data$row_label, levels = row_order),
      col_label = factor(.data$col_label, levels = col_order)
    )
}

compute_limit_from_csvs <- function(csv_files) {
  vals <- unlist(lapply(csv_files, function(path) {
    df <- readr::read_csv(path, show_col_types = FALSE)
    unlist(df[-1], use.names = FALSE)
  }), use.names = FALSE)
  vals <- as.numeric(vals)
  vals <- vals[is.finite(vals)]
  if (!length(vals)) {
    return(1)
  }
  max_abs <- max(abs(vals), na.rm = TRUE)
  if (!is.finite(max_abs) || max_abs <= 0) 1 else max_abs
}

heatmap_specs <- sem_heatmap_specs()
latest_date_dir <- resolve_latest_date_dir(output_root, requested_date = args[[1]] %||% NULL)
scenario_dirs <- list.dirs(latest_date_dir, recursive = FALSE, full.names = TRUE)

if (!length(scenario_dirs)) {
  stop("No scenario folders found under ", latest_date_dir, call. = FALSE)
}

for (scenario_dir in scenario_dirs) {
  heatmap_dir <- file.path(
    scenario_dir,
    "data",
    "structural_equation_models",
    "heatmap_data"
  )
  if (!dir.exists(heatmap_dir)) next

  csv_files <- list.files(heatmap_dir, pattern = "\\.csv$", full.names = TRUE)
  if (!length(csv_files)) next

  plot_dir <- file.path(
    scenario_dir,
    "data",
    "structural_equation_models",
    "heatmap_plots_from_csv"
  )
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  shared_limit <- compute_limit_from_csvs(csv_files)
  message("\nScenario: ", basename(scenario_dir), " | shared limit = ", round(shared_limit, 3))

  for (csv_file in csv_files) {
    base_name <- tools::file_path_sans_ext(basename(csv_file))
    species_name <- sub("-.*$", "", base_name)
    stub <- sub("^[^-]+-", "", base_name)
    panel_title <- heatmap_specs$panel_title[match(stub, heatmap_specs$file_stub)]
    panel_title <- panel_title %||% stub

    panel_df <- read_heatmap_csv(csv_file)
    plot_obj <- plot_sem_heatmap_panel(
      panel_df = panel_df,
      limit = shared_limit,
      title = paste(species_name, "-", panel_title)
    )

    out_png <- file.path(plot_dir, paste0(base_name, ".png"))
    ggsave(out_png, plot_obj, width = 9, height = 5.6, dpi = 300)
    message("Saved: ", out_png)
  }
}
