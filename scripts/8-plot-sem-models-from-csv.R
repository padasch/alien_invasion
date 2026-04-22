#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
})

source("functions/_source.R")
source("functions/4-structural-equation-model.R")

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

for (scenario_dir in scenario_dirs) {
  model_dir <- file.path(scenario_dir, "data", "model-sem")
  if (!dir.exists(model_dir)) next

  node_files <- list.files(
    model_dir,
    pattern = "-graph_nodes\\.csv$",
    full.names = TRUE
  )
  if (!length(node_files)) next

  plot_dir <- file.path(model_dir, "plots_from_csv")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  message("\nScenario: ", basename(scenario_dir))

  for (node_file in node_files) {
    stem <- sub("-graph_nodes\\.csv$", "", basename(node_file))
    edge_file <- file.path(model_dir, paste0(stem, "-graph_edges.csv"))
    metrics_file <- file.path(model_dir, paste0(stem, "-graph_metrics.csv"))

    if (!file.exists(edge_file) || !file.exists(metrics_file)) next

    nodes_df <- readr::read_csv(node_file, show_col_types = FALSE)
    edges_df <- readr::read_csv(edge_file, show_col_types = FALSE)
    metrics_df <- readr::read_csv(metrics_file, show_col_types = FALSE)

    plot_obj <- plot_sem_graph_components(
      nodes_df = nodes_df,
      edges_df = edges_df,
      metrics_df = metrics_df
    )

    out_png <- file.path(plot_dir, paste0(stem, ".png"))
    ggsave(out_png, plot_obj, width = 9, height = 5.6, dpi = 300)
    message("Saved: ", out_png)
  }
}
