# Write small, tidy tables containing the scientific values shown in figures.
# Plot styling, annotations, and ggplot rendering details are intentionally not
# included.

alinv_figure_data_dir <- function(project_root) {
  file.path(project_root, "data", "final", "figure_data")
}

alinv_write_figure_data <- function(data, figure_id, table_name, project_root) {
  if (!is.data.frame(data)) {
    stop("Figure data must be supplied as a data frame.", call. = FALSE)
  }
  if (!nzchar(figure_id) || !nzchar(table_name)) {
    stop("figure_id and table_name must be non-empty.", call. = FALSE)
  }
  if (any(vapply(data, is.list, logical(1)))) {
    stop("Figure data cannot contain list columns.", call. = FALSE)
  }

  output_dir <- alinv_figure_data_dir(project_root)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_path <- file.path(output_dir, paste0(figure_id, "_", table_name, ".csv"))
  readr::write_csv(data, output_path, na = "")
  invisible(output_path)
}
