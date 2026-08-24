#!/usr/bin/env Rscript

# Recreate Figure 6 and Supplementary Figure S3 using the past-only measured
# SWC sensitivity (same-day or latest preceding SWC within seven days).
# Original figure files are not overwritten. The phenology SEM is unchanged.

options(stringsAsFactors = FALSE)

arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(arg)) stop("Run this file with Rscript.", call. = FALSE)
script_file <- normalizePath(
  gsub("~[+]~", " ", sub("^--file=", "", arg[[1]])),
  winslash = "/", mustWork = TRUE
)
project_root <- normalizePath(
  file.path(dirname(script_file), "..", "..", "..", "..", ".."),
  winslash = "/", mustWork = TRUE
)
branch_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
output_dir <- file.path(branch_dir, "recreated_heatmaps")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

past_effects_file <- file.path(output_dir, "..", "output", "past-only-sem-bootstrap-effects.csv")
past_totals_file <- file.path(output_dir, "..", "output", "past-only-figure6-ready-total-effects.csv")
if (!file.exists(past_effects_file) || !file.exists(past_totals_file)) {
  stop("Past-only B=1000 SEM outputs are missing. Run 01_past_only_7d_sem_bootstrap.R first.")
}

# Load the publication Figure 6 functions and defaults.
source(file.path(project_root, "scripts", "main_figures", "00_setup.R"))
source(file.path(project_root, "scripts", "main_figures", "06_fig6_volume_sem.R"))

past_effects <- readr::read_csv(past_effects_file, show_col_types = FALSE)
past_totals <- readr::read_csv(past_totals_file, show_col_types = FALSE)
if (any(past_effects$n_boot != 1000L) || any(past_totals$n_boot_success != 1000L)) {
  stop("Past-only heatmaps require 1,000 successful bootstrap replicates per available cell.")
}

# The new Figure 6 uses past-only totals for repeated responses, while the
# current phenology timing SEM remains unchanged. Retaining the original
# Figure 6 fill limit makes colors directly comparable between PDFs.
past_fig6_totals <- past_totals %>%
  dplyr::transmute(
    species = .data$species,
    treatment = .data$treatment,
    response_label = .data$response_label,
    resp_var = .data$resp_var,
    estimate = .data$estimate,
    p_value = .data$p_boot,
    lower = .data$lower,
    upper = .data$upper,
    n_boot_success = .data$n_boot_success,
    source_file = normalizePath(past_totals_file, winslash = "/", mustWork = TRUE)
  )

baseline_fig6_sem <- prepare_fig6_sem_data()
baseline_fig6_fill_limit <- max(abs(baseline_fig6_sem$estimate), na.rm = TRUE)
if (!is.finite(baseline_fig6_fill_limit) || baseline_fig6_fill_limit <= 0) {
  stop("Could not recover the baseline Figure 6 heatmap color limit.")
}

ALINV_FINAL_OUTPUT_DIR <- output_dir
fig6_file <- make_fig6(
  repeated_sem_totals = past_fig6_totals,
  heatmap_fill_limit = baseline_fig6_fill_limit,
  output_file = "fig6-volume-sem-past-only-swc.pdf"
)

# Reuse the supplementary plotting script through its input/output overrides.
# That script uses the baseline S3 color limit when a sensitivity input is set.
Sys.setenv(
  ALINV_PROJECT_ROOT_OVERRIDE = project_root,
  ALINV_REPEATED_SEM_EFFECTS_FILE = normalizePath(past_effects_file, winslash = "/", mustWork = TRUE),
  ALINV_SUPPLEMENTARY_SEM_OUTPUT_DIR = normalizePath(output_dir, winslash = "/", mustWork = TRUE),
  ALINV_SUPPLEMENTARY_SEM_OUTPUT_FILENAME = "figure-s3-sem-path-decomposition-past-only-swc.pdf"
)
on.exit(
  Sys.unsetenv(c(
    "ALINV_PROJECT_ROOT_OVERRIDE",
    "ALINV_REPEATED_SEM_EFFECTS_FILE",
    "ALINV_SUPPLEMENTARY_SEM_OUTPUT_DIR",
    "ALINV_SUPPLEMENTARY_SEM_OUTPUT_FILENAME"
  )),
  add = TRUE
)
source(file.path(
  project_root, "scripts", "supplementary_figures",
  "03-sem-path-decomposition-heatmaps.R"
))

s3_file <- file.path(output_dir, "figure-s3-sem-path-decomposition-past-only-swc.pdf")
if (!file.exists(fig6_file) || !file.exists(s3_file)) {
  stop("One or more past-only SEM heatmap PDFs were not created.")
}

# Confirm that the repeated-response totals used in Figure 6 are exactly the
# path sums used in the past-only decomposition.
path_qa <- past_effects %>%
  dplyr::filter(.data$component %in% c("direct", "indirect", "total")) %>%
  dplyr::select("species", "resp_var", "response_label",
                "treatment", "component", "estimate") %>%
  tidyr::pivot_wider(names_from = "component", values_from = "estimate") %>%
  dplyr::mutate(error = .data$total - (.data$direct + .data$indirect))
if (max(abs(path_qa$error), na.rm = TRUE) > 1e-12) {
  stop("Past-only total != direct + indirect during heatmap QA.")
}

message("Created: ", normalizePath(fig6_file, winslash = "/", mustWork = TRUE))
message("Created: ", normalizePath(s3_file, winslash = "/", mustWork = TRUE))
