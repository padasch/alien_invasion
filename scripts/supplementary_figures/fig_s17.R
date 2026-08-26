#!/usr/bin/env Rscript

# Exploratory supplementary SWC path decomposition. Figure 6 reports the
# path-summed total c' + ab; this figure separates that total into the
# SWC-adjusted component c' and the SWC-associated component ab. Repeated-
# response models use same-day or latest preceding measured SWC within seven
# days. All inferential cells use a block-stratified container-cluster bootstrap.

options(stringsAsFactors = FALSE)

arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(arg)) stop("Run this file with Rscript.", call. = FALSE)
script_file <- normalizePath(
  gsub("~\\+~", " ", sub("^--file=", "", arg[[1]])),
  winslash = "/",
  mustWork = TRUE
)
project_root_override <- Sys.getenv("ALINV_PROJECT_ROOT_OVERRIDE", unset = "")
project_root <- if (nzchar(project_root_override)) {
  normalizePath(project_root_override, winslash = "/", mustWork = TRUE)
} else {
  normalizePath(file.path(dirname(script_file), "..", ".."), winslash = "/", mustWork = TRUE)
}
setwd(project_root)

renv_lib <- Sys.glob(file.path(project_root, "renv", "library", "*", "R-*", "*"))
if (length(renv_lib)) .libPaths(c(normalizePath(renv_lib[[1]], winslash = "/", mustWork = TRUE), .libPaths()))

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(tidyr)
})
suppressMessages(suppressPackageStartupMessages({
  function_dir <- file.path(project_root, "scripts", "functions")
  source(file.path(function_dir, "bootstrap.R"))
  source(file.path(function_dir, "figure_data.R"))
}))

output_dir <- Sys.getenv(
  "ALINV_SUPPLEMENTARY_SEM_OUTPUT_DIR",
  unset = file.path(project_root, "output", "supplementary")
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_filename <- Sys.getenv(
  "ALINV_SUPPLEMENTARY_SEM_OUTPUT_FILENAME",
  unset = "fig_s17.pdf"
)
repeated_effects_file <- Sys.getenv(
  "ALINV_REPEATED_SEM_EFFECTS_FILE",
  unset = ""
)

species_order <- c("fagus", "quercus")
species_labels <- c(fagus = "Fagus", quercus = "Quercus")
treatment_order <- c("precipitation", "robinia", "culture", "extreme_event")
treatment_labels <- c(
  precipitation = "Precipitation",
  robinia = "Robinia",
  culture = "Culture",
  extreme_event = "Extreme event"
)
response_order <- c(
  "Leaf-out", "Volume (total)", "Volume (phase incr.)", "Chlorophyll",
  "Vitality", "Quantum yield", "Senescence (%)", "Senescence (Chl)"
)
component_order <- c("SWC-adjusted (c')", "Via SWC (ab)")

format_standard_effects <- function(x) {
  x %>%
  dplyr::filter(.data$component %in% c("direct", "indirect")) %>%
  dplyr::transmute(
    species = .data$species,
    response_label = dplyr::recode(
      .data$response_label,
      `Volume (incr.)` = "Volume (phase incr.)"
    ),
    treatment = .data$treatment,
    component = dplyr::recode(
      .data$component,
      direct = "SWC-adjusted (c')",
      indirect = "Via SWC (ab)"
    ),
    estimate = .data$estimate,
    lower = .data$lower,
    upper = .data$upper,
    p_boot = .data$p_boot,
    n_boot_success = .data$n_boot
  )
}

baseline_standard <- alinv_read_past_only_7d_sem_bootstrap_effects(project_root) %>%
  format_standard_effects()

standard <- if (nzchar(repeated_effects_file)) {
  if (!file.exists(repeated_effects_file)) {
    stop("Configured repeated-response SEM effects file does not exist: ", repeated_effects_file)
  }
  readr::read_csv(repeated_effects_file, show_col_types = FALSE) %>%
    format_standard_effects()
} else {
  baseline_standard
}

phenology_bundle <- alinv_read_phenology_sem_bootstrap_effects(project_root)
# The cached phenology SEM uses the display-oriented timing index: negative is
# later and positive is earlier leaf-out.
phenology_di <- phenology_bundle$effects %>%
  dplyr::filter(.data$metric %in% c("direct", "indirect")) %>%
  dplyr::transmute(
    species = .data$species,
    response_label = "Leaf-out",
    treatment = .data$effect,
    component = dplyr::recode(
      .data$metric,
      direct = "SWC-adjusted (c')",
      indirect = "Via SWC (ab)"
    ),
    estimate = .data$estimate,
    lower = .data$lower_95,
    upper = .data$upper_95,
    p_boot = .data$p_boot,
    n_boot_success = .data$n_boot_success
  )
observed <- dplyr::bind_rows(standard, phenology_di)
if (any(observed$n_boot_success < ALINV_BOOTSTRAP_TARGET, na.rm = TRUE)) {
  stop("One or more SEM cells have fewer than 1,000 successful bootstrap refits.", call. = FALSE)
}

plot_data <- tidyr::expand_grid(
  species = species_order,
  response_label = response_order,
  treatment = treatment_order,
  component = component_order
) %>%
  dplyr::left_join(observed, by = c("species", "response_label", "treatment", "component")) %>%
  dplyr::mutate(
    species_label = unname(species_labels[.data$species]),
    treatment_label = factor(
      unname(treatment_labels[.data$treatment]),
      levels = unname(treatment_labels[treatment_order])
    ),
    response_label = factor(.data$response_label, levels = rev(response_order)),
    component = factor(.data$component, levels = component_order),
    significant = is.finite(.data$p_boot) & .data$p_boot < 0.05,
    label = dplyr::if_else(.data$significant, sprintf("%.2f", .data$estimate), "")
  )

fill_reference <- if (nzchar(repeated_effects_file)) {
  dplyr::bind_rows(baseline_standard, phenology_di)
} else {
  observed
}
fill_limit <- max(abs(fill_reference$estimate), na.rm = TRUE)
fill_limit <- ceiling(fill_limit * 10) / 10

theme_sem <- function(show_y = TRUE) {
  ggplot2::theme_classic(base_size = 7.4) +
    ggplot2::theme(
      text = ggplot2::element_text(color = "black"),
      axis.title = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1, vjust = 1, size = 6.6),
      axis.text.y = if (show_y) ggplot2::element_text(size = 6.8) else ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "plain", size = 8.2, margin = ggplot2::margin(b = 2)),
      plot.margin = ggplot2::margin(2, 2, 2, 2),
      legend.position = "bottom"
    )
}

build_panel <- function(species_i, component_i, show_y = TRUE) {
  df <- plot_data %>%
    dplyr::filter(.data$species == species_i, as.character(.data$component) == component_i) %>%
    dplyr::mutate(
      missing = !is.finite(.data$estimate),
      text_color = dplyr::if_else(
        .data$significant & abs(.data$estimate) >= 0.62 * fill_limit,
        "white",
        "black"
      )
    )

  ggplot2::ggplot(df, ggplot2::aes(.data$treatment_label, .data$response_label)) +
    ggplot2::geom_tile(
      ggplot2::aes(fill = .data$estimate),
      width = 1.01,
      height = 1.01,
      color = NA,
      linewidth = 0
    ) +
    ggplot2::geom_text(
      data = dplyr::filter(df, .data$significant),
      ggplot2::aes(label = .data$label, color = .data$text_color),
      size = 1.8,
      show.legend = FALSE
    ) +
    ggplot2::scale_color_identity() +
    ggplot2::scale_fill_gradient2(
      low = "#D65F5F",
      mid = "white",
      high = "#3C6E8F",
      midpoint = 0,
      limits = c(-fill_limit, fill_limit),
      na.value = "#D9D9D9",
      oob = scales::squish,
      name = "Effect size",
      guide = ggplot2::guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        barwidth = grid::unit(28, "mm"),
        barheight = grid::unit(2.5, "mm")
      )
    ) +
    ggplot2::scale_x_discrete(drop = FALSE) +
    ggplot2::scale_y_discrete(drop = FALSE) +
    ggplot2::labs(title = paste0(species_labels[[species_i]], ": ", component_i)) +
    theme_sem(show_y)
}

panels <- unlist(lapply(species_order, function(species_i) {
  list(
    build_panel(species_i, "SWC-adjusted (c')", TRUE),
    build_panel(species_i, "Via SWC (ab)", FALSE)
  )
}), recursive = FALSE)

figure <- patchwork::wrap_plots(panels, ncol = 2, guides = "collect") +
  patchwork::plot_layout(heights = c(1, 1.16)) &
  ggplot2::theme(legend.position = "bottom")

pdf_file <- file.path(output_dir, output_filename)
ggplot2::ggsave(
  pdf_file,
  figure,
  width = 138,
  height = 132,
  units = "mm",
  bg = "white",
  limitsize = FALSE
)
alinv_write_figure_data(
  plot_data %>%
    dplyr::transmute(
      species = .data$species,
      response = as.character(.data$response_label),
      treatment = .data$treatment,
      component = as.character(.data$component),
      estimate = .data$estimate,
      lower_95 = .data$lower,
      upper_95 = .data$upper,
      p_value = .data$p_boot,
      significant = .data$significant
    ),
  figure_id = tools::file_path_sans_ext(basename(output_filename)),
  table_name = "heatmap",
  project_root = project_root
)
message("Created: ", normalizePath(pdf_file, winslash = "/", mustWork = TRUE))
