#!/usr/bin/env Rscript

renv_lib <- Sys.glob(file.path("renv", "library", "*", "R-*", "*"))
if (length(renv_lib)) {
  .libPaths(c(normalizePath(renv_lib[[1]], winslash = "/", mustWork = TRUE), .libPaths()))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(tibble)
})

source("functions/_source.R")
source("functions/1-summary-figures.R")
source("functions/3-effect-size-factorial.R")
source("functions/4-structural-equation-model.R")
source("functions/8-size-trajectories.R")

analysis_date <- as.character(Sys.Date())
proj_root <- .alinv_project_root()
output_root <- file.path(proj_root, "output")

ctx <- alinv_set_analysis_context(
  scenario_label = "both soils (without soil as treatment)",
  soil_filter = "both",
  include_soil_treatment = FALSE,
  analysis_date = analysis_date,
  output_root = output_root
)

comparison_dir <- alinv_data_path("diagnostics", "volume_calculation_comparison", create_dir = TRUE)
comparison_interim_dir <- file.path(comparison_dir, "interim")
dir.create(comparison_interim_dir, recursive = TRUE, showWarnings = FALSE)

base_growth <- readr::read_csv(file.path(proj_root, "data", "interim", "tree_growth.csv"), show_col_types = FALSE)

growth_compare <- base_growth %>%
  mutate(
    volume_allometry = .data$volume,
    volume_allometry_inc_t0 = .data$volume_inc_t0,
    volume_allometry_inc_phase_abs = .data$volume_inc_phase_abs,
    volume_cylinder = pi * (.data$diameter / 20)^2 * .data$height
  ) %>%
  group_by(.data$tree_id) %>%
  arrange(.data$date, .by_group = TRUE) %>%
  mutate(
    first_volume_cylinder = first(.data$volume_cylinder),
    volume_cylinder_inc_t0 = .data$volume_cylinder - .data$first_volume_cylinder,
    phase_volume_cylinder_baseline = dplyr::case_when(
      .data$phase == "until June"  ~ .data$first_volume_cylinder,
      .data$phase == "July-August" ~ dplyr::last(.data$volume_cylinder[.data$phase == "until June" & !is.na(.data$volume_cylinder)], default = NA_real_),
      .data$phase == "September+"  ~ dplyr::last(.data$volume_cylinder[.data$phase == "July-August" & !is.na(.data$volume_cylinder)], default = NA_real_),
      TRUE ~ NA_real_
    ),
    volume_cylinder_inc_phase_abs = .data$volume_cylinder - .data$phase_volume_cylinder_baseline
  ) %>%
  ungroup()

readr::write_csv(growth_compare, file.path(comparison_interim_dir, "tree_growth.csv"))

original_get_data <- get_data
get_data <- function(type = c("tree", "box"),
                     data_name,
                     with_meta = TRUE,
                     path = "./data/interim",
                     swc_source = "measured") {
  type <- match.arg(type)

  if (identical(type, "tree") && identical(data_name, "growth")) {
    return(
      original_get_data(
        type = type,
        data_name = data_name,
        with_meta = with_meta,
        path = comparison_interim_dir,
        swc_source = swc_source
      )
    )
  }

  original_get_data(
    type = type,
    data_name = data_name,
    with_meta = with_meta,
    path = path,
    swc_source = swc_source
  )
}

comparison_metric_specs <- tibble::tribble(
  ~metric_group, ~abs_resp_var, ~phase_inc_resp_var, ~t0_inc_resp_var, ~abs_label, ~phase_inc_label, ~t0_inc_label, ~abs_panel_title, ~phase_panel_title, ~t0_panel_title, ~abs_y_label, ~phase_y_label, ~t0_y_label,
  "Height", "height", "height_inc_phase_abs", "height_inc_t0", "Height", "Height inc. phase", "Height inc. t0", "Height - Absolute size", "Height - Absolute increment within phase", "Height - Absolute increment from first measurement", "Height (cm)", "Height increment within phase (cm)", "Height increment from first measurement (cm)",
  "Diameter", "diameter", "diameter_inc_phase_abs", "diameter_inc_t0", "Diameter", "Diameter inc. phase", "Diameter inc. t0", "Diameter - Absolute size", "Diameter - Absolute increment within phase", "Diameter - Absolute increment from first measurement", "Diameter (mm)", "Diameter increment within phase (mm)", "Diameter increment from first measurement (mm)",
  "Volume (cylinder)", "volume_cylinder", "volume_cylinder_inc_phase_abs", "volume_cylinder_inc_t0", "Volume cylinder", "Volume cylinder inc. phase", "Volume cylinder inc. t0", "Volume (cylinder) - Absolute size", "Volume (cylinder) - Absolute increment within phase", "Volume (cylinder) - Absolute increment from first measurement", "Volume as cylinder (cm^3)", "Cylinder volume increment within phase (cm^3)", "Cylinder volume increment from first measurement (cm^3)",
  "Volume (allometry)", "volume_allometry", "volume_allometry_inc_phase_abs", "volume_allometry_inc_t0", "Volume allometry", "Volume allometry inc. phase", "Volume allometry inc. t0", "Volume (allometry) - Absolute size", "Volume (allometry) - Absolute increment within phase", "Volume (allometry) - Absolute increment from first measurement", "Volume from allometry (proxy, g)", "Allometric proxy increment within phase (g)", "Allometric proxy increment from first measurement (g)"
)

comparison_sections <- list(
  phase = list(
    key = "phase",
    title = "Phase-wise baseline comparison",
    subtitle = "Absolute values and absolute increments within phase",
    plot_cols = c("abs_resp_var", "phase_inc_resp_var"),
    label_cols = c("abs_label", "phase_inc_label"),
    y_cols = c("abs_y_label", "phase_y_label"),
    heatmap_order = c(
      "height", "height_inc_phase_abs",
      "diameter", "diameter_inc_phase_abs",
      "volume_cylinder", "volume_cylinder_inc_phase_abs",
      "volume_allometry", "volume_allometry_inc_phase_abs"
    ),
    file_stub = "volume-calculation-comparison-phase"
  ),
  t0 = list(
    key = "t0",
    title = "First-measurement baseline comparison",
    subtitle = "Absolute values and absolute increments from the first measurement",
    plot_cols = c("abs_resp_var", "t0_inc_resp_var"),
    label_cols = c("abs_label", "t0_inc_label"),
    y_cols = c("abs_y_label", "t0_y_label"),
    heatmap_order = c(
      "height", "height_inc_t0",
      "diameter", "diameter_inc_t0",
      "volume_cylinder", "volume_cylinder_inc_t0",
      "volume_allometry", "volume_allometry_inc_t0"
    ),
    file_stub = "volume-calculation-comparison-t0"
  ),
  absolute = list(
    key = "absolute",
    title = "Absolute size comparison",
    subtitle = "Absolute values only",
    plot_cols = c("abs_resp_var"),
    label_cols = c("abs_label"),
    y_cols = c("abs_y_label"),
    heatmap_order = c(
      "height",
      "diameter",
      "volume_cylinder",
      "volume_allometry"
    ),
    file_stub = "volume-calculation-comparison-absolute"
  )
)

comparison_resp_label_map <- c(
  height = "Height",
  height_inc_phase_abs = "Height inc. phase",
  height_inc_t0 = "Height inc. t0",
  diameter = "Diameter",
  diameter_inc_phase_abs = "Diameter inc. phase",
  diameter_inc_t0 = "Diameter inc. t0",
  volume_cylinder = "Volume cylinder",
  volume_cylinder_inc_phase_abs = "Volume cylinder inc. phase",
  volume_cylinder_inc_t0 = "Volume cylinder inc. t0",
  volume_allometry = "Volume allometry",
  volume_allometry_inc_phase_abs = "Volume allometry inc. phase",
  volume_allometry_inc_t0 = "Volume allometry inc. t0"
)

comparison_plot_data <- function(resp_var, within_phase = FALSE) {
  df_metric <- get_data("tree", "growth") %>%
    alinv_apply_soil_context(
      soil_filter = ctx$soil_filter,
      include_soil_treatment = ctx$include_soil_treatment
    ) %>%
    dplyr::filter(.data$species %in% c("fagus", "quercus")) %>%
    dplyr::filter(!is.na(.data[[resp_var]])) %>%
    normalize_factors_tree() %>%
    dplyr::mutate(
      phase = factor(.data$phase, levels = c("until June", "July-August", "September+"))
    )

  if (!nrow(df_metric)) {
    return(df_metric)
  }

  if (isTRUE(within_phase)) {
    df_metric <- add_growth_phase_restart_rows(df_metric, resp_var)
  } else {
    df_metric <- df_metric %>%
      dplyr::mutate(is_phase_restart = FALSE)
  }

  df_metric
}

comparison_plot_summary <- function(df_metric, resp_var, within_phase = FALSE) {
  keys <- c("species", "robinia", "precipitation", "culture")
  if (isTRUE(within_phase)) {
    keys <- c(keys, "phase")
  }

  summarize_ts(df_metric, !!rlang::sym(resp_var), keys = keys) %>%
    dplyr::mutate(
      line_group = if (isTRUE(within_phase)) {
        interaction(.data$precipitation, .data$culture, .data$phase, drop = TRUE)
      } else {
        interaction(.data$precipitation, .data$culture, drop = TRUE)
      }
    )
}

plot_comparison_timeseries <- function(resp_var,
                                       title,
                                       y_label,
                                       within_phase = FALSE) {
  df_metric <- comparison_plot_data(resp_var = resp_var, within_phase = within_phase)

  if (!nrow(df_metric)) {
    return(alinv_empty_plot(paste("No data for", resp_var)))
  }

  df_sum <- comparison_plot_summary(
    df_metric = df_metric,
    resp_var = resp_var,
    within_phase = within_phase
  )

  p <- ggplot2::ggplot(
    df_sum,
    ggplot2::aes(
      x = .data$date,
      y = .data$mean,
      color = .data$precipitation,
      fill = .data$precipitation,
      linetype = .data$culture,
      group = .data$line_group
    )
  ) +
    ci_layers("band") +
    facet_tree() +
    scale_precip_color() +
    scale_precip_fill() +
    scale_culture_linetype() +
    ggplot2::scale_x_date(date_labels = "%b") +
    ggplot2::labs(
      x = "Date",
      y = y_label,
      title = title
    ) +
    theme_common()

  add_drought_bars_background(p, df_sum, show = FALSE)
}

collect_comparison_sem_matrix_data <- function(resp_vars, species_name) {
  out <- vector("list", length(resp_vars))

  for (i in seq_along(resp_vars)) {
    resp_var <- resp_vars[[i]]
    sem_res <- tryCatch(
      run_sem_for_trait(
        type = "tree",
        data_name = "growth",
        resp_var = resp_var,
        species = species_name,
        soil_type = ctx$soil_filter,
        include_soil_treatment = ctx$include_soil_treatment,
        include_interaction = FALSE,
        phase_window = "all",
        scale_all_numeric = TRUE,
        swc_source = "measured",
        force_run = FALSE
      ),
      error = function(e) {
        warning("Skipping SEM for ", species_name, " / ", resp_var, ": ", conditionMessage(e), call. = FALSE)
        NULL
      }
    )

    if (is.null(sem_res)) {
      next
    }

    matrix_df <- sem_res$matrix_data %||% tibble::tibble()
    if (!nrow(matrix_df) && !is.null(sem_res$effects) && nrow(sem_res$effects)) {
      matrix_df <- build_sem_matrix_data(
        effects_main = sem_res$effects,
        effects_int = sem_res$effects_int,
        resp_var = resp_var,
        species = species_name,
        include_interaction = FALSE,
        swc_source = "measured",
        phase = "all"
      )
    }

    if (!nrow(matrix_df)) {
      next
    }

    out[[i]] <- matrix_df
  }

  dplyr::bind_rows(out)
}

build_comparison_heatmap_page <- function(section_spec) {
  resp_order <- section_spec$heatmap_order
  resp_labels <- comparison_resp_label_map[resp_order]

  matrix_df <- purrr::map_dfr(
    c("fagus", "quercus"),
    ~ collect_comparison_sem_matrix_data(resp_order, species_name = .x)
  )

  if (!nrow(matrix_df)) {
    return(alinv_empty_plot(paste("No SEM matrix data for", section_spec$title)))
  }

  shared_limit <- compute_shared_sem_heatmap_limit(
    matrix_df %>% dplyr::filter(.data$effect_class == "main", .data$path_type == "total")
  )

  species_plots <- purrr::map(
    c("fagus", "quercus"),
    function(species_name) {
      panel_df <- prepare_sem_heatmap_panel(
        matrix_df = matrix_df %>%
          dplyr::filter(.data$species == species_name, .data$effect_class == "main"),
        path_type = "total",
        resp_labels = resp_labels,
        treat_labels = sem_heatmap_treatment_labels(),
        resp_var_order = resp_order
      )

      plot_sem_heatmap_panel(
        panel_df = panel_df,
        limit = shared_limit,
        title = stringr::str_to_title(species_name),
        annotate_values = TRUE,
        value_digits = 2
      )
    }
  )

  patchwork::wrap_plots(species_plots, ncol = 1) +
    patchwork::plot_annotation(
      title = section_spec$title,
      subtitle = paste(section_spec$subtitle, "| Total SEM effects")
    )
}

build_comparison_timeseries_page <- function(section_spec) {
  plots <- vector("list", 0)

  for (i in seq_len(nrow(comparison_metric_specs))) {
    spec_i <- comparison_metric_specs[i, ]

    for (j in seq_along(section_spec$plot_cols)) {
      resp_var <- spec_i[[section_spec$plot_cols[[j]]]][[1]]
      y_label <- spec_i[[section_spec$y_cols[[j]]]][[1]]
      within_phase <- identical(resp_var, spec_i$phase_inc_resp_var[[1]])
      title <- dplyr::case_when(
        identical(section_spec$plot_cols[[j]], "abs_resp_var") ~ spec_i$abs_panel_title[[1]],
        identical(section_spec$plot_cols[[j]], "phase_inc_resp_var") ~ spec_i$phase_panel_title[[1]],
        identical(section_spec$plot_cols[[j]], "t0_inc_resp_var") ~ spec_i$t0_panel_title[[1]],
        TRUE ~ paste(spec_i$metric_group[[1]], "-", spec_i[[section_spec$label_cols[[j]]]][[1]])
      )

      plots[[length(plots) + 1L]] <- plot_comparison_timeseries(
        resp_var = resp_var,
        title = title,
        y_label = y_label,
        within_phase = within_phase
      )
    }
  }

  ncol_layout <- length(section_spec$plot_cols)
  patchwork::wrap_plots(plots, ncol = ncol_layout, guides = "collect") +
    patchwork::plot_annotation(
      title = section_spec$title,
      subtitle = paste(section_spec$subtitle, "| Timeseries (mean ± SE)")
    ) &
    ggplot2::theme(legend.position = "bottom")
}

pdf_path <- file.path(
  alinv_analysis_path("figures", create_dir = TRUE),
  "volume-calculation-comparison.pdf"
)

grDevices::cairo_pdf(pdf_path, width = 18, height = 12, onefile = TRUE)
for (section_name in names(comparison_sections)) {
  section_spec <- comparison_sections[[section_name]]
  print(build_comparison_heatmap_page(section_spec))
  print(build_comparison_timeseries_page(section_spec))
}
grDevices::dev.off()

message("Saved volume calculation comparison report to: ", pdf_path)
