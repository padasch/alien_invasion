growth_metric_catalog <- function() {
  tibble::tribble(
    ~metric_group, ~variant_key, ~resp_var, ~variant_label, ~y_label, ~panel_title, ~within_phase,
    "Height", "absolute", "height", "Absolute size", "Height (cm)", "Height over time", FALSE,
    "Height", "inc_t0_abs", "height_inc_t0", "Absolute increment from first measurement", "Height increment from first measurement (cm)", "Height increment from first measurement", FALSE,
    "Height", "inc_t0_rel", "height_inc_t0_rel", "Relative increment from first measurement", "Relative height increment from first measurement (-)", "Height relative increment to first measurement", FALSE,
    "Height", "inc_phase_abs", "height_inc_phase_abs", "Absolute increment within phase", "Height increment within phase (cm)", "Height increment within phase", TRUE,
    "Height", "inc_phase_rel", "height_inc_phase_rel", "Relative increment within phase", "Relative height increment within phase (-)", "Height relative increment within phase", TRUE,
    "Diameter", "absolute", "diameter", "Absolute size", "Diameter (mm)", "Diameter over time", FALSE,
    "Diameter", "inc_t0_abs", "diameter_inc_t0", "Absolute increment from first measurement", "Diameter increment from first measurement (mm)", "Diameter increment from first measurement", FALSE,
    "Diameter", "inc_t0_rel", "diameter_inc_t0_rel", "Relative increment from first measurement", "Relative diameter increment from first measurement (-)", "Diameter relative increment to first measurement", FALSE,
    "Diameter", "inc_phase_abs", "diameter_inc_phase_abs", "Absolute increment within phase", "Diameter increment within phase (mm)", "Diameter increment within phase", TRUE,
    "Diameter", "inc_phase_rel", "diameter_inc_phase_rel", "Relative increment within phase", "Relative diameter increment within phase (-)", "Diameter relative increment within phase", TRUE,
    "Volume", "absolute", "volume", "Absolute size", "Volume (allometric proxy, g)", "Volume over time", FALSE,
    "Volume", "inc_t0_abs", "volume_inc_t0", "Absolute increment from first measurement", "Volume increment from first measurement (allometric proxy, g)", "Volume increment from first measurement", FALSE,
    "Volume", "inc_t0_rel", "volume_inc_t0_rel", "Relative increment from first measurement", "Relative volume increment from first measurement (-)", "Volume relative increment to first measurement", FALSE,
    "Volume", "inc_phase_abs", "volume_inc_phase_abs", "Absolute increment within phase", "Volume increment within phase (allometric proxy, g)", "Volume increment within phase", TRUE,
    "Volume", "inc_phase_rel", "volume_inc_phase_rel", "Relative increment within phase", "Relative volume increment within phase (-)", "Volume relative increment within phase", TRUE
  )
}

growth_metric_spec <- function(resp_var) {
  growth_metric_catalog() %>%
    dplyr::filter(.data$resp_var == .env$resp_var) %>%
    dplyr::slice_head(n = 1)
}

growth_phase_reset_metrics <- function() {
  growth_metric_catalog() %>%
    dplyr::filter(.data$within_phase) %>%
    dplyr::pull(.data$resp_var)
}

axis_limits_for_growth_metric <- function(resp_var) {
  if (identical(resp_var, "height")) {
    return(get_ylim("height"))
  }
  if (identical(resp_var, "diameter")) {
    return(get_ylim("diameter"))
  }
  NULL
}

add_growth_phase_restart_rows <- function(df_metric, resp_var) {
  phase_levels <- levels(df_metric$phase)
  if (is.null(phase_levels) || !"phase" %in% names(df_metric)) {
    return(df_metric %>% dplyr::mutate(is_phase_restart = FALSE))
  }

  df_base <- df_metric %>%
    dplyr::mutate(is_phase_restart = FALSE)

  restart_rows <- df_base %>%
    dplyr::group_by(.data$tree_id) %>%
    dplyr::group_modify(function(.x, .y) {
      purrr::map_dfr(phase_levels[-1], function(phase_i) {
        current_rows <- .x %>%
          dplyr::filter(as.character(.data$phase) == phase_i) %>%
          dplyr::arrange(.data$date)
        previous_rows <- .x %>%
          dplyr::filter(as.integer(.data$phase) < match(phase_i, phase_levels)) %>%
          dplyr::arrange(.data$date)

        if (!nrow(current_rows) || !nrow(previous_rows)) {
          return(tibble::tibble())
        }

        restart_row <- current_rows[1, , drop = FALSE]
        restart_row$date <- max(previous_rows$date, na.rm = TRUE)
        restart_row[[resp_var]] <- 0
        restart_row$is_phase_restart <- TRUE
        restart_row
      })
    }) %>%
    dplyr::ungroup()

  dplyr::bind_rows(df_base, restart_rows) %>%
    dplyr::arrange(.data$tree_id, .data$date, .data$is_phase_restart)
}

prepare_growth_metric_plot_data <- function(resp_var,
                                            soil_type = NULL,
                                            include_soil_treatment = NULL) {
  ctx <- alinv_get_analysis_context()
  soil_type <- soil_type %||% ctx$soil_filter %||% "both"
  include_soil_treatment <- alinv_resolve_include_soil_treatment(
    include_soil_treatment = include_soil_treatment,
    soil_filter = soil_type
  )

  spec <- growth_metric_spec(resp_var)
  within_phase <- isTRUE(spec$within_phase[[1]])

  df_metric <- get_data("tree", "growth") %>%
    alinv_apply_soil_context(
      soil_filter = soil_type,
      include_soil_treatment = include_soil_treatment
    ) %>%
    dplyr::filter(!is.na(.data[[resp_var]])) %>%
    normalize_factors_tree() %>%
    dplyr::mutate(
      phase = factor(.data$phase, levels = c("until June", "July-August", "September+"))
    )

  if (!nrow(df_metric)) {
    return(df_metric)
  }

  if (within_phase) {
    df_metric <- add_growth_phase_restart_rows(df_metric, resp_var)
  } else {
    df_metric <- df_metric %>%
      dplyr::mutate(is_phase_restart = FALSE)
  }

  df_metric
}

summarize_growth_metric_plot_data <- function(df_metric, resp_var, within_phase = FALSE) {
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

plot_growth_metric_ts <- function(resp_var,
                                  y_label = NULL,
                                  panel_title = NULL,
                                  soil_type = NULL,
                                  include_soil_treatment = NULL,
                                  drought_bars = TRUE,
                                  style = "band") {
  spec <- growth_metric_spec(resp_var)
  if (!nrow(spec)) {
    stop("Unknown growth response variable: ", resp_var, call. = FALSE)
  }

  y_label <- y_label %||% spec$y_label[[1]]
  panel_title <- panel_title %||% spec$panel_title[[1]]
  within_phase <- isTRUE(spec$within_phase[[1]])

  df_metric <- prepare_growth_metric_plot_data(
    resp_var = resp_var,
    soil_type = soil_type,
    include_soil_treatment = include_soil_treatment
  )

  if (!nrow(df_metric)) {
    return(alinv_empty_plot(paste("No data available for", resp_var)))
  }

  df_sum <- summarize_growth_metric_plot_data(
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
    ci_layers(style) +
    facet_tree() +
    scale_precip_color() +
    scale_precip_fill() +
    scale_culture_linetype() +
    ggplot2::scale_x_date(date_labels = "%m/%y") +
    ggplot2::labs(
      x = "Date",
      y = y_label,
      title = panel_title
    ) +
    theme_common()

  y_limits <- axis_limits_for_growth_metric(resp_var)
  if (!is.null(y_limits)) {
    p <- p + ggplot2::coord_cartesian(ylim = y_limits)
  }

  add_drought_bars_background(p, df_sum, show = drought_bars)
}

build_growth_plot_bundle <- function(soil_type = NULL,
                                     include_soil_treatment = NULL) {
  catalog <- growth_metric_catalog()

  raw_list <- vector("list", nrow(catalog))
  summary_list <- vector("list", nrow(catalog))

  for (i in seq_len(nrow(catalog))) {
    spec_i <- catalog[i, ]
    df_metric <- prepare_growth_metric_plot_data(
      resp_var = spec_i$resp_var[[1]],
      soil_type = soil_type,
      include_soil_treatment = include_soil_treatment
    )

    if (!nrow(df_metric)) {
      next
    }

    raw_list[[i]] <- df_metric %>%
      dplyr::transmute(
        tree_id, boxlabel, date, phase, species, robinia, precipitation, culture, soiltype,
        metric_group = spec_i$metric_group[[1]],
        variant_key = spec_i$variant_key[[1]],
        resp_var = spec_i$resp_var[[1]],
        variant_label = spec_i$variant_label[[1]],
        y_label = spec_i$y_label[[1]],
        value = .data[[spec_i$resp_var[[1]]]],
        within_phase = spec_i$within_phase[[1]],
        is_phase_restart
      )

    summary_list[[i]] <- summarize_growth_metric_plot_data(
      df_metric = df_metric,
      resp_var = spec_i$resp_var[[1]],
      within_phase = spec_i$within_phase[[1]]
    ) %>%
      dplyr::mutate(
        metric_group = spec_i$metric_group[[1]],
        variant_key = spec_i$variant_key[[1]],
        resp_var = spec_i$resp_var[[1]],
        variant_label = spec_i$variant_label[[1]],
        y_label = spec_i$y_label[[1]],
        within_phase = spec_i$within_phase[[1]]
      )
  }

  list(
    raw = dplyr::bind_rows(raw_list),
    summary = dplyr::bind_rows(summary_list)
  )
}

export_growth_plot_bundle <- function(bundle,
                                      export_dir = alinv_data_path("size_trajectories", create_dir = TRUE)) {
  dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(bundle$raw %||% tibble::tibble(), file.path(export_dir, "size-trajectories-raw.csv"))
  readr::write_csv(bundle$summary %||% tibble::tibble(), file.path(export_dir, "size-trajectories-summary.csv"))
  invisible(bundle)
}

growth_metric_slug <- function(metric_group, variant_key) {
  metric_slug <- gsub("[^a-z0-9]+", "-", tolower(metric_group))
  variant_slug <- gsub("[^a-z0-9]+", "-", tolower(variant_key))
  paste(metric_slug, variant_slug, sep = "-")
}

save_growth_plot_png <- function(plot_obj,
                                 file_path,
                                 width = 9.5,
                                 height = 5.8,
                                 dpi = 300) {
  dir.create(dirname(file_path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    filename = file_path,
    plot = plot_obj,
    width = width,
    height = height,
    dpi = dpi
  )
  file_path
}

growth_markdown_image <- function(path,
                                  width = "100%") {
  if (is.null(path) || is.na(path) || !nzchar(path)) {
    return("_Figure unavailable._\n\n")
  }

  if (!file.exists(path)) {
    return("_Figure unavailable._\n\n")
  }

  paste0(
    "<img src=\"",
    knitr::image_uri(path),
    "\" style=\"width: ",
    width,
    ";\" />\n\n"
  )
}

build_growth_figure_manifest <- function(soil_type = NULL,
                                         include_soil_treatment = NULL,
                                         species_vec = c("fagus", "quercus"),
                                         export_dir = alinv_data_path("size_trajectories", "figures", create_dir = TRUE)) {
  catalog <- growth_metric_catalog()
  manifest <- vector("list", nrow(catalog) * (1 + 2 * length(species_vec)))
  idx <- 0L

  dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)

  for (i in seq_len(nrow(catalog))) {
    spec_i <- catalog[i, ]
    variant_slug <- growth_metric_slug(spec_i$metric_group[[1]], spec_i$variant_key[[1]])

    ts_plot <- plot_growth_metric_ts(
      resp_var = spec_i$resp_var[[1]],
      y_label = spec_i$y_label[[1]],
      panel_title = spec_i$panel_title[[1]],
      soil_type = soil_type,
      include_soil_treatment = include_soil_treatment,
      drought_bars = TRUE
    )

    ts_path <- file.path(export_dir, "timeseries", paste0(variant_slug, ".png"))
    save_growth_plot_png(ts_plot, ts_path)

    idx <- idx + 1L
    manifest[[idx]] <- tibble::tibble(
      figure_type = "timeseries",
      metric_group = spec_i$metric_group[[1]],
      variant_key = spec_i$variant_key[[1]],
      resp_var = spec_i$resp_var[[1]],
      variant_label = spec_i$variant_label[[1]],
      species = NA_character_,
      path = ts_path,
      path_rel = alinv_project_relative_path(ts_path)
    )

    for (species_i in species_vec) {
      res_i <- make_temporal_sem_combo(
        type = "tree",
        data_name = "growth",
        resp_var = spec_i$resp_var[[1]],
        species = species_i,
        soil_type = soil_type,
        include_soil_treatment = include_soil_treatment,
        include_interaction = FALSE
      )

      temporal_path <- file.path(export_dir, "temporal", paste0(species_i, "-", variant_slug, ".png"))
      sem_path <- file.path(export_dir, "sem", paste0(species_i, "-", variant_slug, ".png"))

      save_growth_plot_png(res_i$temporal$plot, temporal_path)
      save_growth_plot_png(res_i$sem$plot, sem_path)

      idx <- idx + 1L
      manifest[[idx]] <- tibble::tibble(
        figure_type = "temporal",
        metric_group = spec_i$metric_group[[1]],
        variant_key = spec_i$variant_key[[1]],
        resp_var = spec_i$resp_var[[1]],
        variant_label = spec_i$variant_label[[1]],
        species = species_i,
        path = temporal_path,
        path_rel = alinv_project_relative_path(temporal_path)
      )

      idx <- idx + 1L
      manifest[[idx]] <- tibble::tibble(
        figure_type = "sem",
        metric_group = spec_i$metric_group[[1]],
        variant_key = spec_i$variant_key[[1]],
        resp_var = spec_i$resp_var[[1]],
        variant_label = spec_i$variant_label[[1]],
        species = species_i,
        path = sem_path,
        path_rel = alinv_project_relative_path(sem_path)
      )
    }
  }

  dplyr::bind_rows(manifest) %>%
    dplyr::filter(!is.na(.data$figure_type))
}

get_growth_figure_path <- function(manifest,
                                   figure_type,
                                   resp_var,
                                   species = NA_character_) {
  df <- manifest %>%
    dplyr::filter(
      .data$figure_type == .env$figure_type,
      .data$resp_var == .env$resp_var
    )

  if (is.na(species)) {
    df <- df %>% dplyr::filter(is.na(.data$species))
  } else {
    df <- df %>% dplyr::filter(.data$species == .env$species)
  }

  if (!nrow(df)) {
    return(NA_character_)
  }

  df %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::pull(.data$path) %>%
    .[[1]]
}

combined_volume_sem_response_vars <- function() {
  c(
    "chl",
    "condition",
    "volume_inc_phase_abs",
    "stage",
    "remaining_green",
    "chlavg",
    "qy"
  )
}

combined_volume_sem_resp_labels <- function(resp_labels = alinv_response_labels()) {
  keep <- combined_volume_sem_response_vars()
  out <- resp_labels[keep]
  out[!is.na(out)]
}

compute_species_timeseries_limits <- function(df_sum,
                                              species_col = "species",
                                              lower_col = "mean",
                                              se_col = "se",
                                              pad = 0.05) {
  if (!nrow(df_sum)) {
    return(tibble::tibble())
  }

  df_sum %>%
    dplyr::mutate(
      lower_val = .data[[lower_col]] - .data[[se_col]],
      upper_val = .data[[lower_col]] + .data[[se_col]]
    ) %>%
    dplyr::group_by(.data[[species_col]]) %>%
    dplyr::summarise(
      y_lower = min(lower_val, 0, na.rm = TRUE),
      y_upper = max(upper_val, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    {
      out <- .
      names(out)[1] <- "species_value"
      out
    } %>%
    dplyr::mutate(
      y_pad = pmax((.data$y_upper - .data$y_lower) * pad, 0),
      y_lower = .data$y_lower,
      y_upper = .data$y_upper + .data$y_pad
    ) %>%
    dplyr::select(-.data$y_pad)
}

plot_volume_phase_species_panel <- function(species_name,
                                            robinia_level,
                                            soil_type = NULL,
                                            include_soil_treatment = NULL,
                                            y_limits = NULL,
                                            drought_bars = FALSE,
                                            show_y_axis = TRUE,
                                            show_legend = TRUE,
                                            title_prefix = NULL) {
  df_metric <- prepare_growth_metric_plot_data(
    resp_var = "volume_inc_phase_abs",
    soil_type = soil_type,
    include_soil_treatment = include_soil_treatment
  ) %>%
    dplyr::filter(
      as.character(.data$species) == .env$species_name,
      as.character(.data$robinia) == .env$robinia_level
    )

  if (!nrow(df_metric)) {
    return(alinv_empty_plot(paste("No data available for", species_name, robinia_level)))
  }

  df_sum <- summarize_growth_metric_plot_data(
    df_metric = df_metric,
    resp_var = "volume_inc_phase_abs",
    within_phase = TRUE
  )

  title_text <- paste(
    title_prefix %||% stringr::str_to_title(species_name),
    alinv_level_labels("robinia")[[robinia_level]],
    sep = ": "
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
    scale_precip_color() +
    scale_precip_fill() +
    scale_culture_linetype() +
    ggplot2::scale_x_date(date_labels = "%b") +
    ggplot2::labs(
      x = "Date",
      y = "Volume increment within phase (allometric proxy, g)",
      title = title_text
    ) +
    theme_common()

  p <- add_drought_bars_background(p, df_sum, show = drought_bars)

  if (!is.null(y_limits) && length(y_limits) == 2L) {
    p <- p + ggplot2::coord_cartesian(ylim = y_limits)
  }

  if (!isTRUE(show_y_axis)) {
    p <- p + ggplot2::theme(
      axis.title.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank()
    )
  }

  if (!isTRUE(show_legend)) {
    p <- p + ggplot2::theme(legend.position = "none")
  }

  p
}

plot_volume_sem_species_row <- function(species_name,
                                        matrix_df,
                                        heatmap_limit,
                                        soil_type = NULL,
                                        include_soil_treatment = NULL) {
  df_metric <- prepare_growth_metric_plot_data(
    resp_var = "volume_inc_phase_abs",
    soil_type = soil_type,
    include_soil_treatment = include_soil_treatment
  ) %>%
    dplyr::filter(as.character(.data$species) == .env$species_name)

  if (!nrow(df_metric)) {
    return(alinv_empty_plot(paste("No volume phase data for", species_name)))
  }

  df_sum <- summarize_growth_metric_plot_data(
    df_metric = df_metric,
    resp_var = "volume_inc_phase_abs",
    within_phase = TRUE
  )

  y_limit_tbl <- compute_species_timeseries_limits(df_sum)
  y_limits <- y_limit_tbl %>%
    dplyr::filter(as.character(.data$species_value) == .env$species_name) %>%
    dplyr::select(.data$y_lower, .data$y_upper)

  y_limits <- if (nrow(y_limits)) as.numeric(y_limits[1, ]) else NULL

  p_without <- plot_volume_phase_species_panel(
    species_name = species_name,
    robinia_level = "without-robinia",
    soil_type = soil_type,
    include_soil_treatment = include_soil_treatment,
    y_limits = y_limits,
    drought_bars = FALSE,
    show_y_axis = TRUE,
    show_legend = TRUE,
    title_prefix = stringr::str_to_title(species_name)
  )

  p_with <- plot_volume_phase_species_panel(
    species_name = species_name,
    robinia_level = "with-robinia",
    soil_type = soil_type,
    include_soil_treatment = include_soil_treatment,
    y_limits = y_limits,
    drought_bars = FALSE,
    show_y_axis = FALSE,
    show_legend = FALSE,
    title_prefix = NULL
  )

  panel_df <- prepare_sem_heatmap_panel(
    matrix_df = matrix_df %>%
      dplyr::filter(
        .data$species == .env$species_name,
        .data$effect_class == "main",
        .data$response_var %in% combined_volume_sem_response_vars()
      ),
    path_type = "total",
    resp_labels = combined_volume_sem_resp_labels(),
    treat_labels = sem_heatmap_treatment_labels(),
    resp_var_order = combined_volume_sem_response_vars()
  )

  p_heatmap <- plot_sem_heatmap_panel(
    panel_df = panel_df,
    limit = heatmap_limit,
    title = "Total SEM effects",
    annotate_values = TRUE,
    value_digits = 2,
    value_text_size = 2.7
  )

  (p_without | p_with | p_heatmap) +
    patchwork::plot_layout(widths = c(1, 1, 1.15), guides = "collect") &
    ggplot2::theme(legend.position = "bottom")
}

plot_volume_sem_summary_figure <- function(matrix_df,
                                           soil_type = NULL,
                                           include_soil_treatment = NULL,
                                           species_vec = c("fagus", "quercus")) {
  if (is.null(matrix_df) || !nrow(matrix_df)) {
    return(alinv_empty_plot("No SEM matrix data available for combined volume summary figure"))
  }

  heatmap_limit <- compute_shared_sem_heatmap_limit(matrix_df)

  row_plots <- purrr::map(
    species_vec,
    ~ plot_volume_sem_species_row(
      species_name = .x,
      matrix_df = matrix_df,
      heatmap_limit = heatmap_limit,
      soil_type = soil_type,
      include_soil_treatment = include_soil_treatment
    )
  )

  patchwork::wrap_plots(row_plots, ncol = 1, guides = "collect")
}
