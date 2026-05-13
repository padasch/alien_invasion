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
    "Volume", "absolute", "volume", "Absolute size", "Volume (cm^3)", "Volume over time", FALSE,
    "Volume", "inc_t0_abs", "volume_inc_t0", "Absolute increment from first measurement", "Volume increment from first measurement (cm^3)", "Volume increment from first measurement", FALSE,
    "Volume", "inc_t0_rel", "volume_inc_t0_rel", "Relative increment from first measurement", "Relative volume increment from first measurement (-)", "Volume relative increment to first measurement", FALSE,
    "Volume", "inc_phase_abs", "volume_inc_phase_abs", "Absolute increment within phase", "Volume increment within phase (cm^3)", "Volume increment within phase", TRUE,
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
