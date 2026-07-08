season_background_layer <- function(alpha = 0.45) {
  ggplot2::geom_rect(
    data = SEASON_WINDOWS,
    ggplot2::aes(xmin = .data$start, xmax = .data$end, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    fill = SEASON_WINDOWS$fill,
    alpha = alpha,
    color = NA
  )
}

FIG2_X_LIMITS <- c(as.Date("2025-03-01"), as.Date("2025-12-31"))

format_fig2_axis <- function(p, show_x = FALSE) {
  p +
    ggplot2::scale_x_date(
      date_breaks = "2 months",
      date_labels = "%b"
    ) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = if (isTRUE(show_x)) ggplot2::element_text() else ggplot2::element_blank(),
      axis.ticks.x = if (isTRUE(show_x)) ggplot2::element_line(linewidth = 0.25) else ggplot2::element_blank(),
      plot.title.position = "plot",
      plot.margin = ggplot2::margin(1.5, 4, 1.5, 4)
    )
}

summarize_fig2_tree <- function(data_name, value_col) {
  df <- get_data("tree", data_name) %>%
    dplyr::filter(.data$species %in% c("fagus", "quercus"), !is.na(.data[[value_col]])) %>%
    dplyr::mutate(
      species = factor(.data$species, levels = c("fagus", "quercus")),
      robinia = factor(.data$robinia, levels = names(ROBINIA_COLORS)),
      precipitation = factor(.data$precipitation, levels = c("control", "drought")),
      culture = factor(.data$culture, levels = c("mono", "mixed"))
    )

  thin <- df %>%
    summary_se(c("date", "species", "robinia", "precipitation", "culture"), value_col) %>%
    dplyr::mutate(line_group = interaction(.data$robinia, .data$precipitation, .data$culture, drop = TRUE))

  plot_level <- df %>%
    dplyr::group_by(.data$date, .data$species, .data$robinia, .data$boxlabel) %>%
    dplyr::summarise(value = mean(.data[[value_col]], na.rm = TRUE), .groups = "drop")

  thick <- plot_level %>%
    summary_se(c("date", "species", "robinia"), "value")

  list(thin = thin, thick = thick)
}

summarize_fig2_swc <- function() {
  df <- get_data("box", "soilwater") %>%
    dplyr::filter(!is.na(.data$swc), lubridate::year(.data$date) == 2025) %>%
    dplyr::mutate(
      robinia = factor(.data$robinia, levels = names(ROBINIA_COLORS)),
      precipitation = factor(.data$precipitation, levels = c("control", "drought")),
      culture = factor(.data$culture, levels = c("mono", "mixed"))
    )

  thin <- df %>%
    summary_se(c("date", "robinia", "precipitation", "culture"), "swc") %>%
    dplyr::mutate(line_group = interaction(.data$robinia, .data$precipitation, .data$culture, drop = TRUE))

  plot_level <- df %>%
    dplyr::group_by(.data$date, .data$robinia, .data$boxlabel) %>%
    dplyr::summarise(value = mean(.data$swc, na.rm = TRUE), .groups = "drop")

  thick <- plot_level %>%
    summary_se(c("date", "robinia"), "value")

  list(thin = thin, thick = thick)
}

build_fig2_response_panel <- function(summary_obj, species_key = NULL, y_lab, title, y_limits = NULL, show_x = FALSE) {
  thin <- summary_obj$thin
  thick <- summary_obj$thick

  if (!is.null(species_key)) {
    thin <- dplyr::filter(thin, .data$species == .env$species_key)
    thick <- dplyr::filter(thick, .data$species == .env$species_key)
  }

  y_min <- if (!is.null(y_limits)) y_limits[[1]] else min(thick$mean - thick$se, na.rm = TRUE)
  drought_y <- y_min

  p <- ggplot2::ggplot() +
    season_background_layer() +
    ggplot2::geom_line(
      data = thin,
      ggplot2::aes(
        x = .data$date,
        y = .data$mean,
        color = .data$robinia,
        linetype = .data$precipitation,
        alpha = .data$culture,
        group = .data$line_group
      ),
      linewidth = 0.3
    ) +
    ggplot2::geom_errorbar(
      data = thick,
      ggplot2::aes(x = .data$date, ymin = .data$mean - .data$se, ymax = .data$mean + .data$se, color = .data$robinia),
      width = 2.5,
      linewidth = 0.28,
      alpha = 0.9
    ) +
    ggplot2::geom_line(
      data = thick,
      ggplot2::aes(x = .data$date, y = .data$mean, color = .data$robinia, group = .data$robinia),
      linewidth = 0.75
    ) +
    ggplot2::geom_point(
      data = thick,
      ggplot2::aes(x = .data$date, y = .data$mean, color = .data$robinia),
      size = 0.8
    ) +
    ggplot2::scale_color_manual(values = ROBINIA_COLORS, labels = ROBINIA_LABELS, name = "Robinia", drop = FALSE) +
    ggplot2::scale_linetype_manual(
      values = c(control = "solid", drought = "42"),
      labels = c(control = "control", drought = "reduced"),
      name = "Precipitation",
      drop = FALSE
    ) +
    ggplot2::scale_alpha_manual(
      values = c(mono = 0.35, mixed = 0.75),
      labels = c(mono = "mono", mixed = "mixed"),
      name = "Culture",
      drop = FALSE
    ) +
    ggplot2::labs(x = NULL, y = y_lab, title = title) +
    theme_alinv_pub(base_size = 6.6)

  if (!is.null(y_limits)) {
    p <- p + ggplot2::coord_cartesian(xlim = FIG2_X_LIMITS, ylim = y_limits)
  } else {
    p <- p + ggplot2::coord_cartesian(xlim = FIG2_X_LIMITS)
  }

  p <- add_drought_segments(p, drought_y)
  format_fig2_axis(p, show_x = show_x)
}

read_daily_temperature <- function() {
  meteo_file <- file.path(ALINV_PROJECT_ROOT, "data", "raw", "sensor_data", "meteo_10min.dat")
  readr::read_csv(
    meteo_file,
    skip = 4,
    col_names = c("timestamp", "record", "ta", "rh", "ap", "vpd", "globrad"),
    show_col_types = FALSE
  ) %>%
    dplyr::mutate(
      timestamp = lubridate::ymd_hms(.data$timestamp, quiet = TRUE),
      date = as.Date(.data$timestamp)
    ) %>%
    dplyr::filter(lubridate::year(.data$date) == 2025) %>%
    dplyr::group_by(.data$date) %>%
    dplyr::summarise(ta = mean(.data$ta, na.rm = TRUE), .groups = "drop")
}

build_fig2_weather_panel <- function() {
  ta <- read_daily_temperature()
  precip <- readr::read_csv(file.path(ALINV_PROJECT_ROOT, "data", "interim", "site_precipitation_daily.csv"), show_col_types = FALSE) %>%
    dplyr::filter(lubridate::year(.data$date) == 2025) %>%
    dplyr::select(date, precip_mm)

  weather <- dplyr::full_join(ta, precip, by = "date") %>%
    dplyr::arrange(.data$date)

  y_min <- -10
  y_max <- 32
  precip_max <- max(weather$precip_mm, na.rm = TRUE)
  if (!is.finite(precip_max) || precip_max <= 0) precip_max <- 1

  weather <- weather %>%
    dplyr::mutate(precip_scaled = y_min + (.data$precip_mm / precip_max) * (y_max - y_min))

  p <- ggplot2::ggplot() +
    season_background_layer() +
    ggplot2::geom_rect(
      data = dplyr::filter(weather, !is.na(.data$precip_mm), .data$precip_mm > 0),
      ggplot2::aes(
        xmin = .data$date - 0.45,
        xmax = .data$date + 0.45,
        ymin = y_min,
        ymax = .data$precip_scaled
      ),
      inherit.aes = FALSE,
      fill = "#7DBAD6",
      alpha = 0.75,
      color = NA
    ) +
    ggplot2::geom_line(
      data = weather,
      ggplot2::aes(x = .data$date, y = .data$ta),
      color = "#E64B35",
      linewidth = 0.35
    ) +
    ggplot2::scale_y_continuous(
      name = "Ta (deg C)",
      sec.axis = ggplot2::sec_axis(
        transform = ~ (. - y_min) / (y_max - y_min) * precip_max,
        name = "Precipitation (mm)"
      )
    ) +
    ggplot2::labs(x = NULL, title = "Site climate") +
    theme_alinv_pub(base_size = 6.6) +
    ggplot2::theme(
      legend.position = "none",
      plot.margin = ggplot2::margin(1.5, 4, 1.5, 4)
    )

  p <- add_drought_segments(p, y_min)
  p <- p + ggplot2::coord_cartesian(xlim = FIG2_X_LIMITS, ylim = c(y_min, y_max))
  format_fig2_axis(p, show_x = TRUE)
}

make_fig2 <- function() {
  growth <- summarize_fig2_tree("growth", "diameter_inc_t0")
  qy <- summarize_fig2_tree("quantum_yield", "qy")
  swc <- summarize_fig2_swc()

  p_a <- build_fig2_response_panel(growth, "fagus", "Diameter incr. (mm)", "Fagus", y_limits = c(0, 8))
  p_b <- build_fig2_response_panel(growth, "quercus", "Diameter incr. (mm)", "Quercus", y_limits = c(0, 8))
  p_c <- build_fig2_response_panel(qy, "fagus", "Fv/Fm", "Fagus", y_limits = c(0, 1.02))
  p_d <- build_fig2_response_panel(qy, "quercus", "Fv/Fm", "Quercus", y_limits = c(0, 1.02))
  p_e <- build_fig2_response_panel(swc, NULL, "SWC (%)", "Soil water content", y_limits = c(0, 31))
  p_f <- build_fig2_weather_panel()

  fig <- (p_a / p_b / p_c / p_d / p_e / p_f) +
    patchwork::plot_layout(heights = c(1, 1, 1, 1, 1, 1), guides = "collect") +
    patchwork::plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
    ggplot2::theme(
      legend.position = "bottom",
      legend.title = ggplot2::element_text(size = 6.5),
      legend.text = ggplot2::element_text(size = 6.2)
    )

  alinv_save_pdf(fig, "fig2_variation_timeseries.pdf", height_mm = 220)
}
