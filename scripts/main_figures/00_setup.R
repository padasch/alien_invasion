alinv_current_source_file <- function() {
  for (i in rev(seq_along(sys.frames()))) {
    file_i <- sys.frames()[[i]]$ofile %||% NULL
    if (!is.null(file_i) && nzchar(file_i)) {
      return(file_i)
    }
  }

  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) {
    return(gsub("~[+]~", " ", sub("^--file=", "", file_arg[[1]])))
  }

  stop("Could not resolve the current source file. Run from make_all_figures.R.", call. = FALSE)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

ALINV_FINAL_SCRIPT_DIR <- dirname(normalizePath(alinv_current_source_file(), winslash = "/", mustWork = TRUE))
ALINV_PROJECT_ROOT <- normalizePath(file.path(ALINV_FINAL_SCRIPT_DIR, "..", ".."), winslash = "/", mustWork = TRUE)
ALINV_FINAL_OUTPUT_DIR <- file.path(ALINV_PROJECT_ROOT, "output", "main_figures")
ALINV_FINAL_CACHE_DIR <- file.path(tempdir(), "alinv-final-figures-cache")
dir.create(ALINV_FINAL_OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(ALINV_FINAL_CACHE_DIR, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(XDG_CACHE_HOME = ALINV_FINAL_CACHE_DIR)

setwd(ALINV_PROJECT_ROOT)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(lubridate)
  library(patchwork)
  library(scales)
})

suppressMessages(suppressPackageStartupMessages({
  auxiliary_functions <- file.path(ALINV_PROJECT_ROOT, "scripts", "auxiliary", "functions")
  source(file.path(auxiliary_functions, "_source.R"))
  source(file.path(auxiliary_functions, "1-summary-figures.R"))
  source(file.path(auxiliary_functions, "8-size-trajectories.R"))
  source(file.path(auxiliary_functions, "11-bootstrap-inference.R"))
}))

alinv_set_analysis_context(
  scenario_label = "both soils (without soil as treatment)",
  soil_filter = "both",
  include_soil_treatment = FALSE,
  analysis_date = Sys.Date(),
  output_root = file.path(ALINV_PROJECT_ROOT, "output"),
  create_dirs = TRUE
)

FIG_WIDTH_MM <- 160
FIG_WIDTH_IN <- FIG_WIDTH_MM / 25.4

SPECIES_LABELS <- c(fagus = "Fagus", quercus = "Quercus")
ROBINIA_LABELS <- c(`without-robinia` = "without robinia", `with-robinia` = "with robinia")
ROBINIA_COLORS <- c(`without-robinia` = "#4F6674", `with-robinia` = "#FF5A45")
PRECIP_COLORS <- c(control = "#4F6674", drought = "#D65F5F")
EFFECT_COLORS <- c(
  culture = "#1B9E77",
  precipitation = "#D95F02",
  robinia = "#7570B3",
  extreme_event = "#666666"
)
EFFECT_LABELS_SHORT <- c(
  culture = "Culture (mono -> mixed)",
  precipitation = "Precipitation (control -> reduced)",
  robinia = "Robinia (without -> with)",
  extreme_event = "Extreme event (no -> yes)"
)

SUMMER_START <- as.Date("2025-06-20")
SUMMER_END <- as.Date("2025-09-01")
DROUGHT_WINDOWS <- tibble::tribble(
  ~start, ~end,
  as.Date("2025-06-20"), as.Date("2025-07-02"),
  as.Date("2025-08-12"), as.Date("2025-08-20")
)
SEASON_WINDOWS <- tibble::tribble(
  ~season, ~start, ~end, ~fill,
  "Spring", as.Date("2025-03-01"), SUMMER_START, "#C9F7C5",
  "Summer", SUMMER_START, SUMMER_END, "#A9C9AA",
  "Autumn", SUMMER_END, as.Date("2025-12-15"), "#F8F6D6"
)

theme_alinv_pub <- function(base_size = 7) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      text = ggplot2::element_text(family = "Helvetica", color = "black"),
      axis.text = ggplot2::element_text(color = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.25, color = "black"),
      axis.line = ggplot2::element_line(linewidth = 0.25, color = "black"),
      panel.border = ggplot2::element_rect(fill = NA, color = "grey55", linewidth = 0.25),
      panel.grid.major.y = ggplot2::element_line(color = "grey90", linewidth = 0.18),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.key.height = grid::unit(3.5, "mm"),
      legend.key.width = grid::unit(7, "mm"),
      legend.box.spacing = grid::unit(1, "mm"),
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 1, margin = ggplot2::margin(b = 2)),
      plot.tag = ggplot2::element_text(face = "bold", size = base_size + 1),
      strip.background = ggplot2::element_rect(fill = "grey12", color = "grey12", linewidth = 0.25),
      strip.text = ggplot2::element_text(color = "white", face = "bold", size = base_size)
    )
}

alinv_save_pdf <- function(plot, filename, height_mm, width_mm = FIG_WIDTH_MM) {
  if (width_mm > FIG_WIDTH_MM) {
    stop("Requested figure width exceeds 160 mm: ", width_mm, " mm", call. = FALSE)
  }

  out_path <- file.path(ALINV_FINAL_OUTPUT_DIR, filename)
  ggplot2::ggsave(
    filename = out_path,
    plot = plot,
    width = width_mm / 25.4,
    height = height_mm / 25.4,
    units = "in",
    device = grDevices::cairo_pdf
  )
  normalizePath(out_path, winslash = "/", mustWork = TRUE)
}

add_season_rects <- function(p, ymin = -Inf, ymax = Inf, alpha = 0.45) {
  p +
    ggplot2::geom_rect(
      data = SEASON_WINDOWS,
      ggplot2::aes(xmin = .data$start, xmax = .data$end, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE,
      fill = SEASON_WINDOWS$fill,
      alpha = alpha,
      color = NA
    )
}

add_summer_lines <- function(p) {
  p +
    ggplot2::geom_vline(xintercept = SUMMER_START, linetype = "22", linewidth = 0.35) +
    ggplot2::geom_vline(xintercept = SUMMER_END, linetype = "22", linewidth = 0.35)
}

add_drought_segments <- function(p, y) {
  p +
    ggplot2::geom_segment(
      data = DROUGHT_WINDOWS,
      ggplot2::aes(x = .data$start, xend = .data$end, y = y, yend = y),
      inherit.aes = FALSE,
      color = "#E69F00",
      linewidth = 1.25,
      lineend = "round"
    )
}

summary_se <- function(data, group_cols, value_col) {
  data %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::summarise(
      mean = mean(.data[[value_col]], na.rm = TRUE),
      sd = stats::sd(.data[[value_col]], na.rm = TRUE),
      n = sum(!is.na(.data[[value_col]])),
      se = dplyr::if_else(.data$n > 1, .data$sd / sqrt(.data$n), 0),
      .groups = "drop"
    )
}

read_latest_temporal_effects <- function(data_name, resp_var, species) {
  alinv_read_temporal_bootstrap_effects(
    data_name = data_name,
    resp_var = resp_var,
    species = species,
    project_root = ALINV_PROJECT_ROOT
  )
}

compute_symmetric_limits <- function(df, cols = c("lower", "upper", "estimate"), pad = 0.08, floor = 0.5) {
  vals <- unlist(df[intersect(cols, names(df))], use.names = FALSE)
  vals <- vals[is.finite(vals)]
  if (!length(vals)) {
    return(c(-floor, floor))
  }
  lim <- max(abs(vals), na.rm = TRUE) * (1 + pad)
  lim <- max(lim, floor)
  c(-lim, lim)
}

compute_range_limits <- function(df, cols = c("lower", "upper", "estimate"), pad = 0.08, floor = 0.25) {
  vals <- unlist(df[intersect(cols, names(df))], use.names = FALSE)
  vals <- vals[is.finite(vals)]
  if (!length(vals)) {
    return(c(-floor, floor))
  }
  rng <- range(c(vals, 0), na.rm = TRUE)
  span <- diff(rng)
  if (!is.finite(span) || span <= 0) span <- floor
  rng + c(-1, 1) * span * pad
}

build_temporal_effect_panel <- function(effects_df, panel_title, y_limits, show_x = TRUE, show_y = TRUE) {
  effect_order <- c("culture", "precipitation", "robinia")
  effects_plot <- effects_df %>%
    dplyr::filter(.data$effect %in% effect_order) %>%
    dplyr::mutate(
      effect = factor(.data$effect, levels = effect_order),
      effect_label = dplyr::recode(as.character(.data$effect), !!!EFFECT_LABELS_SHORT)
    )

  drought_y <- y_limits[[1]] + diff(y_limits) * 0.04

  p <- ggplot2::ggplot(
    effects_plot,
    ggplot2::aes(
      x = .data$date,
      y = .data$estimate,
      ymin = .data$lower,
      ymax = .data$upper,
      color = .data$effect,
      fill = .data$effect,
      group = .data$effect
    )
  ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "22", linewidth = 0.25, color = "grey35") +
    ggplot2::geom_vline(xintercept = SUMMER_START, linetype = "42", linewidth = 0.35) +
    ggplot2::geom_vline(xintercept = SUMMER_END, linetype = "42", linewidth = 0.35) +
    ggplot2::geom_ribbon(alpha = 0.16, color = NA) +
    ggplot2::geom_line(linewidth = 0.55) +
    ggplot2::geom_point(size = 0.8) +
    ggplot2::scale_color_manual(
      values = EFFECT_COLORS[effect_order],
      breaks = effect_order,
      labels = unname(EFFECT_LABELS_SHORT[effect_order]),
      name = NULL
    ) +
    ggplot2::scale_fill_manual(
      values = EFFECT_COLORS[effect_order],
      breaks = effect_order,
      labels = unname(EFFECT_LABELS_SHORT[effect_order]),
      name = NULL
    ) +
    ggplot2::scale_x_date(date_breaks = "1 month", date_labels = "%b") +
    ggplot2::coord_cartesian(ylim = y_limits) +
    ggplot2::labs(
      x = NULL,
      y = if (isTRUE(show_y)) "Effect size (SD units; 95% bootstrap CI)" else NULL,
      title = panel_title
    ) +
    theme_alinv_pub(base_size = 7) +
    ggplot2::theme(
      axis.text.x = if (isTRUE(show_x)) ggplot2::element_text() else ggplot2::element_blank(),
      axis.ticks.x = if (isTRUE(show_x)) ggplot2::element_line(linewidth = 0.25) else ggplot2::element_blank(),
      axis.text.y = if (isTRUE(show_y)) ggplot2::element_text() else ggplot2::element_blank(),
      axis.ticks.y = if (isTRUE(show_y)) ggplot2::element_line(linewidth = 0.25) else ggplot2::element_blank(),
      legend.position = "bottom",
      plot.margin = ggplot2::margin(2, 4, 2, 16)
    )

  add_drought_segments(p, drought_y)
}

make_temporal_effect_figure <- function(data_name, resp_var, panel_suffix, output_file, height_mm) {
  effects <- dplyr::bind_rows(
    read_latest_temporal_effects(data_name, resp_var, "fagus"),
    read_latest_temporal_effects(data_name, resp_var, "quercus")
  )

  y_limits <- compute_range_limits(effects, pad = 0.08)

  p_fagus <- effects %>%
    dplyr::filter(.data$species == "fagus") %>%
    build_temporal_effect_panel(paste0("Fagus (", panel_suffix, ")"), y_limits, show_x = FALSE)

  p_quercus <- effects %>%
    dplyr::filter(.data$species == "quercus") %>%
    build_temporal_effect_panel(paste0("Quercus (", panel_suffix, ")"), y_limits, show_x = TRUE)

  fig <- (p_fagus / p_quercus) +
    patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 6.3),
      legend.key.width = grid::unit(7, "mm")
    )

  alinv_save_pdf(fig, output_file, height_mm = height_mm)
}

alinv_final_setup_loaded <- TRUE
