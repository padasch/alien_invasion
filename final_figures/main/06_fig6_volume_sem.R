FIG6_HEATMAP_NONSIGNIFICANT_MODE <- "no_coefficient"
FIG6_HEATMAP_NONSIGNIFICANT_MODES <- c("white", "no_coefficient")
FIG6_HEATMAP_AXIS_LAYOUT <- "variables_factors"
FIG6_HEATMAP_AXIS_LAYOUTS <- c("factors_variables", "variables_factors")
FIG6_MISSING_CELL_FILL <- "#D9D9D9"
FIG6_HEIGHT_MM <- 142
FIG6_TIMELINE_X_LIMITS <- c(as.Date("2025-03-01"), as.Date("2025-12-31"))
FIG6_TIMELINE_X_BREAKS <- as.Date(c("2025-05-01", "2025-08-01", "2025-11-01"))
FIG6_PHASE_FILLS <- c(
  `until June` = "#C9F7C5",
  `July-August` = "#A9C9AA",
  `September+` = "#F8F6D6"
)
FIG6_HEATMAP_TREATMENT_LABELS <- c(
  extreme_event = "Extreme event",
  culture = "Culture",
  robinia = "Robinia",
  precipitation = "Precipitation"
)

prepare_fig6_volume_summary <- function() {
  df_metric <- prepare_growth_metric_plot_data(
    resp_var = "volume_inc_phase_rel",
    soil_type = "both",
    include_soil_treatment = FALSE
  ) %>%
    dplyr::filter(
      .data$species %in% c("fagus", "quercus"),
      !is.na(.data$volume_inc_phase_rel)
    ) %>%
    dplyr::mutate(volume_inc_phase_rel_pct = .data$volume_inc_phase_rel * 100)

  summarize_growth_metric_plot_data(
    df_metric = df_metric,
    resp_var = "volume_inc_phase_rel_pct",
    within_phase = TRUE
  )
}

prepare_fig6_phase_windows <- function(summary_df) {
  phase_levels <- levels(summary_df$phase)
  if (is.null(phase_levels)) {
    phase_levels <- unique(as.character(summary_df$phase))
  }

  phase_starts <- summary_df %>%
    dplyr::filter(!is.na(.data$phase), !is.na(.data$date)) %>%
    dplyr::mutate(phase_chr = as.character(.data$phase)) %>%
    dplyr::group_by(.data$phase_chr) %>%
    dplyr::summarise(start = min(.data$date), .groups = "drop") %>%
    dplyr::filter(.data$phase_chr %in% phase_levels) %>%
    dplyr::arrange(match(.data$phase_chr, phase_levels))

  if (nrow(phase_starts) < 2L) {
    return(SEASON_WINDOWS)
  }

  transition_dates <- phase_starts$start[-1]
  tibble::tibble(
    phase = factor(phase_levels, levels = phase_levels),
    start = c(FIG6_TIMELINE_X_LIMITS[[1]], transition_dates),
    end = c(transition_dates, FIG6_TIMELINE_X_LIMITS[[2]]),
    fill = unname(FIG6_PHASE_FILLS[phase_levels])
  )
}

build_fig6_volume_panel <- function(summary_df, species_key, robinia_key, y_limits, show_y = TRUE, show_legend = TRUE) {
  df_plot <- summary_df %>%
    dplyr::filter(.data$species == .env$species_key, .data$robinia == .env$robinia_key)
  phase_windows <- prepare_fig6_phase_windows(summary_df)
  phase_transitions <- phase_windows$start[-1]

  ggplot2::ggplot(
    df_plot,
    ggplot2::aes(
      x = .data$date,
      y = .data$mean,
      ymin = .data$mean - .data$se,
      ymax = .data$mean + .data$se,
      color = .data$precipitation,
      fill = .data$precipitation,
      linetype = .data$culture,
      group = .data$line_group
    )
  ) +
    ggplot2::geom_rect(
      data = phase_windows,
      ggplot2::aes(xmin = .data$start, xmax = .data$end, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill = phase_windows$fill,
      alpha = 0.42,
      color = NA,
      show.legend = FALSE
    ) +
    ggplot2::geom_vline(xintercept = phase_transitions, linewidth = 0.25, linetype = "42", color = "grey55") +
    ggplot2::geom_ribbon(alpha = 0.14, color = NA, show.legend = FALSE) +
    ggplot2::geom_line(linewidth = 0.48, lineend = "round") +
    ggplot2::scale_color_manual(
      values = PRECIP_COLORS,
      labels = c(control = "Precipitation: control", drought = "Precipitation: reduced"),
      name = "Treatment",
      drop = FALSE
    ) +
    ggplot2::scale_fill_manual(
      values = PRECIP_COLORS,
      labels = c(control = "control", drought = "reduced"),
      name = NULL,
      guide = "none",
      drop = FALSE
    ) +
    ggplot2::scale_linetype_manual(
      values = c(mono = "solid", mixed = "42"),
      labels = c(mono = "Culture: monoculture", mixed = "Culture: mixed"),
      name = "Treatment",
      drop = FALSE
    ) +
    ggplot2::scale_x_date(
      limits = FIG6_TIMELINE_X_LIMITS,
      breaks = FIG6_TIMELINE_X_BREAKS,
      date_labels = "%b",
      expand = ggplot2::expansion(mult = 0, add = c(6, 6))
    ) +
    ggplot2::coord_cartesian(ylim = y_limits, clip = "off") +
    ggplot2::labs(
      x = NULL,
      y = if (isTRUE(show_y)) "Volume increment per phase (%)" else NULL,
      title = paste0(SPECIES_LABELS[[species_key]], ": ", ROBINIA_LABELS[[robinia_key]])
    ) +
    theme_alinv_pub(base_size = 6.4) +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        ncol = 2,
        byrow = TRUE,
        override.aes = list(linewidth = 0.8, alpha = 1)
      ),
      linetype = ggplot2::guide_legend(
        ncol = 2,
        byrow = TRUE,
        override.aes = list(linewidth = 0.8, color = "black", alpha = 1)
      )
    ) +
    ggplot2::theme(
      axis.text.y = if (isTRUE(show_y)) ggplot2::element_text() else ggplot2::element_blank(),
      axis.ticks.y = if (isTRUE(show_y)) ggplot2::element_line(linewidth = 0.25) else ggplot2::element_blank(),
      legend.position = if (isTRUE(show_legend)) "bottom" else "none",
      legend.box = "vertical",
      legend.spacing.y = grid::unit(0.6, "mm"),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      plot.margin = ggplot2::margin(2, 2.5, 2, 2.5)
    )
}

find_sem_matrix_file <- function(species, resp_var) {
  root <- file.path(ALINV_PROJECT_ROOT, "output")
  pattern <- paste0(
    "sem-tree-.*-", resp_var, "-", species,
    "-soil-both_without_soil_treatment-noInt-scaled-.*-all-swcMeas.*-matrix_data[.]csv$"
  )
  files <- list.files(root, pattern = pattern, recursive = TRUE, full.names = TRUE)
  files <- files[file.exists(files)]
  if (!length(files)) return(NA_character_)

  date_part <- stringr::str_extract(files, "[0-9]{4}-[0-9]{2}-[0-9]{2}")
  prefer_rfe <- grepl("rfeAIC2", files)
  ord <- order(as.Date(date_part), prefer_rfe, file.info(files)$mtime, decreasing = TRUE, na.last = TRUE)
  files <- files[ord]

  for (file in files) {
    cols <- tryCatch(
      names(readr::read_csv(file, show_col_types = FALSE, n_max = 0)),
      error = function(e) character()
    )
    if ("path_type" %in% cols) {
      return(normalizePath(file, winslash = "/", mustWork = TRUE))
    }
  }

  NA_character_
}

read_sem_total_for_response <- function(species, resp_var, response_label) {
  file <- find_sem_matrix_file(species, resp_var)
  if (is.na(file)) {
    return(tibble::tibble())
  }

  readr::read_csv(file, show_col_types = FALSE) %>%
    dplyr::filter(.data$path_type == "total", .data$treatment != "swc") %>%
    dplyr::transmute(
      species = .data$species,
      treatment = .data$treatment,
      response_label = response_label,
      estimate = .data$estimate,
      p_value = .data$p_value,
      source_file = file
    )
}

ensure_phenology_average_transition_sem <- function() {
  file <- find_latest_file("phenology-average-transition-sem-matrix-data[.]csv$")
  if (!is.na(file)) {
    cols <- tryCatch(
      names(readr::read_csv(file, show_col_types = FALSE, n_max = 0)),
      error = function(e) character()
    )
    if ("effect_scale" %in% cols) {
      scale_values <- tryCatch(
        readr::read_csv(file, show_col_types = FALSE, col_select = "effect_scale") %>%
          dplyr::pull(.data$effect_scale) %>%
          unique(),
        error = function(e) character()
      )
      if ("standardized_y_swc" %in% scale_values) {
        return(file)
      }
    }
  }

  source(file.path(ALINV_PROJECT_ROOT, "scripts", "10-plot-phenology-average-transition-sem.R"), local = new.env(parent = globalenv()))
  file <- find_latest_file("phenology-average-transition-sem-matrix-data[.]csv$")
  if (!is.na(file)) {
    return(file)
  }
  stop("Phenology average-transition SEM did not create matrix data.", call. = FALSE)
}

read_phenology_transition_sem_totals <- function() {
  file <- ensure_phenology_average_transition_sem()
  readr::read_csv(file, show_col_types = FALSE) %>%
    dplyr::filter(.data$path_type == "total", .data$treatment != "swc") %>%
    dplyr::transmute(
      species = .data$species,
      treatment = .data$treatment,
      response_label = "Phenology timing",
      estimate = .data$estimate,
      p_value = .data$p_value,
      source_file = file
    )
}

prepare_fig6_sem_data <- function() {
  response_specs <- tibble::tribble(
    ~resp_var, ~response_label,
    "volume", "Volume (total)",
    "volume_inc_phase_rel", "Volume (incr.)",
    "chl", "Chlorophyll",
    "condition", "Vitality",
    "qy", "Quantum yield",
    "remaining_green", "Senescence (%)",
    "chlavg", "Senescence (Chl)"
  )

  standard <- purrr::map_dfr(c("fagus", "quercus"), function(species_i) {
    purrr::pmap_dfr(response_specs, function(resp_var, response_label) {
      read_sem_total_for_response(species_i, resp_var, response_label)
    })
  })

  sem_df <- dplyr::bind_rows(read_phenology_transition_sem_totals(), standard)

  treatment_order <- c("extreme_event", "culture", "robinia", "precipitation")
  response_order <- c(
    "Phenology timing",
    "Volume (total)",
    "Volume (incr.)",
    "Chlorophyll",
    "Vitality",
    "Quantum yield",
    "Senescence (%)",
    "Senescence (Chl)"
  )

  sem_df %>%
    dplyr::filter(.data$treatment %in% treatment_order) %>%
    dplyr::mutate(
      response_label = dplyr::if_else(.data$response_label == "Phenology timing", "Leaf-out", .data$response_label)
    ) %>%
    tidyr::complete(
      species = c("fagus", "quercus"),
      treatment = treatment_order,
      response_label = c("Leaf-out", response_order[-1])
    ) %>%
    dplyr::mutate(
      treatment = factor(.data$treatment, levels = rev(treatment_order)),
      treatment_label = dplyr::recode(as.character(.data$treatment), !!!FIG6_HEATMAP_TREATMENT_LABELS),
      treatment_label = factor(
        .data$treatment_label,
        levels = unname(FIG6_HEATMAP_TREATMENT_LABELS[rev(treatment_order)])
      ),
      response_label = factor(
        .data$response_label,
        levels = c("Leaf-out", response_order[-1])
      ),
      significant = !is.na(.data$p_value) & .data$p_value < 0.05,
      alpha_value = dplyr::case_when(
        is.na(.data$estimate) ~ 0,
        .data$significant ~ 1,
        TRUE ~ 0.28
      )
    )
}

build_fig6_sem_heatmap <- function(sem_df,
                                   species_key,
                                   fill_limit,
                                   show_y = TRUE,
                                   show_legend = TRUE,
                                   nonsignificant_mode = FIG6_HEATMAP_NONSIGNIFICANT_MODE,
                                   axis_layout = FIG6_HEATMAP_AXIS_LAYOUT) {
  nonsignificant_mode <- match.arg(nonsignificant_mode, FIG6_HEATMAP_NONSIGNIFICANT_MODES)
  axis_layout <- match.arg(axis_layout, FIG6_HEATMAP_AXIS_LAYOUTS)
  df_plot <- sem_df %>% dplyr::filter(.data$species == .env$species_key)

  if (identical(axis_layout, "factors_variables")) {
    x_limits <- levels(sem_df$response_label)
    y_limits <- levels(sem_df$treatment_label)
    x_var <- "response_label"
    y_var <- "treatment_label"
    x_angle <- 42
    x_size <- 5.5
    y_size <- 5.7
  } else {
    x_limits <- levels(sem_df$treatment_label)
    y_limits <- rev(levels(sem_df$response_label))
    x_var <- "treatment_label"
    y_var <- "response_label"
    x_angle <- 35
    x_size <- 5.7
    y_size <- 5.5
  }

  df_plot <- df_plot %>%
    dplyr::mutate(
      missing_cell = !is.finite(.data$estimate),
      display_estimate = dplyr::case_when(
        .data$missing_cell ~ NA_real_,
        .env$nonsignificant_mode == "white" & !.data$significant ~ NA_real_,
        TRUE ~ .data$estimate
      ),
      heatmap_x = factor(.data[[x_var]], levels = x_limits),
      heatmap_y = factor(.data[[y_var]], levels = y_limits)
    )

  label_df <- df_plot %>%
    dplyr::filter(is.finite(.data$estimate), .data$significant) %>%
    dplyr::mutate(
      label = sprintf("%.2f", .data$estimate)
    )

  ggplot2::ggplot(
    df_plot,
    ggplot2::aes(x = .data$heatmap_x, y = .data$heatmap_y)
  ) +
    ggplot2::geom_tile(
      data = df_plot %>% dplyr::filter(.data$missing_cell),
      fill = FIG6_MISSING_CELL_FILL,
      color = NA
    ) +
    ggplot2::geom_tile(
      data = df_plot %>% dplyr::filter(!.data$missing_cell),
      fill = "white",
      color = NA
    ) +
    ggplot2::geom_tile(
      data = df_plot %>% dplyr::filter(is.finite(.data$display_estimate)),
      ggplot2::aes(fill = .data$display_estimate),
      color = NA
    ) +
    ggplot2::geom_text(
      data = label_df,
      ggplot2::aes(label = .data$label),
      color = "black",
      size = 1.75,
      fontface = "bold",
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#D65F5F",
      mid = "white",
      high = "#3C6E8F",
      midpoint = 0,
      limits = c(-fill_limit, fill_limit),
      na.value = "white",
      name = "Standardized Effect",
      guide = ggplot2::guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        barwidth = grid::unit(24, "mm"),
        barheight = grid::unit(2.5, "mm"),
        label.position = "bottom"
      )
    ) +
    ggplot2::scale_x_discrete(limits = x_limits, drop = FALSE) +
    ggplot2::scale_y_discrete(limits = y_limits, drop = FALSE) +
    ggplot2::labs(x = NULL, y = NULL, title = paste0(SPECIES_LABELS[[species_key]], ": total SEM effects")) +
    theme_alinv_pub(base_size = 6.2) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = x_angle, hjust = 1, vjust = 1, size = x_size),
      axis.text.y = if (isTRUE(show_y)) ggplot2::element_text(size = y_size) else ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_line(linewidth = 0.25, color = "black"),
      axis.ticks.y = if (isTRUE(show_y)) {
        ggplot2::element_line(linewidth = 0.25, color = "black")
      } else {
        ggplot2::element_blank()
      },
      axis.ticks.length = grid::unit(1.2, "mm"),
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      axis.line.x.bottom = ggplot2::element_line(linewidth = 0.25, color = "black"),
      axis.line.y.left = if (isTRUE(show_y)) {
        ggplot2::element_line(linewidth = 0.25, color = "black")
      } else {
        ggplot2::element_blank()
      },
      legend.position = if (isTRUE(show_legend)) "bottom" else "none",
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      plot.margin = ggplot2::margin(2, 2, 2, 2)
    )
}

build_fig6_line_legend <- function() {
  line_length <- 0.039
  label_gap <- 0.018
  legend_df <- tibble::tibble(
    x = c(0.18, 0.58, 0.18, 0.58),
    y = c(0.38, 0.38, 0.12, 0.12),
    label = c(
      "Culture: monoculture",
      "Culture: mixed",
      "Precipitation: control",
      "Precipitation: reduced"
    ),
    color = c("black", "black", unname(PRECIP_COLORS["control"]), unname(PRECIP_COLORS["drought"])),
    linetype = c("solid", "42", "solid", "solid")
  ) %>%
    dplyr::mutate(
      xend = .data$x + line_length,
      label_x = .data$xend + label_gap
    )

  ggplot2::ggplot(legend_df) +
    ggplot2::annotate(
      "text",
      x = 0.5,
      y = 0.66,
      label = "Treatment",
      fontface = "bold",
      size = 2.15
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = .data$x,
        xend = .data$xend,
        y = .data$y,
        yend = .data$y,
        color = .data$color,
        linetype = .data$linetype
      ),
      linewidth = 0.65,
      lineend = "round"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = .data$label_x, y = .data$y, label = .data$label),
      hjust = 0,
      vjust = 0.5,
      size = 2.05,
      color = "black"
    ) +
    ggplot2::scale_color_identity() +
    ggplot2::scale_linetype_identity() +
    ggplot2::coord_cartesian(xlim = c(0, 1.08), ylim = c(0, 0.76), clip = "off") +
    ggplot2::theme_void(base_size = 6.2) +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0))
}

build_fig6_heatmap_legend <- function(fill_limit) {
  legend_limit <- max(fill_limit, 1)
  bar_limit <- legend_limit * 0.4
  bar_x <- seq(-bar_limit, bar_limit, length.out = 220)
  bar_effect <- seq(-legend_limit, legend_limit, length.out = 220)
  bar_width <- diff(range(bar_x)) / (length(bar_x) - 1)
  tick_values <- c(-1, -0.5, 0, 0.5, 1)
  tick_values <- tick_values[tick_values >= -legend_limit & tick_values <= legend_limit]
  na_x <- bar_limit + 0.13 * legend_limit
  na_width <- 0.075 * legend_limit

  bar_df <- tibble::tibble(x = bar_x, effect = bar_effect, y = 1)
  tick_df <- tibble::tibble(
    effect = tick_values,
    x = tick_values / legend_limit * bar_limit,
    label = sprintf("%.1f", tick_values)
  )

  ggplot2::ggplot() +
    ggplot2::geom_tile(
      data = bar_df,
      ggplot2::aes(x = .data$x, y = .data$y, fill = .data$effect),
      width = bar_width * 1.05,
      height = 0.24,
      color = NA
    ) +
    ggplot2::geom_segment(
      data = tick_df,
      ggplot2::aes(x = .data$x, xend = .data$x, y = 0.88, yend = 1.12),
      color = "white",
      linewidth = 0.18
    ) +
    ggplot2::geom_text(
      data = tick_df,
      ggplot2::aes(x = .data$x, y = 0.67, label = .data$label),
      size = 1.95,
      color = "black"
    ) +
    ggplot2::annotate(
      "text",
      x = 0,
      y = 1.36,
      label = "Standardized Effect",
      size = 2.05,
      color = "black"
    ) +
    ggplot2::annotate(
      "tile",
      x = na_x,
      y = 1,
      width = na_width,
      height = 0.24,
      fill = FIG6_MISSING_CELL_FILL,
      color = NA
    ) +
    ggplot2::annotate(
      "text",
      x = na_x + 0.06 * legend_limit,
      y = 1,
      label = "NA",
      hjust = 0,
      vjust = 0.5,
      size = 2.05,
      color = "black"
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#D65F5F",
      mid = "white",
      high = "#3C6E8F",
      midpoint = 0,
      limits = c(-fill_limit, fill_limit),
      guide = "none"
    ) +
    ggplot2::theme_void(base_size = 6.2) +
    ggplot2::coord_cartesian(
      xlim = c(-legend_limit, legend_limit),
      ylim = c(0.48, 1.5),
      clip = "off"
    ) +
    ggplot2::theme(
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    )
}

make_fig6 <- function(nonsignificant_mode = FIG6_HEATMAP_NONSIGNIFICANT_MODE,
                      axis_layout = FIG6_HEATMAP_AXIS_LAYOUT) {
  nonsignificant_mode <- match.arg(nonsignificant_mode, FIG6_HEATMAP_NONSIGNIFICANT_MODES)
  axis_layout <- match.arg(axis_layout, FIG6_HEATMAP_AXIS_LAYOUTS)
  volume_summary <- prepare_fig6_volume_summary()
  y_limits <- compute_range_limits(
    volume_summary %>%
      dplyr::mutate(lower = .data$mean - .data$se, upper = .data$mean + .data$se, estimate = .data$mean),
    cols = c("lower", "upper", "estimate"),
    pad = 0.07,
    floor = 0.15
  )

  sem_df <- prepare_fig6_sem_data()
  fill_limit <- max(abs(sem_df$estimate), na.rm = TRUE)
  if (!is.finite(fill_limit) || fill_limit <= 0) fill_limit <- 1

  p_a <- build_fig6_volume_panel(volume_summary, "fagus", "without-robinia", y_limits, show_y = TRUE, show_legend = FALSE)
  p_b <- build_fig6_volume_panel(volume_summary, "fagus", "with-robinia", y_limits, show_y = FALSE, show_legend = FALSE)
  p_c <- build_fig6_sem_heatmap(sem_df, "fagus", fill_limit, show_y = TRUE, show_legend = FALSE, nonsignificant_mode = nonsignificant_mode, axis_layout = axis_layout)
  p_d <- build_fig6_volume_panel(volume_summary, "quercus", "without-robinia", y_limits, show_y = TRUE, show_legend = FALSE)
  p_e <- build_fig6_volume_panel(volume_summary, "quercus", "with-robinia", y_limits, show_y = FALSE, show_legend = FALSE)
  p_f <- build_fig6_sem_heatmap(sem_df, "quercus", fill_limit, show_y = TRUE, show_legend = FALSE, nonsignificant_mode = nonsignificant_mode, axis_layout = axis_layout)

  p_a <- p_a + ggplot2::labs(tag = "A")
  p_b <- p_b + ggplot2::labs(tag = "B")
  p_c <- p_c + ggplot2::labs(tag = "C")
  p_d <- p_d + ggplot2::labs(tag = "D")
  p_e <- p_e + ggplot2::labs(tag = "E")
  p_f <- p_f + ggplot2::labs(tag = "F")

  line_legend <- build_fig6_line_legend()
  heatmap_legend <- build_fig6_heatmap_legend(fill_limit)

  fig <- p_a + p_b + p_c + p_d + p_e + p_f +
    patchwork::wrap_elements(full = line_legend) +
    patchwork::wrap_elements(full = heatmap_legend) +
    patchwork::plot_layout(
      design = "
ABCC
DEFF
GGHH
",
      widths = c(1, 1, 0.59, 0.59),
      heights = c(1, 1, 0.13)
    )

  alinv_save_pdf(fig, "fig6_volume_sem.pdf", height_mm = FIG6_HEIGHT_MM)
}
