#!/usr/bin/env Rscript

# Supplementary Figure S1: observed spring phenology progression.

options(stringsAsFactors = FALSE)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this file with Rscript.", call. = FALSE)
SUPP_SCRIPT_FILE <- normalizePath(
  sub("^--file=", "", script_arg[[1]]),
  winslash = "/",
  mustWork = TRUE
)
source(file.path(dirname(SUPP_SCRIPT_FILE), "_phenology-setup.R"))

phenology_observed <- get_data("tree", "phenology_transitions") %>%
  alinv_filter_by_soil(soil_filter = SOIL_FILTER) %>%
  dplyr::filter(
    .data$species %in% SPECIES,
    .data$stage %in% 1:4,
    !is.na(.data$doy),
    !is.na(.data$stage_date)
  ) %>%
  dplyr::mutate(
    species = factor(
      .data$species,
      levels = SPECIES,
      labels = unname(SPECIES_LABELS[SPECIES])
    ),
    robinia = factor(
      .data$robinia,
      levels = alinv_factor_levels("robinia"),
      labels = c("Without Robinia", "With Robinia")
    ),
    precipitation = factor(
      .data$precipitation,
      levels = alinv_factor_levels("precipitation")
    ),
    culture = factor(
      .data$culture,
      levels = alinv_factor_levels("culture")
    ),
    stage = as.integer(.data$stage)
  )

progression_source_data <- phenology_observed %>%
  dplyr::group_by(
    .data$species,
    .data$robinia,
    .data$precipitation,
    .data$culture,
    .data$stage
  ) %>%
  dplyr::summarise(
    mean_doy = mean(.data$doy),
    sd_doy = stats::sd(.data$doy),
    n = dplyr::n(),
    se_doy = .data$sd_doy / sqrt(.data$n),
    minimum_doy = min(.data$doy),
    maximum_doy = max(.data$doy),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    .data$species,
    .data$robinia,
    .data$precipitation,
    .data$culture,
    .data$stage
  )

figure_s1 <- ggplot2::ggplot(
  progression_source_data,
  ggplot2::aes(
    x = .data$mean_doy,
    y = .data$stage,
    color = .data$precipitation,
    linetype = .data$culture,
    group = interaction(.data$precipitation, .data$culture)
  )
) +
  ggplot2::geom_segment(
    ggplot2::aes(
      x = .data$mean_doy - .data$se_doy,
      xend = .data$mean_doy + .data$se_doy,
      yend = .data$stage
    ),
    linewidth = 0.55,
    alpha = 0.55,
    lineend = "round"
  ) +
  ggplot2::geom_line(linewidth = 0.7) +
  ggplot2::geom_point(size = 1.8) +
  ggplot2::facet_grid(.data$species ~ .data$robinia) +
  ggplot2::scale_color_manual(
    values = PRECIPITATION_COLORS,
    labels = c(control = "Control", drought = "Reduced precipitation")
  ) +
  ggplot2::scale_linetype_manual(
    values = c(mono = "solid", mixed = "dotted"),
    labels = c(mono = "Monoculture", mixed = "Mixed culture")
  ) +
  ggplot2::scale_x_continuous(breaks = seq(90, 130, by = 10)) +
  ggplot2::scale_y_continuous(
    breaks = 1:4,
    labels = paste("Stage", 1:4),
    limits = c(0.8, 4.2)
  ) +
  ggplot2::labs(
    title = "Observed spring phenology progression",
    subtitle = "Treatment-group mean transition day of year (DOY) ± SE",
    x = "Day of year (DOY)",
    y = "Phenological stage",
    color = "Precipitation",
    linetype = "Culture",
    caption = "Stage 1 is shown descriptively; statistical inference uses transitions into stages 2–4."
  ) +
  theme_supplement(base_size = 8) +
  ggplot2::theme(
    panel.grid.major.y = ggplot2::element_line(color = "grey90", linewidth = 0.2),
    panel.grid.major.x = ggplot2::element_line(color = "grey93", linewidth = 0.18),
    legend.box = "horizontal"
  )

save_supplementary_plot(
  figure_s1,
  stem = "figure-s1-phenology-progression",
  width_mm = 160,
  height_mm = 112
)

message("Saved Figure S1 in: ", normalizePath(SUPP_OUTPUT_DIR, winslash = "/", mustWork = TRUE))
