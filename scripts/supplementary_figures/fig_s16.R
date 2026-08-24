#!/usr/bin/env Rscript

# Supplementary Figure S2: primary common-shift phenology LMM with
# block-stratified container-cluster bootstrap uncertainty.

options(stringsAsFactors = FALSE)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this file with Rscript.", call. = FALSE)
SUPP_SCRIPT_FILE <- normalizePath(
  gsub("~\\+~", " ", sub("^--file=", "", script_arg[[1]])),
  winslash = "/",
  mustWork = TRUE
)
source(file.path(dirname(SUPP_SCRIPT_FILE), "_phenology-setup.R"))

primary_effects <- alinv_read_phenology_bootstrap_primary(PROJECT_ROOT) %>%
  dplyr::mutate(
    species_label = factor(
      unname(SPECIES_LABELS[.data$species]),
      levels = unname(SPECIES_LABELS[SPECIES])
    ),
    treatment_label = factor(
      .data$treatment_label,
      levels = c("Drought", "With robinia", "Mixed culture")
    ),
    interval_excludes_zero = .data$lower_95_oriented > 0 |
      .data$upper_95_oriented < 0
  ) %>%
  dplyr::arrange(.data$species, .data$plot_order)

figure_s2 <- ggplot2::ggplot(
  primary_effects,
  ggplot2::aes(
    x = .data$estimate_oriented,
    y = .data$treatment_label,
    color = .data$effect
  )
) +
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey40",
    linewidth = 0.4
  ) +
  ggplot2::geom_segment(
    ggplot2::aes(
      x = .data$lower_95_oriented,
      xend = .data$upper_95_oriented,
      yend = .data$treatment_label
    ),
    linewidth = 0.8,
    lineend = "round"
  ) +
  ggplot2::geom_point(
    ggplot2::aes(shape = .data$interval_excludes_zero),
    size = 2.5,
    fill = "white",
    stroke = 0.8
  ) +
  ggplot2::facet_wrap(~species_label, nrow = 1) +
  ggplot2::scale_color_manual(values = TREATMENT_COLORS, guide = "none") +
  ggplot2::scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), guide = "none") +
  ggplot2::labs(
    x = "Timing effect (days; negative = delayed, positive = earlier)",
    y = NULL
  ) +
  theme_supplement(base_size = 8) +
  ggplot2::theme(
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_line(color = "grey91", linewidth = 0.2)
  )

save_supplementary_plot(
  figure_s2,
  stem = "fig_s16",
  width_mm = 160,
  height_mm = 82
)

print(tibble::as_tibble(primary_effects), n = Inf)
message("Saved Figure S16 in: ", normalizePath(SUPP_OUTPUT_DIR, winslash = "/", mustWork = TRUE))
