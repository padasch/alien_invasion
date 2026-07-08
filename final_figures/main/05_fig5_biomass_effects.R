latest_biomass_workbook <- function() {
  candidates <- list.files(
    file.path(ALINV_PROJECT_ROOT, "data", "raw"),
    pattern = "^data_[0-9]+[.]xlsx$",
    full.names = TRUE
  )
  if (!length(candidates)) {
    stop("No biomass workbook matching data_*.xlsx found under data/raw.", call. = FALSE)
  }
  candidates[order(basename(candidates), decreasing = TRUE)][[1]]
}

prepare_fig5_effects <- function() {
  biomass <- suppressWarnings(wrangle_tree_biomass(latest_biomass_workbook()))
  effects <- dplyr::bind_rows(
    suppressMessages(run_biomass_glmm_species(
      biomass,
      species_keep = "fagus",
      soil_type = "both",
      include_soil_treatment = FALSE
    )$effects),
    suppressMessages(run_biomass_glmm_species(
      biomass,
      species_keep = "quercus",
      soil_type = "both",
      include_soil_treatment = FALSE
    )$effects)
  )

  term_labels <- c(
    `robiniawith-robinia` = "Robinia: without -> with",
    precipitationdrought = "Precipitation: control -> reduced",
    culturemixed = "Culture: mono -> mixed"
  )
  metric_labels <- c(
    shoot_biomass = "Shoot biomass",
    root_biomass = "Root biomass",
    root_shoot_biomass = "Root:shoot biomass"
  )

  effects %>%
    dplyr::mutate(
      term_label = dplyr::recode(.data$term, !!!term_labels),
      term_label = factor(
        .data$term_label,
        levels = c(
          "Culture: mono -> mixed",
          "Precipitation: control -> reduced",
          "Robinia: without -> with"
        )
      ),
      metric_label = factor(
        dplyr::recode(.data$metric, !!!metric_labels),
        levels = c("Shoot biomass", "Root biomass", "Root:shoot biomass")
      ),
      species_label = factor(
        dplyr::recode(.data$species, !!!SPECIES_LABELS),
        levels = c("Fagus", "Quercus")
      ),
      significant = factor(.data$significant, levels = c(FALSE, TRUE))
    )
}

make_fig5 <- function() {
  effects <- prepare_fig5_effects()
  x_limits <- compute_symmetric_limits(
    dplyr::rename(effects, lower = ci_lo, upper = ci_hi),
    cols = c("lower", "upper", "estimate"),
    pad = 0.08,
    floor = 1
  )

  fig <- ggplot2::ggplot(
    effects,
    ggplot2::aes(x = .data$estimate, y = .data$term_label, color = .data$significant)
  ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "42", linewidth = 0.3, color = "grey45") +
    ggplot2::geom_segment(
      ggplot2::aes(x = .data$ci_lo, xend = .data$ci_hi, y = .data$term_label, yend = .data$term_label),
      linewidth = 0.5
    ) +
    ggplot2::geom_point(size = 1.65) +
    ggplot2::facet_grid(.data$metric_label ~ .data$species_label) +
    ggplot2::coord_cartesian(xlim = x_limits) +
    ggplot2::scale_x_continuous(breaks = scales::breaks_pretty(n = 5)) +
    ggplot2::scale_color_manual(
      values = c(`FALSE` = "#7A7A7A", `TRUE` = "#D65F5F"),
      labels = c(`FALSE` = "95% CI includes 0", `TRUE` = "95% CI excludes 0"),
      name = NULL,
      drop = FALSE
    ) +
    ggplot2::labs(x = "Standardized effect size (SD units, 95% CI)", y = NULL) +
    theme_alinv_pub(base_size = 7) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.major.y = ggplot2::element_blank(),
      strip.text.y = ggplot2::element_text(angle = 90),
      axis.text.y = ggplot2::element_text(size = 6.7),
      plot.margin = ggplot2::margin(3, 5, 3, 4)
    )

  alinv_save_pdf(fig, "fig5_biomass_effects.pdf", height_mm = 118)
}
