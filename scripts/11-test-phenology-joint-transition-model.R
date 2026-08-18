#!/usr/bin/env Rscript

# Standalone sensitivity analysis for spring phenology transition dates.
#
# This script intentionally does not modify the production phenology workflow.
# It replaces the separate stage-specific models plus additive pooled model with
# one joint model per focal species:
#
#   DOY ~ stage * (precipitation + robinia + culture)
#       + (1 | container) + (1 | tree)
#
# Stage-specific and equal-stage averaged treatment contrasts are derived from
# the same fitted model. Likelihood-ratio tests assess whether each treatment
# effect differs among stages (the stage-by-treatment interactions).

options(stringsAsFactors = FALSE)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) {
  stop("Run this file with Rscript so the project root can be resolved.", call. = FALSE)
}

script_file <- normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = TRUE)
project_root <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/", mustWork = TRUE)
setwd(project_root)

renv_lib <- Sys.glob(file.path("renv", "library", "*", "R-*", "*"))
if (length(renv_lib)) {
  .libPaths(c(normalizePath(renv_lib[[1]], winslash = "/", mustWork = TRUE), .libPaths()))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(emmeans)
  library(ggplot2)
  library(lme4)
  library(patchwork)
  library(readr)
  library(tidyr)
})

suppressMessages(suppressPackageStartupMessages({
  source(file.path(project_root, "functions", "_source.R"))
  source(file.path(project_root, "functions", "1-summary-figures.R"))
  source(file.path(project_root, "functions", "7-phenology-transition-models.R"))
}))

analysis_date <- as.character(Sys.Date())
species_vec <- c("fagus", "quercus")
stages_keep <- 2:4
soil_type <- "both"
include_soil_treatment <- FALSE

alinv_set_analysis_context(
  scenario_label = "both soils (without soil as treatment)",
  soil_filter = soil_type,
  include_soil_treatment = include_soil_treatment,
  analysis_date = analysis_date,
  output_root = file.path(project_root, "output"),
  create_dirs = TRUE
)

output_dir <- file.path(
  project_root,
  "output",
  analysis_date,
  "phenology-joint-transition-model-test"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

treatment_specs <- ALINV_TREATMENT_CONFIG %>%
  dplyr::filter(.data$effect %in% c("precipitation", "robinia", "culture")) %>%
  dplyr::arrange(.data$plot_order) %>%
  dplyr::select(
    "effect",
    "baseline_level",
    "treatment_level",
    treatment_label = "short_label",
    "contrast_label",
    "plot_order"
  )

species_labels <- c(fagus = "Fagus", quercus = "Quercus")
treatment_colors <- c(
  precipitation = "#1B9E77",
  robinia = "#D95F02",
  culture = "#7570B3"
)
precipitation_colors <- c(control = "#4F6674", drought = "#D65F5F")

joint_formula <- stats::as.formula(
  "doy ~ stage_label * (precipitation + robinia + culture) + (1 | boxlabel) + (1 | tree_id)"
)

model_control <- lme4::lmerControl(
  optimizer = "bobyqa",
  optCtrl = list(maxfun = 200000)
)

fit_joint_model <- function(species_keep) {
  df <- prepare_phenology_transition_data(
    species_keep = species_keep,
    soil_type = soil_type,
    include_soil_treatment = include_soil_treatment,
    stages_keep = stages_keep
  ) %>%
    dplyr::mutate(
      stage_label = factor(
        as.character(.data$stage_label),
        levels = paste0("Stage ", stages_keep)
      )
    ) %>%
    droplevels()

  fit_reml <- lme4::lmer(
    joint_formula,
    data = df,
    REML = TRUE,
    control = model_control
  )

  fit_ml <- lme4::lmer(
    joint_formula,
    data = df,
    REML = FALSE,
    control = model_control
  )

  list(data = df, fit_reml = fit_reml, fit_ml = fit_ml)
}

extract_contrast <- function(fit, species_keep, spec_row) {
  effect_i <- spec_row$effect[[1]]
  expected_levels <- c(
    spec_row$baseline_level[[1]],
    spec_row$treatment_level[[1]]
  )

  model_levels <- fit@frame[[effect_i]] %>% levels()
  if (!identical(model_levels, expected_levels)) {
    stop(
      "Unexpected factor order for ", effect_i, ": ",
      paste(model_levels, collapse = ", "),
      call. = FALSE
    )
  }

  stage_formula <- stats::as.formula(paste("~", effect_i, "| stage_label"))
  emm_stage <- emmeans::emmeans(
    fit,
    specs = stage_formula,
    weights = "equal",
    lmer.df = "kenward-roger"
  )
  stage_contrast <- emmeans::contrast(
    emm_stage,
    method = "revpairwise",
    by = "stage_label",
    adjust = "none"
  ) %>%
    summary(infer = c(TRUE, TRUE), level = 0.95, adjust = "none") %>%
    as.data.frame()

  average_formula <- stats::as.formula(paste("~", effect_i))
  emm_average <- suppressMessages(
    emmeans::emmeans(
      fit,
      specs = average_formula,
      weights = "equal",
      lmer.df = "kenward-roger"
    )
  )
  average_contrast <- emmeans::contrast(
    emm_average,
    method = "revpairwise",
    adjust = "none"
  ) %>%
    summary(infer = c(TRUE, TRUE), level = 0.95, adjust = "none") %>%
    as.data.frame() %>%
    dplyr::mutate(stage_label = "Average (Stages 2-4)")

  dplyr::bind_rows(stage_contrast, average_contrast) %>%
    dplyr::transmute(
      species = species_keep,
      stage_label = as.character(.data$stage_label),
      effect = effect_i,
      treatment_label = spec_row$treatment_label[[1]],
      contrast_label = spec_row$contrast_label[[1]],
      estimate_days = .data$estimate,
      se = .data$SE,
      df = .data$df,
      lower_95 = .data$lower.CL,
      upper_95 = .data$upper.CL,
      test_statistic = .data$t.ratio,
      p_value = .data$p.value,
      timing_direction = dplyr::case_when(
        .data$estimate_days < 0 ~ "earlier",
        .data$estimate_days > 0 ~ "later",
        TRUE ~ "no change"
      ),
      plot_order = spec_row$plot_order[[1]]
    )
}

extract_all_contrasts <- function(model_result, species_keep) {
  dplyr::bind_rows(lapply(seq_len(nrow(treatment_specs)), function(i) {
    extract_contrast(
      fit = model_result$fit_reml,
      species_keep = species_keep,
      spec_row = treatment_specs[i, ]
    )
  }))
}

test_stage_interactions <- function(model_result, species_keep) {
  full_ml <- model_result$fit_ml
  df <- model_result$data

  dplyr::bind_rows(lapply(seq_len(nrow(treatment_specs)), function(i) {
    spec_i <- treatment_specs[i, ]
    effect_i <- spec_i$effect[[1]]
    interaction_term <- paste0("stage_label:", effect_i)
    reduced_formula <- stats::update(
      stats::formula(full_ml),
      paste(". ~ . -", interaction_term)
    )
    reduced_ml <- lme4::lmer(
      reduced_formula,
      data = df,
      REML = FALSE,
      control = model_control
    )
    comparison <- stats::anova(reduced_ml, full_ml, refit = FALSE)

    tibble::tibble(
      species = species_keep,
      effect = effect_i,
      treatment_label = spec_i$treatment_label[[1]],
      tested_term = interaction_term,
      df = comparison$Df[[2]],
      chi_square = comparison$Chisq[[2]],
      p_value = comparison$`Pr(>Chisq)`[[2]],
      aic_reduced = comparison$AIC[[1]],
      aic_full = comparison$AIC[[2]],
      plot_order = spec_i$plot_order[[1]]
    )
  }))
}

extract_fixed_effects <- function(model_result, species_keep) {
  coef_table <- as.data.frame(summary(model_result$fit_reml)$coefficients)
  coef_table$term <- rownames(coef_table)
  rownames(coef_table) <- NULL

  coef_table %>%
    dplyr::transmute(
      species = species_keep,
      term = .data$term,
      estimate = .data$Estimate,
      se = .data$`Std. Error`,
      t_value = .data$`t value`
    )
}

extract_model_summary <- function(model_result, species_keep) {
  fit <- model_result$fit_reml
  df <- model_result$data
  convergence_messages <- fit@optinfo$conv$lme4$messages
  if (is.null(convergence_messages)) {
    convergence_messages <- "none"
  } else {
    convergence_messages <- paste(convergence_messages, collapse = "; ")
  }

  tibble::tibble(
    species = species_keep,
    formula = paste(deparse(stats::formula(fit)), collapse = " "),
    n_observations = stats::nobs(fit),
    n_trees = dplyr::n_distinct(df$tree_id),
    n_containers = dplyr::n_distinct(df$boxlabel),
    n_stage_2 = sum(df$stage_label == "Stage 2"),
    n_stage_3 = sum(df$stage_label == "Stage 3"),
    n_stage_4 = sum(df$stage_label == "Stage 4"),
    singular = lme4::isSingular(fit, tol = 1e-4),
    convergence_messages = convergence_messages,
    sigma = stats::sigma(fit),
    REML_criterion = stats::deviance(fit, REML = TRUE)
  )
}

extract_kenward_roger_tests <- function(model_result, species_keep) {
  emmeans::joint_tests(
    model_result$fit_reml,
    lmer.df = "kenward-roger"
  ) %>%
    as.data.frame() %>%
    dplyr::transmute(
      species = species_keep,
      model_term = .data$`model term`,
      numerator_df = .data$df1,
      denominator_df = .data$df2,
      f_ratio = .data$F.ratio,
      p_value = .data$p.value
    )
}

model_results <- stats::setNames(lapply(species_vec, fit_joint_model), species_vec)

contrast_results <- dplyr::bind_rows(lapply(species_vec, function(species_i) {
  extract_all_contrasts(model_results[[species_i]], species_i)
}))

interaction_tests <- dplyr::bind_rows(lapply(species_vec, function(species_i) {
  test_stage_interactions(model_results[[species_i]], species_i)
}))

fixed_effects <- dplyr::bind_rows(lapply(species_vec, function(species_i) {
  extract_fixed_effects(model_results[[species_i]], species_i)
}))

model_summary <- dplyr::bind_rows(lapply(species_vec, function(species_i) {
  extract_model_summary(model_results[[species_i]], species_i)
}))

kenward_roger_tests <- dplyr::bind_rows(lapply(species_vec, function(species_i) {
  extract_kenward_roger_tests(model_results[[species_i]], species_i)
}))

diagnostic_data <- dplyr::bind_rows(lapply(species_vec, function(species_i) {
  fit <- model_results[[species_i]]$fit_reml
  df <- model_results[[species_i]]$data
  tibble::tibble(
    species = species_i,
    stage_label = df$stage_label,
    fitted = stats::fitted(fit),
    residual = stats::residuals(fit, type = "pearson")
  )
}))

residual_spread <- diagnostic_data %>%
  dplyr::group_by(.data$species, .data$stage_label) %>%
  dplyr::summarise(
    n = dplyr::n(),
    residual_sd = stats::sd(.data$residual),
    residual_mad = stats::mad(.data$residual),
    .groups = "drop"
  )

model_cell_counts <- dplyr::bind_rows(lapply(species_vec, function(species_i) {
  model_results[[species_i]]$data %>%
    dplyr::count(
      .data$stage_label,
      .data$precipitation,
      .data$robinia,
      .data$culture,
      name = "n"
    ) %>%
    dplyr::mutate(species = species_i, .before = 1)
}))

# Panel A: descriptive mean transition dates for context. Stage 1 is displayed
# but was deliberately excluded from the joint inferential model.
descriptive_data <- get_data("tree", "phenology_transitions") %>%
  alinv_filter_by_soil(soil_filter = soil_type) %>%
  dplyr::filter(
    .data$species %in% species_vec,
    .data$stage %in% 1:4,
    !is.na(.data$doy),
    !is.na(.data$stage_date)
  ) %>%
  dplyr::mutate(
    species = factor(.data$species, levels = species_vec, labels = unname(species_labels[species_vec])),
    robinia = factor(.data$robinia, levels = alinv_factor_levels("robinia")),
    precipitation = factor(.data$precipitation, levels = alinv_factor_levels("precipitation")),
    culture = factor(.data$culture, levels = alinv_factor_levels("culture")),
    stage = factor(as.integer(.data$stage), levels = 1:4, labels = paste("Stage", 1:4))
  ) %>%
  dplyr::group_by(.data$species, .data$robinia, .data$precipitation, .data$culture, .data$stage) %>%
  dplyr::summarise(
    mean_doy = mean(.data$doy),
    sd_doy = stats::sd(.data$doy),
    n = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    se_doy = .data$sd_doy / sqrt(.data$n),
    stage_number = as.numeric(.data$stage),
    robinia = factor(
      .data$robinia,
      levels = alinv_factor_levels("robinia"),
      labels = c("without Robinia", "with Robinia")
    )
  )

plot_descriptive <- ggplot2::ggplot(
  descriptive_data,
  ggplot2::aes(
    x = .data$mean_doy,
    y = .data$stage_number,
    color = .data$precipitation,
    linetype = .data$culture,
    group = interaction(.data$precipitation, .data$culture)
  )
) +
  ggplot2::geom_segment(
    ggplot2::aes(
      x = .data$mean_doy - .data$se_doy,
      xend = .data$mean_doy + .data$se_doy,
      yend = .data$stage_number
    ),
    linewidth = 0.65,
    alpha = 0.55,
    lineend = "round"
  ) +
  ggplot2::geom_line(linewidth = 0.75) +
  ggplot2::geom_point(size = 2.0) +
  ggplot2::facet_grid(.data$species ~ .data$robinia) +
  ggplot2::scale_color_manual(
    values = precipitation_colors,
    labels = c(control = "Control", drought = "Reduced precipitation")
  ) +
  ggplot2::scale_linetype_manual(
    values = c(mono = "solid", mixed = "dotted"),
    labels = c(mono = "Monoculture", mixed = "Mixed culture")
  ) +
  ggplot2::scale_y_continuous(
    breaks = 1:4,
    labels = paste("Stage", 1:4),
    limits = c(0.8, 4.2)
  ) +
  ggplot2::labs(
    title = "Observed phenological transition timing",
    subtitle = "Treatment-group means ± SE; stage 1 is descriptive only and is not included in the joint model.",
    x = "Day of year (DOY)",
    y = "Phenological stage",
    color = "Precipitation",
    linetype = "Culture"
  ) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(
    panel.grid.major = ggplot2::element_line(color = "grey90", linewidth = 0.25),
    panel.grid.minor = ggplot2::element_blank(),
    strip.background = ggplot2::element_rect(fill = "grey12", color = "grey12"),
    strip.text = ggplot2::element_text(color = "white", face = "bold"),
    legend.position = "bottom",
    plot.title = ggplot2::element_text(face = "bold")
  )

stage_levels <- c(paste0("Stage ", stages_keep), "Average (Stages 2-4)")
treatment_levels <- treatment_specs$treatment_label

contrast_plot_data <- contrast_results %>%
  dplyr::mutate(
    species = factor(.data$species, levels = species_vec, labels = unname(species_labels[species_vec])),
    stage_label = factor(.data$stage_label, levels = stage_levels),
    treatment_label = factor(.data$treatment_label, levels = treatment_levels),
    significant = .data$p_value < 0.05
  )

effect_limit <- max(abs(c(contrast_plot_data$lower_95, contrast_plot_data$upper_95)), na.rm = TRUE)
effect_limit <- ceiling(effect_limit * 2) / 2

plot_contrasts <- ggplot2::ggplot(
  contrast_plot_data,
  ggplot2::aes(
    x = .data$estimate_days,
    y = .data$treatment_label,
    color = .data$effect
  )
) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.45) +
  ggplot2::geom_segment(
    ggplot2::aes(
      x = .data$lower_95,
      xend = .data$upper_95,
      yend = .data$treatment_label
    ),
    linewidth = 0.9,
    alpha = 0.8,
    lineend = "round"
  ) +
  ggplot2::geom_point(
    ggplot2::aes(shape = .data$significant),
    size = 2.7,
    stroke = 0.9,
    fill = "white"
  ) +
  ggplot2::facet_grid(.data$species ~ .data$stage_label) +
  ggplot2::coord_cartesian(xlim = c(-effect_limit, effect_limit)) +
  ggplot2::scale_color_manual(values = treatment_colors, guide = "none") +
  ggplot2::scale_shape_manual(
    values = c(`FALSE` = 16, `TRUE` = 21),
    labels = c(`FALSE` = "P ≥ 0.05", `TRUE` = "P < 0.05")
  ) +
  ggplot2::labs(
    title = "Treatment contrasts from one joint transition-date model per species",
    subtitle = paste(
      "Points are treatment minus reference estimates with unadjusted 95% Kenward–Roger intervals.",
      "The Average column is an equal-stage marginal contrast from the same model."
    ),
    x = "Shift in transition DOY under treatment (days)",
    y = NULL,
    shape = "Stage contrast"
  ) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(
    panel.grid.major.y = ggplot2::element_line(color = "grey92", linewidth = 0.25),
    strip.background = ggplot2::element_rect(fill = "grey12", color = "grey12"),
    strip.text = ggplot2::element_text(color = "white", face = "bold"),
    legend.position = "bottom",
    plot.title = ggplot2::element_text(face = "bold")
  )

interaction_plot_data <- interaction_tests %>%
  dplyr::mutate(
    species = factor(.data$species, levels = species_vec, labels = unname(species_labels[species_vec])),
    treatment_label = factor(.data$treatment_label, levels = treatment_levels),
    evidence = factor(
      ifelse(.data$p_value < 0.05, "P < 0.05", "P ≥ 0.05"),
      levels = c("P ≥ 0.05", "P < 0.05")
    ),
    cell_label = sprintf("χ²(%d) = %.2f\nP = %s", .data$df, .data$chi_square, format.pval(.data$p_value, digits = 2, eps = 0.001))
  )

plot_interactions <- ggplot2::ggplot(
  interaction_plot_data,
  ggplot2::aes(x = .data$treatment_label, y = .data$species, fill = .data$evidence)
) +
  ggplot2::geom_tile(color = "grey70", linewidth = 0.45) +
  ggplot2::geom_text(ggplot2::aes(label = .data$cell_label), size = 3.4, lineheight = 0.95) +
  ggplot2::scale_fill_manual(
    values = c("P ≥ 0.05" = "white", "P < 0.05" = "#E6B44A"),
    drop = FALSE
  ) +
  ggplot2::labs(
    title = "Do treatment effects differ among stages?",
    subtitle = "Likelihood-ratio tests compare the full ML model with a model omitting each stage × treatment interaction.",
    x = NULL,
    y = NULL,
    fill = "Interaction test"
  ) +
  ggplot2::coord_fixed(ratio = 0.48) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.text = ggplot2::element_text(color = "black"),
    legend.position = "bottom",
    plot.title = ggplot2::element_text(face = "bold")
  )

summary_plot <- (
  plot_descriptive /
    plot_contrasts /
    patchwork::wrap_elements(full = plot_interactions)
) +
  patchwork::plot_layout(heights = c(1.05, 1.15, 0.55)) +
  patchwork::plot_annotation(
    title = "Phenology treatment assessment using one joint transition-date model",
    subtitle = paste(
      "Separate models are fitted for Fagus and Quercus; within each species, stages 2–4 are analysed together.",
      "Negative contrasts indicate earlier and positive contrasts later stage attainment."
    ),
    tag_levels = "A",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 17),
      plot.subtitle = ggplot2::element_text(size = 11)
    )
  ) &
  ggplot2::theme(plot.margin = ggplot2::margin(t = 8, r = 10, b = 8, l = 10))

diagnostic_plot <- (
  ggplot2::ggplot(
    diagnostic_data,
    ggplot2::aes(x = .data$fitted, y = .data$residual, color = .data$stage_label)
  ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    ggplot2::geom_point(alpha = 0.28, size = 1.0) +
    ggplot2::facet_wrap(~species, scales = "free_x", labeller = ggplot2::as_labeller(species_labels)) +
    ggplot2::scale_color_manual(values = c("Stage 2" = "#4F6674", "Stage 3" = "#D95F02", "Stage 4" = "#7570B3")) +
    ggplot2::labs(title = "Pearson residuals versus fitted transition DOY", x = "Fitted DOY", y = "Pearson residual", color = "Stage") +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"), legend.position = "bottom")
) / (
  ggplot2::ggplot(diagnostic_data, ggplot2::aes(sample = .data$residual, color = .data$stage_label)) +
    ggplot2::stat_qq(alpha = 0.28, size = 1.0) +
    ggplot2::stat_qq_line(color = "grey25", linewidth = 0.55) +
    ggplot2::facet_wrap(~species, labeller = ggplot2::as_labeller(species_labels)) +
    ggplot2::scale_color_manual(values = c("Stage 2" = "#4F6674", "Stage 3" = "#D95F02", "Stage 4" = "#7570B3")) +
    ggplot2::labs(title = "Normal Q–Q plot of Pearson residuals", x = "Theoretical quantile", y = "Observed quantile", color = "Stage") +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"), legend.position = "bottom")
) +
  patchwork::plot_annotation(
    title = "Joint phenology transition model diagnostics",
    subtitle = "Colour identifies the phenological stage; models were fitted separately by species."
  )

summary_png <- file.path(output_dir, "phenology-joint-transition-model-summary.png")
summary_pdf <- file.path(output_dir, "phenology-joint-transition-model-summary.pdf")
diagnostic_png <- file.path(output_dir, "phenology-joint-transition-model-diagnostics.png")

ggplot2::ggsave(summary_png, summary_plot, width = 18, height = 16, dpi = 300, limitsize = FALSE)
ggplot2::ggsave(summary_pdf, summary_plot, width = 18, height = 16, device = grDevices::cairo_pdf, limitsize = FALSE)
ggplot2::ggsave(diagnostic_png, diagnostic_plot, width = 13, height = 10, dpi = 300, limitsize = FALSE)

readr::write_csv(contrast_results, file.path(output_dir, "joint-model-contrasts.csv"))
readr::write_csv(interaction_tests, file.path(output_dir, "joint-model-stage-interaction-tests.csv"))
readr::write_csv(fixed_effects, file.path(output_dir, "joint-model-fixed-effects.csv"))
readr::write_csv(model_summary, file.path(output_dir, "joint-model-summary.csv"))
readr::write_csv(diagnostic_data, file.path(output_dir, "joint-model-diagnostic-data.csv"))
readr::write_csv(kenward_roger_tests, file.path(output_dir, "joint-model-kenward-roger-joint-tests.csv"))
readr::write_csv(residual_spread, file.path(output_dir, "joint-model-residual-spread-by-stage.csv"))
readr::write_csv(model_cell_counts, file.path(output_dir, "joint-model-cell-counts.csv"))

saveRDS(
  list(
    formula = joint_formula,
    models = model_results,
    contrasts = contrast_results,
    stage_interaction_tests = interaction_tests,
    kenward_roger_joint_tests = kenward_roger_tests,
    model_summary = model_summary
  ),
  file.path(output_dir, "phenology-joint-transition-model-test.rds")
)

print(tibble::as_tibble(model_summary), n = Inf)
print(
  interaction_tests %>%
    dplyr::select("species", "treatment_label", "chi_square", "df", "p_value") %>%
    tibble::as_tibble(),
  n = Inf
)
print(
  contrast_results %>%
    dplyr::select("species", "stage_label", "treatment_label", "estimate_days", "lower_95", "upper_95", "p_value") %>%
    tibble::as_tibble(),
  n = Inf
)

message("Saved joint-model test outputs in: ", normalizePath(output_dir, winslash = "/", mustWork = TRUE))
