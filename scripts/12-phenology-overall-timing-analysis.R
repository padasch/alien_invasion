#!/usr/bin/env Rscript

# Complete phenology analysis centred on overall spring leaf-out timing.
#
# Primary question:
#   Do reduced precipitation, Robinia presence, and mixed culture shift the
#   overall timing of stages 2--4?
#
# Primary model (one Gaussian LMM per focal species):
#   transition DOY ~ stage + precipitation + Robinia + culture + block
#                  + (1 | container) + (1 | tree)
#
# The additive treatment coefficients are common horizontal shifts of the
# stage-2--4 transition profile. A stage-by-treatment model checks whether that
# common-shift simplification is reasonable. Duration from stage 2 to stage 4
# is analysed separately. Finally, an associative piecewise path analysis uses
# a stage-centred overall timing index and container-level mean SWC during a
# fixed pre-leaf-out window that ends before the earliest focal transition.

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
species_labels <- c(fagus = "Fagus sylvatica", quercus = "Quercus ilex")
species_short <- c(fagus = "Fagus", quercus = "Quercus")
stages_keep <- 2:4
soil_type <- "both"
include_soil_treatment <- FALSE
spring_swc_start <- as.Date("2025-03-04")
spring_swc_end <- as.Date("2025-04-02")

args <- commandArgs(trailingOnly = TRUE)
bootstrap_arg <- grep("^--bootstrap=", args, value = TRUE)
n_boot <- if (length(bootstrap_arg)) {
  as.integer(sub("^--bootstrap=", "", bootstrap_arg[[1]]))
} else {
  1000L
}
if (is.na(n_boot) || n_boot < 100L) {
  stop("--bootstrap must be an integer of at least 100.", call. = FALSE)
}
cores_arg <- grep("^--cores=", args, value = TRUE)
n_cores <- if (length(cores_arg)) as.integer(sub("^--cores=", "", cores_arg[[1]])) else 4L
if (is.na(n_cores) || n_cores < 1L) stop("--cores must be a positive integer.", call. = FALSE)

alinv_set_analysis_context(
  scenario_label = "both soils (without soil as treatment)",
  soil_filter = soil_type,
  include_soil_treatment = include_soil_treatment,
  analysis_date = analysis_date,
  output_root = file.path(project_root, "output"),
  create_dirs = TRUE
)

output_dir <- file.path(project_root, "output", analysis_date, "phenology-overall-timing-analysis")
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

treatment_colors <- c(
  precipitation = "#1B9E77",
  robinia = "#D95F02",
  culture = "#7570B3"
)
precipitation_colors <- c(control = "#4F6674", drought = "#D65F5F")
component_colors <- c(
  "Direct" = "#3B6FB6",
  "Indirect via pre-leaf-out SWC" = "#D95F02",
  "Path-summed total" = "#222222",
  "Reduced-form total" = "#8A8A8A"
)
stage_colors <- c("Stage 2" = "#4F6674", "Stage 3" = "#D95F02", "Stage 4" = "#7570B3")

model_control <- lme4::lmerControl(
  optimizer = "bobyqa",
  optCtrl = list(maxfun = 200000),
  check.conv.singular = lme4::.makeCC(action = "ignore", tol = 1e-4)
)

add_block <- function(df) {
  df %>%
    dplyr::mutate(
      block = factor(sub("-.*$", "", as.character(.data$boxlabel)), levels = c("b1", "b2", "b3")),
      boxlabel = factor(.data$boxlabel),
      tree_id = factor(.data$tree_id),
      precipitation = factor(.data$precipitation, levels = alinv_factor_levels("precipitation")),
      robinia = factor(.data$robinia, levels = alinv_factor_levels("robinia")),
      culture = factor(.data$culture, levels = alinv_factor_levels("culture"))
    ) %>%
    droplevels()
}

transition_data <- stats::setNames(lapply(species_vec, function(species_i) {
  prepare_phenology_transition_data(
    species_keep = species_i,
    soil_type = soil_type,
    include_soil_treatment = include_soil_treatment,
    stages_keep = stages_keep
  ) %>%
    add_block()
}), species_vec)

primary_formula <- stats::as.formula(
  "doy ~ stage_label + precipitation + robinia + culture + block + (1 | boxlabel) + (1 | tree_id)"
)
interaction_formula <- stats::as.formula(
  "doy ~ stage_label * (precipitation + robinia + culture) + block + (1 | boxlabel) + (1 | tree_id)"
)

fit_transition_models <- function(df) {
  list(
    primary_reml = lme4::lmer(primary_formula, data = df, REML = TRUE, control = model_control),
    primary_ml = lme4::lmer(primary_formula, data = df, REML = FALSE, control = model_control),
    interaction_reml = lme4::lmer(interaction_formula, data = df, REML = TRUE, control = model_control),
    interaction_ml = lme4::lmer(interaction_formula, data = df, REML = FALSE, control = model_control)
  )
}

transition_models <- stats::setNames(lapply(transition_data, fit_transition_models), species_vec)

emmeans_contrast <- function(fit, species_i, spec_row, by_stage = FALSE) {
  effect_i <- spec_row$effect[[1]]
  expected_levels <- c(spec_row$baseline_level[[1]], spec_row$treatment_level[[1]])
  if (!identical(levels(fit@frame[[effect_i]]), expected_levels)) {
    stop("Unexpected factor level order for ", effect_i, ".", call. = FALSE)
  }

  specs_formula <- if (by_stage) {
    stats::as.formula(paste("~", effect_i, "| stage_label"))
  } else {
    stats::as.formula(paste("~", effect_i))
  }

  out <- emmeans::emmeans(
    fit,
    specs = specs_formula,
    weights = "equal",
    lmer.df = "kenward-roger"
  ) %>%
    emmeans::contrast(method = "revpairwise", by = if (by_stage) "stage_label" else NULL, adjust = "none") %>%
    summary(infer = c(TRUE, TRUE), level = 0.95, adjust = "none") %>%
    as.data.frame()

  if (!by_stage) out$stage_label <- "Stages 2–4 overall"

  out %>%
    dplyr::transmute(
      species = species_i,
      stage_label = as.character(.data$stage_label),
      effect = effect_i,
      treatment_label = spec_row$treatment_label[[1]],
      contrast_label = spec_row$contrast_label[[1]],
      estimate_days_raw = .data$estimate,
      lower_95_raw = .data$lower.CL,
      upper_95_raw = .data$upper.CL,
      estimate_oriented = -.data$estimate,
      lower_95_oriented = -.data$upper.CL,
      upper_95_oriented = -.data$lower.CL,
      se = .data$SE,
      df = .data$df,
      test_statistic = .data$t.ratio,
      p_value = .data$p.value,
      raw_interpretation = dplyr::case_when(
        .data$estimate_days_raw < 0 ~ "earlier",
        .data$estimate_days_raw > 0 ~ "later",
        TRUE ~ "no shift"
      ),
      plot_order = spec_row$plot_order[[1]]
    )
}

extract_treatment_contrasts <- function(fit, species_i, by_stage = FALSE) {
  dplyr::bind_rows(lapply(seq_len(nrow(treatment_specs)), function(i) {
    emmeans_contrast(fit, species_i, treatment_specs[i, ], by_stage = by_stage)
  }))
}

primary_effects <- dplyr::bind_rows(lapply(species_vec, function(species_i) {
  extract_treatment_contrasts(transition_models[[species_i]]$primary_reml, species_i, by_stage = FALSE)
}))

complete_case_data <- stats::setNames(lapply(species_vec, function(species_i) {
  df <- transition_data[[species_i]]
  complete_tree_ids <- df %>%
    dplyr::count(.data$tree_id, name = "n_stages") %>%
    dplyr::filter(.data$n_stages == length(stages_keep)) %>%
    dplyr::pull(.data$tree_id)
  df %>% dplyr::filter(.data$tree_id %in% complete_tree_ids) %>% droplevels()
}), species_vec)

complete_case_models <- stats::setNames(lapply(complete_case_data, function(df) {
  lme4::lmer(primary_formula, data = df, REML = TRUE, control = model_control)
}), species_vec)

complete_case_effects <- dplyr::bind_rows(lapply(species_vec, function(species_i) {
  extract_treatment_contrasts(complete_case_models[[species_i]], species_i, by_stage = FALSE)
}))

# A heteroscedastic-residual sensitivity model allows each stage to have a
# different residual SD. It is retained as a robustness check because the raw
# DOY residuals are discrete and visibly differ in spread among stages.
heteroscedastic_models <- stats::setNames(lapply(species_vec, function(species_i) {
  nlme::lme(
    fixed = doy ~ stage_label + precipitation + robinia + culture + block,
    random = ~1 | boxlabel/tree_id,
    weights = nlme::varIdent(form = ~1 | stage_label),
    data = transition_data[[species_i]],
    method = "REML",
    control = nlme::lmeControl(opt = "optim", maxIter = 200, msMaxIter = 200)
  )
}), species_vec)

heteroscedastic_effects <- dplyr::bind_rows(lapply(species_vec, function(species_i) {
  dplyr::bind_rows(lapply(seq_len(nrow(treatment_specs)), function(i) {
    spec_i <- treatment_specs[i, ]
    effect_i <- spec_i$effect[[1]]
    emm <- emmeans::emmeans(
      heteroscedastic_models[[species_i]],
      specs = stats::as.formula(paste("~", effect_i)),
      weights = "equal"
    ) %>%
      emmeans::contrast(method = "revpairwise", adjust = "none") %>%
      summary(infer = c(TRUE, TRUE), level = 0.95, adjust = "none") %>%
      as.data.frame()
    emm %>%
      dplyr::transmute(
        species = species_i,
        effect = effect_i,
        treatment_label = spec_i$treatment_label[[1]],
        estimate_days_raw = .data$estimate,
        lower_95_raw = .data$lower.CL,
        upper_95_raw = .data$upper.CL,
        estimate_oriented = -.data$estimate,
        lower_95_oriented = -.data$upper.CL,
        upper_95_oriented = -.data$lower.CL,
        se = .data$SE,
        df = .data$df,
        test_statistic = .data$t.ratio,
        p_value = .data$p.value,
        plot_order = spec_i$plot_order[[1]]
      )
  }))
}))

heteroscedastic_residual_multipliers <- dplyr::bind_rows(lapply(species_vec, function(species_i) {
  multipliers <- stats::coef(
    heteroscedastic_models[[species_i]]$modelStruct$varStruct,
    unconstrained = FALSE
  )
  tibble::tibble(
    species = species_i,
    stage_label = c("Stage 2", names(multipliers)),
    residual_sd_multiplier = c(1, unname(multipliers))
  )
}))

stage_effects <- dplyr::bind_rows(lapply(species_vec, function(species_i) {
  extract_treatment_contrasts(transition_models[[species_i]]$interaction_reml, species_i, by_stage = TRUE)
}))

test_stage_interactions <- function(species_i) {
  full_ml <- transition_models[[species_i]]$interaction_ml
  df <- transition_data[[species_i]]

  individual <- dplyr::bind_rows(lapply(seq_len(nrow(treatment_specs)), function(i) {
    spec_i <- treatment_specs[i, ]
    effect_i <- spec_i$effect[[1]]
    reduced_formula <- stats::update(
      stats::formula(full_ml),
      paste(". ~ . - stage_label:", effect_i)
    )
    reduced_fit <- lme4::lmer(reduced_formula, data = df, REML = FALSE, control = model_control)
    cmp <- stats::anova(reduced_fit, full_ml, refit = FALSE)
    tibble::tibble(
      species = species_i,
      test = paste0("stage × ", effect_i),
      effect = effect_i,
      treatment_label = spec_i$treatment_label[[1]],
      df = cmp$Df[[2]],
      chi_square = cmp$Chisq[[2]],
      p_value = cmp$`Pr(>Chisq)`[[2]],
      aic_reduced = cmp$AIC[[1]],
      aic_full = cmp$AIC[[2]]
    )
  }))

  overall_cmp <- stats::anova(
    transition_models[[species_i]]$primary_ml,
    transition_models[[species_i]]$interaction_ml,
    refit = FALSE
  )
  overall <- tibble::tibble(
    species = species_i,
    test = "all stage × treatment interactions",
    effect = "all",
    treatment_label = "All treatments",
    df = overall_cmp$Df[[2]],
    chi_square = overall_cmp$Chisq[[2]],
    p_value = overall_cmp$`Pr(>Chisq)`[[2]],
    aic_reduced = overall_cmp$AIC[[1]],
    aic_full = overall_cmp$AIC[[2]]
  )
  dplyr::bind_rows(overall, individual)
}

stage_interaction_tests <- dplyr::bind_rows(lapply(species_vec, test_stage_interactions))

make_duration_data <- function(df) {
  df %>%
    dplyr::filter(.data$stage %in% c(2L, 4L)) %>%
    dplyr::select(
      "tree_id", "boxlabel", "block", "precipitation", "robinia", "culture", "stage", "doy"
    ) %>%
    tidyr::pivot_wider(names_from = "stage", values_from = "doy", names_prefix = "stage_") %>%
    dplyr::filter(!is.na(.data$stage_2), !is.na(.data$stage_4)) %>%
    dplyr::mutate(duration_days = .data$stage_4 - .data$stage_2) %>%
    droplevels()
}

duration_data <- stats::setNames(lapply(transition_data, make_duration_data), species_vec)
duration_formula <- stats::as.formula(
  "duration_days ~ precipitation + robinia + culture + block + (1 | boxlabel)"
)
duration_models <- stats::setNames(lapply(duration_data, function(df) {
  lme4::lmer(duration_formula, data = df, REML = TRUE, control = model_control)
}), species_vec)

duration_effects <- dplyr::bind_rows(lapply(species_vec, function(species_i) {
  dplyr::bind_rows(lapply(seq_len(nrow(treatment_specs)), function(i) {
    spec_i <- treatment_specs[i, ]
    effect_i <- spec_i$effect[[1]]
    emm <- emmeans::emmeans(
      duration_models[[species_i]],
      specs = stats::as.formula(paste("~", effect_i)),
      weights = "equal",
      lmer.df = "kenward-roger"
    ) %>%
      emmeans::contrast(method = "revpairwise", adjust = "none") %>%
      summary(infer = c(TRUE, TRUE), level = 0.95, adjust = "none") %>%
      as.data.frame()
    emm %>%
      dplyr::transmute(
        species = species_i,
        effect = effect_i,
        treatment_label = spec_i$treatment_label[[1]],
        contrast_label = spec_i$contrast_label[[1]],
        estimate_days_raw = .data$estimate,
        lower_95_raw = .data$lower.CL,
        upper_95_raw = .data$upper.CL,
        estimate_oriented = -.data$estimate,
        lower_95_oriented = -.data$upper.CL,
        upper_95_oriented = -.data$lower.CL,
        se = .data$SE,
        df = .data$df,
        test_statistic = .data$t.ratio,
        p_value = .data$p.value,
        plot_order = spec_i$plot_order[[1]]
      )
  }))
}))

# Stage-centred timing index: a tree's mean deviation from the species-specific
# mean DOY of each observed stage. At least two of stages 2--4 are required.
make_timing_index <- function(df, species_i) {
  df %>%
    dplyr::group_by(.data$stage_label) %>%
    dplyr::mutate(stage_mean_doy = mean(.data$doy), stage_centered_doy = .data$doy - .data$stage_mean_doy) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(
      .data$tree_id, .data$boxlabel, .data$block,
      .data$precipitation, .data$robinia, .data$culture
    ) %>%
    dplyr::summarise(
      timing_index_days = mean(.data$stage_centered_doy),
      n_transitions = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::filter(.data$n_transitions >= 2L) %>%
    dplyr::mutate(
      timing_z = as.numeric(scale(.data$timing_index_days)),
      phenology_oriented_z = -.data$timing_z,
      species = species_i
    ) %>%
    droplevels()
}

timing_index_data <- stats::setNames(Map(make_timing_index, transition_data, species_vec), species_vec)

spring_swc <- get_data("box", "soilwater", swc_source = "measured") %>%
  dplyr::filter(.data$date >= spring_swc_start, .data$date <= spring_swc_end, !is.na(.data$swc)) %>%
  dplyr::group_by(
    .data$boxlabel, .data$precipitation, .data$robinia, .data$culture
  ) %>%
  dplyr::summarise(
    spring_swc = mean(.data$swc),
    spring_swc_sd = stats::sd(.data$swc),
    n_swc_dates = dplyr::n_distinct(.data$date),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    block = factor(sub("-.*$", "", .data$boxlabel), levels = c("b1", "b2", "b3")),
    precipitation = factor(.data$precipitation, levels = alinv_factor_levels("precipitation")),
    robinia = factor(.data$robinia, levels = alinv_factor_levels("robinia")),
    culture = factor(.data$culture, levels = alinv_factor_levels("culture"))
  )

prepare_sem_data <- function(species_i) {
  timing_i <- timing_index_data[[species_i]]
  boxes_i <- spring_swc %>%
    dplyr::filter(.data$boxlabel %in% as.character(timing_i$boxlabel)) %>%
    dplyr::mutate(
      swc_z = as.numeric(scale(.data$spring_swc)),
      boxlabel = factor(.data$boxlabel)
    ) %>%
    droplevels()
  trees_i <- timing_i %>%
    dplyr::mutate(boxlabel = as.character(.data$boxlabel)) %>%
    dplyr::left_join(
      boxes_i %>% dplyr::select("boxlabel", "swc_z", "spring_swc", "n_swc_dates"),
      by = "boxlabel"
    ) %>%
    dplyr::filter(!is.na(.data$swc_z)) %>%
    dplyr::mutate(boxlabel = factor(.data$boxlabel)) %>%
    droplevels()
  list(box = boxes_i, tree = trees_i)
}

sem_data <- stats::setNames(lapply(species_vec, prepare_sem_data), species_vec)

mediator_formula <- stats::as.formula("swc_z ~ precipitation + robinia + culture + block")
outcome_formula <- stats::as.formula(
  "phenology_oriented_z ~ swc_z + precipitation + robinia + culture + block + (1 | boxlabel)"
)
total_formula <- stats::as.formula(
  "phenology_oriented_z ~ precipitation + robinia + culture + block + (1 | boxlabel)"
)

coef_name_for <- function(spec_row) paste0(spec_row$effect[[1]], spec_row$treatment_level[[1]])

fit_sem_paths <- function(data_i, box_group = "boxlabel") {
  outcome_formula_i <- stats::as.formula(
    paste(
      "phenology_oriented_z ~ swc_z + precipitation + robinia + culture + block +",
      paste0("(1 | ", box_group, ")")
    )
  )
  total_formula_i <- stats::as.formula(
    paste(
      "phenology_oriented_z ~ precipitation + robinia + culture + block +",
      paste0("(1 | ", box_group, ")")
    )
  )
  mediator_fit <- stats::lm(mediator_formula, data = data_i$box)
  outcome_fit <- suppressWarnings(lme4::lmer(
    outcome_formula_i, data = data_i$tree, REML = TRUE, control = model_control
  ))
  total_fit <- suppressWarnings(lme4::lmer(
    total_formula_i, data = data_i$tree, REML = TRUE, control = model_control
  ))

  beta_a <- stats::coef(mediator_fit)
  beta_out <- lme4::fixef(outcome_fit)
  beta_total <- lme4::fixef(total_fit)
  b_path <- unname(beta_out[["swc_z"]])

  effects <- dplyr::bind_rows(lapply(seq_len(nrow(treatment_specs)), function(i) {
    spec_i <- treatment_specs[i, ]
    coef_i <- coef_name_for(spec_i)
    a_path <- unname(beta_a[[coef_i]])
    direct <- unname(beta_out[[coef_i]])
    reduced_total <- unname(beta_total[[coef_i]])
    tibble::tibble(
      effect = spec_i$effect[[1]],
      treatment_label = spec_i$treatment_label[[1]],
      plot_order = spec_i$plot_order[[1]],
      a_treatment_to_swc = a_path,
      b_swc_to_phenology = b_path,
      direct = direct,
      indirect = a_path * b_path,
      total_path = direct + a_path * b_path,
      total_reduced = reduced_total,
      total_discrepancy = direct + a_path * b_path - reduced_total
    )
  }))

  list(mediator = mediator_fit, outcome = outcome_fit, total = total_fit, effects = effects)
}

sem_fits <- stats::setNames(lapply(sem_data, fit_sem_paths), species_vec)

bootstrap_one_sem <- function(data_i, replicate_i, seed_i) {
  set.seed(seed_i)
  box_ids <- as.character(data_i$box$boxlabel)
  sampled_ids <- sample(box_ids, length(box_ids), replace = TRUE)
  sampled_copies <- seq_along(sampled_ids)

  boot_box <- dplyr::bind_rows(lapply(sampled_copies, function(j) {
    data_i$box %>%
      dplyr::filter(as.character(.data$boxlabel) == sampled_ids[[j]]) %>%
      dplyr::mutate(boot_box = factor(paste0("bootstrap-box-", j)))
  }))
  boot_tree <- dplyr::bind_rows(lapply(sampled_copies, function(j) {
    data_i$tree %>%
      dplyr::filter(as.character(.data$boxlabel) == sampled_ids[[j]]) %>%
      dplyr::mutate(boot_box = factor(paste0("bootstrap-box-", j)))
  }))

  if (nrow(boot_box) != length(box_ids) || !nrow(boot_tree)) return(NULL)
  if (any(vapply(c("precipitation", "robinia", "culture", "block"), function(x) {
    dplyr::n_distinct(boot_box[[x]]) < 2L
  }, logical(1)))) return(NULL)

  boot_data <- list(box = boot_box, tree = boot_tree)
  fit <- tryCatch(
    suppressWarnings(fit_sem_paths(boot_data, box_group = "boot_box")),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  fit$effects %>% dplyr::mutate(replicate = replicate_i)
}

bootstrap_sem <- function(species_i, n_boot_i, seed = 20260818L) {
  data_i <- sem_data[[species_i]]
  set.seed(seed + match(species_i, species_vec) * 10000L)
  seeds <- sample.int(.Machine$integer.max, n_boot_i)
  cores <- if (.Platform$OS.type == "windows") 1L else n_cores
  results <- parallel::mclapply(
    seq_len(n_boot_i),
    function(i) bootstrap_one_sem(data_i, i, seeds[[i]]),
    mc.cores = cores,
    mc.preschedule = TRUE
  )
  dplyr::bind_rows(results) %>% dplyr::mutate(species = species_i, .before = 1)
}

message("Running ", n_boot, " container-cluster bootstrap replicates per species for the SEM...")
sem_bootstrap <- dplyr::bind_rows(lapply(species_vec, bootstrap_sem, n_boot_i = n_boot))

point_sem_effects <- dplyr::bind_rows(lapply(species_vec, function(species_i) {
  sem_fits[[species_i]]$effects %>% dplyr::mutate(species = species_i, .before = 1)
}))

summarise_boot_metric <- function(metric) {
  sem_bootstrap %>%
    dplyr::group_by(.data$species, .data$effect, .data$treatment_label, .data$plot_order) %>%
    dplyr::summarise(
      lower_95 = stats::quantile(.data[[metric]], 0.025, na.rm = TRUE, names = FALSE),
      upper_95 = stats::quantile(.data[[metric]], 0.975, na.rm = TRUE, names = FALSE),
      p_boot = min(1, 2 * min(
        (sum(.data[[metric]] <= 0, na.rm = TRUE) + 1) / (sum(!is.na(.data[[metric]])) + 1),
        (sum(.data[[metric]] >= 0, na.rm = TRUE) + 1) / (sum(!is.na(.data[[metric]])) + 1)
      )),
      n_boot_success = sum(!is.na(.data[[metric]])),
      .groups = "drop"
    ) %>%
    dplyr::left_join(
      point_sem_effects %>%
        dplyr::select("species", "effect", "treatment_label", estimate = dplyr::all_of(metric)),
      by = c("species", "effect", "treatment_label")
    )
}

metric_labels <- c(
  direct = "Direct",
  indirect = "Indirect via pre-leaf-out SWC",
  total_path = "Path-summed total",
  total_reduced = "Reduced-form total"
)

sem_effects <- dplyr::bind_rows(lapply(names(metric_labels), function(metric_i) {
  summarise_boot_metric(metric_i) %>%
    dplyr::mutate(component = unname(metric_labels[[metric_i]]), metric = metric_i)
}))

sem_paths <- dplyr::bind_rows(
  summarise_boot_metric("a_treatment_to_swc") %>%
    dplyr::mutate(path = "Treatment → pre-leaf-out SWC", metric = "a_treatment_to_swc"),
  summarise_boot_metric("direct") %>%
    dplyr::mutate(path = "Treatment → phenology (direct)", metric = "direct")
)

sem_b_path <- sem_bootstrap %>%
  dplyr::distinct(.data$species, .data$replicate, .data$b_swc_to_phenology) %>%
  dplyr::group_by(.data$species) %>%
  dplyr::summarise(
    lower_95 = stats::quantile(.data$b_swc_to_phenology, 0.025, na.rm = TRUE, names = FALSE),
    upper_95 = stats::quantile(.data$b_swc_to_phenology, 0.975, na.rm = TRUE, names = FALSE),
    p_boot = min(1, 2 * min(
      (sum(.data$b_swc_to_phenology <= 0, na.rm = TRUE) + 1) / (sum(!is.na(.data$b_swc_to_phenology)) + 1),
      (sum(.data$b_swc_to_phenology >= 0, na.rm = TRUE) + 1) / (sum(!is.na(.data$b_swc_to_phenology)) + 1)
    )),
    n_boot_success = sum(!is.na(.data$b_swc_to_phenology)),
    .groups = "drop"
  ) %>%
  dplyr::left_join(
    point_sem_effects %>%
      dplyr::distinct(.data$species, estimate = .data$b_swc_to_phenology),
    by = "species"
  ) %>%
  dplyr::mutate(
    effect = "swc",
    treatment_label = "Pre-leaf-out SWC",
    plot_order = 4L,
    path = "Pre-leaf-out SWC → phenology",
    metric = "b_swc_to_phenology"
  )

sem_paths <- dplyr::bind_rows(sem_paths, sem_b_path)

model_status <- function(fit, species_i, model_name, data, n_trees = NA_integer_) {
  convergence_messages <- fit@optinfo$conv$lme4$messages
  if (is.null(convergence_messages)) convergence_messages <- "none"
  tibble::tibble(
    species = species_i,
    model = model_name,
    formula = paste(deparse(stats::formula(fit)), collapse = " "),
    n_observations = stats::nobs(fit),
    n_trees = n_trees,
    n_containers = dplyr::n_distinct(data$boxlabel),
    singular = lme4::isSingular(fit, tol = 1e-4),
    convergence_messages = paste(convergence_messages, collapse = "; "),
    sigma = stats::sigma(fit)
  )
}

model_summary <- dplyr::bind_rows(lapply(species_vec, function(species_i) {
  dplyr::bind_rows(
    model_status(
      transition_models[[species_i]]$primary_reml, species_i, "primary common-shift LMM",
      transition_data[[species_i]], dplyr::n_distinct(transition_data[[species_i]]$tree_id)
    ),
    model_status(
      transition_models[[species_i]]$interaction_reml, species_i, "stage-interaction sensitivity LMM",
      transition_data[[species_i]], dplyr::n_distinct(transition_data[[species_i]]$tree_id)
    ),
    model_status(
      complete_case_models[[species_i]], species_i, "all-three-transitions sensitivity LMM",
      complete_case_data[[species_i]], dplyr::n_distinct(complete_case_data[[species_i]]$tree_id)
    ),
    model_status(
      duration_models[[species_i]], species_i, "stage 2–4 duration LMM",
      duration_data[[species_i]], nrow(duration_data[[species_i]])
    ),
    model_status(
      sem_fits[[species_i]]$outcome, species_i, "SEM phenology outcome LMM",
      sem_data[[species_i]]$tree, nrow(sem_data[[species_i]]$tree)
    ),
    model_status(
      sem_fits[[species_i]]$total, species_i, "SEM reduced-form total LMM",
      sem_data[[species_i]]$tree, nrow(sem_data[[species_i]]$tree)
    )
  )
}))

phenology_all <- get_data("tree", "phenology_transitions") %>%
  alinv_filter_by_soil(soil_filter = soil_type) %>%
  dplyr::filter(.data$species %in% species_vec, .data$stage %in% 1:4)

transition_missingness <- phenology_all %>%
  dplyr::mutate(observed = !is.na(.data$doy) & !is.na(.data$stage_date)) %>%
  dplyr::group_by(
    .data$species, .data$stage, .data$precipitation, .data$robinia, .data$culture
  ) %>%
  dplyr::summarise(
    n_expected = dplyr::n(),
    n_observed = sum(.data$observed),
    proportion_observed = mean(.data$observed),
    .groups = "drop"
  )

missingness_model_data <- phenology_all %>%
  dplyr::filter(.data$stage %in% 2:3) %>%
  dplyr::mutate(
    observed = !is.na(.data$doy) & !is.na(.data$stage_date),
    block = factor(sub("-.*$", "", .data$boxlabel), levels = c("b1", "b2", "b3")),
    precipitation = factor(.data$precipitation, levels = alinv_factor_levels("precipitation")),
    robinia = factor(.data$robinia, levels = alinv_factor_levels("robinia")),
    culture = factor(.data$culture, levels = alinv_factor_levels("culture"))
  )

missingness_models <- list()
missingness_effects <- dplyr::bind_rows(lapply(species_vec, function(species_i) {
  dplyr::bind_rows(lapply(2:3, function(stage_i) {
    df <- missingness_model_data %>%
      dplyr::filter(.data$species == species_i, .data$stage == stage_i) %>%
      droplevels()
    fit <- stats::glm(
      observed ~ precipitation + robinia + culture + block,
      data = df,
      family = stats::binomial()
    )
    missingness_models[[paste(species_i, stage_i, sep = "-")]] <<- fit
    vcov_cluster <- sandwich::vcovCL(fit, cluster = df$boxlabel, type = "HC1")
    beta <- stats::coef(fit)

    dplyr::bind_rows(lapply(seq_len(nrow(treatment_specs)), function(i) {
      spec_i <- treatment_specs[i, ]
      coef_i <- coef_name_for(spec_i)
      estimate_i <- unname(beta[[coef_i]])
      se_i <- sqrt(unname(vcov_cluster[coef_i, coef_i]))
      tibble::tibble(
        species = species_i,
        stage = stage_i,
        effect = spec_i$effect[[1]],
        treatment_label = spec_i$treatment_label[[1]],
        estimate_log_odds = estimate_i,
        se_cluster = se_i,
        lower_95_log_odds = estimate_i - 1.96 * se_i,
        upper_95_log_odds = estimate_i + 1.96 * se_i,
        odds_ratio = exp(estimate_i),
        lower_95_odds_ratio = exp(estimate_i - 1.96 * se_i),
        upper_95_odds_ratio = exp(estimate_i + 1.96 * se_i),
        z_value = estimate_i / se_i,
        p_value = 2 * stats::pnorm(-abs(estimate_i / se_i)),
        n_expected = nrow(df),
        n_observed = sum(df$observed),
        n_containers = dplyr::n_distinct(df$boxlabel),
        plot_order = spec_i$plot_order[[1]]
      )
    }))
  }))
}))

descriptive_data <- phenology_all %>%
  dplyr::filter(!is.na(.data$doy), !is.na(.data$stage_date)) %>%
  dplyr::mutate(
    species = factor(.data$species, levels = species_vec, labels = unname(species_short[species_vec])),
    robinia = factor(
      .data$robinia,
      levels = alinv_factor_levels("robinia"),
      labels = c("Without Robinia", "With Robinia")
    ),
    precipitation = factor(.data$precipitation, levels = alinv_factor_levels("precipitation")),
    culture = factor(.data$culture, levels = alinv_factor_levels("culture")),
    stage = as.integer(.data$stage)
  ) %>%
  dplyr::group_by(.data$species, .data$robinia, .data$precipitation, .data$culture, .data$stage) %>%
  dplyr::summarise(
    mean_doy = mean(.data$doy),
    se_doy = stats::sd(.data$doy) / sqrt(dplyr::n()),
    n = dplyr::n(),
    .groups = "drop"
  )

plot_timeseries <- ggplot2::ggplot(
  descriptive_data,
  ggplot2::aes(
    x = .data$mean_doy,
    y = .data$stage,
    color = .data$precipitation,
    linetype = .data$culture,
    group = interaction(.data$precipitation, .data$culture)
  )
) +
  ggplot2::geom_segment(
    ggplot2::aes(x = .data$mean_doy - .data$se_doy, xend = .data$mean_doy + .data$se_doy, yend = .data$stage),
    linewidth = 0.65, alpha = 0.5, lineend = "round"
  ) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_point(size = 2.1) +
  ggplot2::facet_grid(.data$species ~ .data$robinia) +
  ggplot2::scale_color_manual(
    values = precipitation_colors,
    labels = c(control = "Control", drought = "Reduced precipitation")
  ) +
  ggplot2::scale_linetype_manual(
    values = c(mono = "solid", mixed = "dotted"),
    labels = c(mono = "Monoculture", mixed = "Mixed culture")
  ) +
  ggplot2::scale_y_continuous(breaks = 1:4, labels = paste("Stage", 1:4), limits = c(0.8, 4.2)) +
  ggplot2::labs(
    title = "Observed spring phenology progression",
    subtitle = "Treatment-group mean transition DOY ± SE. Stage 1 is shown descriptively; inference uses stages 2–4.",
    x = "Day of year (DOY)", y = "Phenological stage", color = "Precipitation", linetype = "Culture"
  ) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(
    panel.grid.major = ggplot2::element_line(color = "grey90", linewidth = 0.25),
    panel.grid.minor = ggplot2::element_blank(),
    strip.background = ggplot2::element_rect(fill = "grey12", color = "grey12"),
    strip.text = ggplot2::element_text(color = "white", face = "bold"),
    legend.position = "bottom", plot.title = ggplot2::element_text(face = "bold")
  )

prep_effect_plot <- function(df) {
  df %>%
    dplyr::mutate(
      species = factor(.data$species, levels = species_vec, labels = unname(species_short[species_vec])),
      treatment_label = factor(.data$treatment_label, levels = treatment_specs$treatment_label),
      significant = .data$p_value < 0.05
    )
}

primary_plot_data <- prep_effect_plot(primary_effects)
plot_primary <- ggplot2::ggplot(
  primary_plot_data,
  ggplot2::aes(x = .data$estimate_oriented, y = .data$treatment_label, color = .data$effect)
) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey45") +
  ggplot2::geom_segment(
    ggplot2::aes(x = .data$lower_95_oriented, xend = .data$upper_95_oriented, yend = .data$treatment_label),
    linewidth = 1.0, lineend = "round"
  ) +
  ggplot2::geom_point(ggplot2::aes(shape = .data$significant), size = 3.0, fill = "white") +
  ggplot2::facet_wrap(~species, nrow = 1) +
  ggplot2::scale_color_manual(values = treatment_colors, guide = "none") +
  ggplot2::scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), guide = "none") +
  ggplot2::labs(
    title = "Overall shift of the stage-2–4 transition profile",
    subtitle = "One additive Gaussian LMM per species; points are treatment contrasts with 95% Kenward–Roger intervals.",
    x = "Oriented timing effect (days): negative = delayed, positive = earlier", y = NULL
  ) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(
    panel.grid.major.y = ggplot2::element_line(color = "grey92"),
    strip.background = ggplot2::element_rect(fill = "grey12", color = "grey12"),
    strip.text = ggplot2::element_text(color = "white", face = "bold"),
    plot.title = ggplot2::element_text(face = "bold")
  )

stage_plot_data <- prep_effect_plot(stage_effects) %>%
  dplyr::mutate(stage_label = factor(.data$stage_label, levels = paste("Stage", stages_keep)))
plot_stage_check <- ggplot2::ggplot(
  stage_plot_data,
  ggplot2::aes(x = .data$estimate_oriented, y = .data$treatment_label, color = .data$effect)
) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey45") +
  ggplot2::geom_segment(
    ggplot2::aes(x = .data$lower_95_oriented, xend = .data$upper_95_oriented, yend = .data$treatment_label),
    linewidth = 0.9, lineend = "round"
  ) +
  ggplot2::geom_point(ggplot2::aes(shape = .data$significant), size = 2.7, fill = "white") +
  ggplot2::facet_grid(.data$species ~ .data$stage_label) +
  ggplot2::scale_color_manual(values = treatment_colors, guide = "none") +
  ggplot2::scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), guide = "none") +
  ggplot2::labs(
    title = "Sensitivity check: treatment shifts estimated separately at each stage",
    subtitle = "Derived from one stage × treatment LMM per species; this tests the common-shift assumption rather than replacing the primary model.",
    x = "Oriented timing effect (days): negative = delayed, positive = earlier", y = NULL
  ) +
  ggplot2::theme_classic(base_size = 10.5) +
  ggplot2::theme(
    panel.grid.major.y = ggplot2::element_line(color = "grey92"),
    strip.background = ggplot2::element_rect(fill = "grey12", color = "grey12"),
    strip.text = ggplot2::element_text(color = "white", face = "bold"),
    plot.title = ggplot2::element_text(face = "bold")
  )

interaction_plot_data <- stage_interaction_tests %>%
  dplyr::filter(.data$effect != "all") %>%
  dplyr::mutate(
    species = factor(.data$species, levels = species_vec, labels = unname(species_short[species_vec])),
    treatment_label = factor(.data$treatment_label, levels = treatment_specs$treatment_label),
    evidence = factor(ifelse(.data$p_value < 0.05, "P < 0.05", "P ≥ 0.05"), levels = c("P ≥ 0.05", "P < 0.05")),
    label = sprintf("χ²(%d) = %.2f\nP = %s", .data$df, .data$chi_square, format.pval(.data$p_value, digits = 2, eps = 0.001))
  )
plot_interactions <- ggplot2::ggplot(
  interaction_plot_data,
  ggplot2::aes(x = .data$treatment_label, y = .data$species, fill = .data$evidence)
) +
  ggplot2::geom_tile(color = "grey70") +
  ggplot2::geom_text(ggplot2::aes(label = .data$label), size = 3.2, lineheight = 0.95) +
  ggplot2::scale_fill_manual(values = c("P ≥ 0.05" = "white", "P < 0.05" = "#E6B44A"), drop = FALSE) +
  ggplot2::labs(
    title = "Stage × treatment likelihood-ratio tests",
    subtitle = "P < 0.05 indicates that a treatment shift is not constant across stages.",
    x = NULL, y = NULL, fill = NULL
  ) +
  ggplot2::theme_minimal(base_size = 10.5) +
  ggplot2::theme(panel.grid = ggplot2::element_blank(), legend.position = "bottom", plot.title = ggplot2::element_text(face = "bold"))

duration_plot_data <- duration_effects %>%
  dplyr::mutate(
    species = factor(.data$species, levels = species_vec, labels = unname(species_short[species_vec])),
    treatment_label = factor(.data$treatment_label, levels = treatment_specs$treatment_label),
    significant = .data$p_value < 0.05
  )
plot_duration <- ggplot2::ggplot(
  duration_plot_data,
  ggplot2::aes(x = .data$estimate_oriented, y = .data$treatment_label, color = .data$effect)
) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey45") +
  ggplot2::geom_segment(
    ggplot2::aes(x = .data$lower_95_oriented, xend = .data$upper_95_oriented, yend = .data$treatment_label),
    linewidth = 1.0, lineend = "round"
  ) +
  ggplot2::geom_point(ggplot2::aes(shape = .data$significant), size = 3.0, fill = "white") +
  ggplot2::facet_wrap(~species, nrow = 1) +
  ggplot2::scale_color_manual(values = treatment_colors, guide = "none") +
  ggplot2::scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), guide = "none") +
  ggplot2::labs(
    title = "Change in the duration from stage 2 to stage 4",
    subtitle = "Secondary LMM among trees with both transitions observed; 95% Kenward–Roger intervals.",
    x = "Oriented duration effect (days): negative = slower/longer, positive = faster/shorter", y = NULL
  ) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(
    panel.grid.major.y = ggplot2::element_line(color = "grey92"),
    strip.background = ggplot2::element_rect(fill = "grey12", color = "grey12"),
    strip.text = ggplot2::element_text(color = "white", face = "bold"),
    plot.title = ggplot2::element_text(face = "bold")
  )

sem_plot_data <- sem_effects %>%
  dplyr::filter(.data$component %in% c("Direct", "Indirect via pre-leaf-out SWC", "Path-summed total")) %>%
  dplyr::mutate(
    species = factor(.data$species, levels = species_vec, labels = unname(species_short[species_vec])),
    treatment_label = factor(.data$treatment_label, levels = treatment_specs$treatment_label),
    component = factor(.data$component, levels = c("Direct", "Indirect via pre-leaf-out SWC", "Path-summed total"))
  )
plot_sem_effects <- ggplot2::ggplot(
  sem_plot_data,
  ggplot2::aes(x = .data$estimate, y = .data$treatment_label, color = .data$component)
) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey45") +
  ggplot2::geom_segment(
    ggplot2::aes(x = .data$lower_95, xend = .data$upper_95, yend = .data$treatment_label),
    position = ggplot2::position_dodge(width = 0.55), linewidth = 0.85, lineend = "round"
  ) +
  ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.55), size = 2.6) +
  ggplot2::facet_wrap(~species, nrow = 1) +
  ggplot2::scale_color_manual(values = component_colors, drop = FALSE) +
  ggplot2::labs(
    title = "Associative SEM decomposition of overall phenology timing",
    subtitle = paste0(
      "Standardized effects with container-cluster bootstrap intervals (", n_boot,
      " replicates). Pre-leaf-out SWC is the container mean from 4 March to 2 April 2025."
    ),
    x = "Oriented standardized effect: negative = delayed, positive = earlier",
    y = NULL, color = "Effect component"
  ) +
  ggplot2::theme_classic(base_size = 10.5) +
  ggplot2::theme(
    panel.grid.major.y = ggplot2::element_line(color = "grey92"),
    strip.background = ggplot2::element_rect(fill = "grey12", color = "grey12"),
    strip.text = ggplot2::element_text(color = "white", face = "bold"),
    legend.position = "bottom", plot.title = ggplot2::element_text(face = "bold")
  )

sem_heatmap_data <- dplyr::bind_rows(
  sem_effects %>%
    dplyr::filter(.data$component %in% c("Direct", "Indirect via pre-leaf-out SWC", "Path-summed total")) %>%
    dplyr::transmute(
      species = .data$species,
      treatment_label = .data$treatment_label,
      component = dplyr::recode(
        .data$component,
        "Indirect via pre-leaf-out SWC" = "Indirect via SWC",
        "Path-summed total" = "Total"
      ),
      estimate = .data$estimate,
      lower_95 = .data$lower_95,
      upper_95 = .data$upper_95,
      p_boot = .data$p_boot
    ),
  sem_paths %>%
    dplyr::filter(.data$path == "Treatment → pre-leaf-out SWC") %>%
    dplyr::transmute(
      species = .data$species,
      treatment_label = .data$treatment_label,
      component = "Treatment → SWC",
      estimate = .data$estimate,
      lower_95 = .data$lower_95,
      upper_95 = .data$upper_95,
      p_boot = .data$p_boot
    )
) %>%
  dplyr::mutate(
    species = factor(.data$species, levels = species_vec, labels = unname(species_short[species_vec])),
    treatment_label = factor(.data$treatment_label, levels = treatment_specs$treatment_label),
    component = factor(
      .data$component,
      levels = c("Direct", "Indirect via SWC", "Total", "Treatment → SWC")
    ),
    interval_excludes_zero = .data$lower_95 > 0 | .data$upper_95 < 0,
    cell_label = paste0(sprintf("%.2f", .data$estimate), ifelse(.data$interval_excludes_zero, "*", ""))
  )

sem_heatmap_limit <- max(abs(sem_heatmap_data$estimate), na.rm = TRUE)
sem_heatmap_limit <- max(0.1, ceiling(sem_heatmap_limit * 10) / 10)
sem_heatmap_data <- sem_heatmap_data %>%
  dplyr::mutate(
    text_color = ifelse(abs(.data$estimate) >= 0.62 * sem_heatmap_limit, "white", "grey15"),
    fontface = ifelse(.data$interval_excludes_zero, "bold", "plain")
  )

plot_sem_heatmap <- ggplot2::ggplot(
  sem_heatmap_data,
  ggplot2::aes(x = .data$component, y = .data$treatment_label, fill = .data$estimate)
) +
  ggplot2::geom_tile(color = "grey82", linewidth = 0.45) +
  ggplot2::geom_text(
    ggplot2::aes(label = .data$cell_label, color = .data$text_color, fontface = .data$fontface),
    size = 3.4
  ) +
  ggplot2::facet_grid(.data$species ~ ., switch = "y") +
  ggplot2::scale_fill_gradient2(
    low = "#C94F50", mid = "white", high = "#4C78A8", midpoint = 0,
    limits = c(-sem_heatmap_limit, sem_heatmap_limit), name = "Std. effect"
  ) +
  ggplot2::scale_color_identity() +
  ggplot2::labs(
    title = "Associative SEM heatmap",
    subtitle = paste(
      "Phenology-directed effects are oriented so negative values indicate delayed timing;",
      "negative treatment → SWC effects indicate lower pre-leaf-out SWC."
    ),
    x = NULL, y = NULL,
    caption = "Cell values are standardized effects; * indicates a 95% container-bootstrap interval excluding zero."
  ) +
  ggplot2::theme_minimal(base_size = 10.5) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    strip.placement = "outside",
    strip.background = ggplot2::element_rect(fill = "grey12", color = "grey12"),
    strip.text.y.left = ggplot2::element_text(color = "white", face = "bold", angle = 90),
    axis.text.x = ggplot2::element_text(face = "bold", angle = 20, hjust = 1),
    axis.text.y = ggplot2::element_text(color = "grey20"),
    legend.position = "right",
    plot.title = ggplot2::element_text(face = "bold"),
    plot.caption = ggplot2::element_text(color = "grey35")
  )

sem_path_plot_data <- sem_paths %>%
  dplyr::mutate(
    species = factor(.data$species, levels = species_vec, labels = unname(species_short[species_vec])),
    path = factor(
      .data$path,
      levels = c("Treatment → pre-leaf-out SWC", "Pre-leaf-out SWC → phenology", "Treatment → phenology (direct)")
    ),
    label = dplyr::case_when(
      .data$effect == "precipitation" ~ "Reduced precipitation",
      .data$effect == "robinia" ~ "With Robinia",
      .data$effect == "culture" ~ "Mixed culture",
      TRUE ~ "Pre-leaf-out SWC"
    )
  )
plot_sem_paths <- ggplot2::ggplot(
  sem_path_plot_data,
  ggplot2::aes(x = .data$estimate, y = .data$label, color = .data$effect)
) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey45") +
  ggplot2::geom_segment(
    ggplot2::aes(x = .data$lower_95, xend = .data$upper_95, yend = .data$label),
    linewidth = 0.9, lineend = "round"
  ) +
  ggplot2::geom_point(size = 2.7) +
  ggplot2::facet_wrap(~species + path, ncol = 3, scales = "free_y") +
  ggplot2::scale_color_manual(values = c(treatment_colors, swc = "#4F6674"), guide = "none") +
  ggplot2::labs(
    title = "Constituent standardized paths in the phenology SEM",
    subtitle = "Intervals are from the same container-cluster bootstrap. Phenology-directed paths use the oriented timing index.",
    x = "Standardized path coefficient", y = NULL
  ) +
  ggplot2::theme_classic(base_size = 10) +
  ggplot2::theme(
    panel.grid.major.y = ggplot2::element_line(color = "grey92"),
    strip.background = ggplot2::element_rect(fill = "grey12", color = "grey12"),
    strip.text = ggplot2::element_text(color = "white", face = "bold"),
    plot.title = ggplot2::element_text(face = "bold")
  )

diagnostic_data <- dplyr::bind_rows(lapply(species_vec, function(species_i) {
  fit <- transition_models[[species_i]]$primary_reml
  tibble::tibble(
    species = species_i,
    stage_label = transition_data[[species_i]]$stage_label,
    fitted = stats::fitted(fit),
    residual = stats::residuals(fit, type = "pearson")
  )
}))

plot_diagnostics <- (
  ggplot2::ggplot(diagnostic_data, ggplot2::aes(x = .data$fitted, y = .data$residual, color = .data$stage_label)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    ggplot2::geom_point(alpha = 0.3, size = 1.0) +
    ggplot2::geom_smooth(method = "loess", se = FALSE, linewidth = 0.7) +
    ggplot2::facet_wrap(~species, scales = "free_x", labeller = ggplot2::as_labeller(species_short)) +
    ggplot2::scale_color_manual(values = stage_colors) +
    ggplot2::labs(title = "Residuals versus fitted values", x = "Fitted transition DOY", y = "Pearson residual", color = "Stage") +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(legend.position = "bottom", plot.title = ggplot2::element_text(face = "bold"))
) / (
  ggplot2::ggplot(diagnostic_data, ggplot2::aes(sample = .data$residual, color = .data$stage_label)) +
    ggplot2::stat_qq(alpha = 0.3, size = 1.0) +
    ggplot2::stat_qq_line(color = "grey25", linewidth = 0.55) +
    ggplot2::facet_wrap(~species, labeller = ggplot2::as_labeller(species_short)) +
    ggplot2::scale_color_manual(values = stage_colors) +
    ggplot2::labs(title = "Normal Q–Q plots", x = "Theoretical quantile", y = "Observed quantile", color = "Stage") +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(legend.position = "bottom", plot.title = ggplot2::element_text(face = "bold"))
) +
  patchwork::plot_annotation(
    title = "Primary common-shift LMM diagnostics",
    subtitle = "Models are fitted separately by species; colour identifies the observed stage."
  )

missingness_plot_data <- missingness_effects %>%
  dplyr::mutate(
    species = factor(.data$species, levels = species_vec, labels = unname(species_short[species_vec])),
    stage_label = factor(paste("Stage", .data$stage), levels = c("Stage 2", "Stage 3")),
    treatment_label = factor(.data$treatment_label, levels = treatment_specs$treatment_label),
    significant = .data$p_value < 0.05
  )
plot_missingness <- ggplot2::ggplot(
  missingness_plot_data,
  ggplot2::aes(x = .data$estimate_log_odds, y = .data$treatment_label, color = .data$effect)
) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey45") +
  ggplot2::geom_segment(
    ggplot2::aes(x = .data$lower_95_log_odds, xend = .data$upper_95_log_odds, yend = .data$treatment_label),
    linewidth = 0.9, lineend = "round"
  ) +
  ggplot2::geom_point(ggplot2::aes(shape = .data$significant), size = 2.7) +
  ggplot2::facet_grid(.data$species ~ .data$stage_label) +
  ggplot2::scale_color_manual(values = treatment_colors, guide = "none") +
  ggplot2::scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), guide = "none") +
  ggplot2::labs(
    title = "Were transition dates equally likely to be observed?",
    subtitle = "Binomial models with container-clustered SEs. Negative estimates indicate lower observation probability under treatment.",
    x = "Treatment effect on log-odds that the transition date was observed", y = NULL
  ) +
  ggplot2::theme_classic(base_size = 10.5) +
  ggplot2::theme(
    panel.grid.major.y = ggplot2::element_line(color = "grey92"),
    strip.background = ggplot2::element_rect(fill = "grey12", color = "grey12"),
    strip.text = ggplot2::element_text(color = "white", face = "bold"),
    plot.title = ggplot2::element_text(face = "bold")
  )

complete_case_plot_data <- dplyr::bind_rows(
  primary_effects %>% dplyr::mutate(dataset = "All available transitions"),
  complete_case_effects %>% dplyr::mutate(dataset = "Trees with all three transitions")
) %>%
  dplyr::mutate(
    species = factor(.data$species, levels = species_vec, labels = unname(species_short[species_vec])),
    treatment_label = factor(.data$treatment_label, levels = treatment_specs$treatment_label),
    dataset = factor(.data$dataset, levels = c("All available transitions", "Trees with all three transitions"))
  )
plot_complete_case <- ggplot2::ggplot(
  complete_case_plot_data,
  ggplot2::aes(x = .data$estimate_oriented, y = .data$treatment_label, color = .data$dataset)
) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey45") +
  ggplot2::geom_segment(
    ggplot2::aes(x = .data$lower_95_oriented, xend = .data$upper_95_oriented, yend = .data$treatment_label),
    position = ggplot2::position_dodge(width = 0.48), linewidth = 0.85, lineend = "round"
  ) +
  ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.48), size = 2.6) +
  ggplot2::facet_wrap(~species, nrow = 1) +
  ggplot2::scale_color_manual(values = c("All available transitions" = "#222222", "Trees with all three transitions" = "#C44E52")) +
  ggplot2::labs(
    title = "Sensitivity of the overall timing estimate to incomplete transition records",
    subtitle = "Both analyses use the same common-shift LMM; the strict subset requires stages 2, 3, and 4 for every included tree.",
    x = "Oriented timing effect (days): negative = delayed, positive = earlier", y = NULL, color = "Dataset"
  ) +
  ggplot2::theme_classic(base_size = 10.5) +
  ggplot2::theme(
    panel.grid.major.y = ggplot2::element_line(color = "grey92"),
    strip.background = ggplot2::element_rect(fill = "grey12", color = "grey12"),
    strip.text = ggplot2::element_text(color = "white", face = "bold"),
    legend.position = "bottom", plot.title = ggplot2::element_text(face = "bold")
  )

plot_missingness_sensitivity <- plot_missingness / plot_complete_case +
  patchwork::plot_layout(heights = c(1.05, 0.95)) +
  patchwork::plot_annotation(
    title = "Transition-date completeness and complete-case sensitivity",
    tag_levels = "A"
  )

summary_plot <- (
  plot_timeseries /
    plot_primary /
    plot_duration /
    plot_sem_effects /
    plot_sem_heatmap
) +
  patchwork::plot_layout(heights = c(1.35, 0.8, 0.8, 1.0, 1.0), guides = "keep") +
  patchwork::plot_annotation(
    title = "Phenology: overall timing, progression duration, and associative SWC pathways",
    subtitle = paste(
      "Primary inference is the common shift of stages 2–4; duration and stage-specific estimates are sensitivity analyses.",
      "For cross-response consistency, negative oriented effects mean delayed/slower phenology."
    ),
    tag_levels = "A",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 16),
      plot.subtitle = ggplot2::element_text(size = 10.5)
    )
  ) &
  ggplot2::theme(plot.margin = ggplot2::margin(t = 7, r = 9, b = 7, l = 9))

save_plot_pair <- function(filename, plot, width, height) {
  ggplot2::ggsave(file.path(output_dir, paste0(filename, ".png")), plot, width = width, height = height, dpi = 300, limitsize = FALSE)
  ggplot2::ggsave(file.path(output_dir, paste0(filename, ".pdf")), plot, width = width, height = height, device = grDevices::cairo_pdf, limitsize = FALSE)
}

save_plot_pair("01-phenology-transition-progression", plot_timeseries, 14, 8)
save_plot_pair("02-primary-common-shift-lmm", plot_primary, 11, 5.5)
save_plot_pair("03-stage-interaction-sensitivity", plot_stage_check / plot_interactions + patchwork::plot_layout(heights = c(1.8, 0.8)), 14, 10)
save_plot_pair("04-stage-2-to-4-duration", plot_duration, 11, 5.5)
save_plot_pair("05-sem-effect-decomposition", plot_sem_effects, 13, 6.5)
save_plot_pair("06-sem-constituent-paths", plot_sem_paths, 15, 7.5)
save_plot_pair("06b-sem-effect-heatmap", plot_sem_heatmap, 13, 7.5)
save_plot_pair("07-primary-lmm-diagnostics", plot_diagnostics, 12, 9)
save_plot_pair("08-transition-completeness-sensitivity", plot_missingness_sensitivity, 14, 12)
save_plot_pair("phenology-overall-timing-summary", summary_plot, 16, 27)

readr::write_csv(primary_effects, file.path(output_dir, "primary-common-shift-effects.csv"))
readr::write_csv(complete_case_effects, file.path(output_dir, "all-three-transitions-sensitivity-effects.csv"))
readr::write_csv(heteroscedastic_effects, file.path(output_dir, "heteroscedastic-lmm-sensitivity-effects.csv"))
readr::write_csv(heteroscedastic_residual_multipliers, file.path(output_dir, "heteroscedastic-residual-sd-multipliers.csv"))
readr::write_csv(stage_effects, file.path(output_dir, "stage-specific-sensitivity-effects.csv"))
readr::write_csv(stage_interaction_tests, file.path(output_dir, "stage-interaction-tests.csv"))
readr::write_csv(duration_effects, file.path(output_dir, "stage-2-to-4-duration-effects.csv"))
readr::write_csv(dplyr::bind_rows(timing_index_data), file.path(output_dir, "stage-centred-timing-index.csv"))
readr::write_csv(spring_swc, file.path(output_dir, "spring-swc-container-summary.csv"))
readr::write_csv(sem_effects, file.path(output_dir, "sem-effect-decomposition.csv"))
readr::write_csv(sem_paths, file.path(output_dir, "sem-constituent-paths.csv"))
readr::write_csv(sem_heatmap_data, file.path(output_dir, "sem-effect-heatmap-data.csv"))
readr::write_csv(sem_bootstrap, file.path(output_dir, "sem-container-cluster-bootstrap.csv"))
readr::write_csv(model_summary, file.path(output_dir, "model-fit-summary.csv"))
readr::write_csv(transition_missingness, file.path(output_dir, "transition-missingness-by-treatment-cell.csv"))
readr::write_csv(missingness_effects, file.path(output_dir, "transition-observation-logistic-effects.csv"))
readr::write_csv(diagnostic_data, file.path(output_dir, "primary-lmm-diagnostic-data.csv"))

saveRDS(
  list(
    settings = list(
      analysis_date = analysis_date,
      stages = stages_keep,
      spring_swc_start = spring_swc_start,
      spring_swc_end = spring_swc_end,
      bootstrap_replicates = n_boot,
      soil_type = soil_type,
      include_soil_treatment = include_soil_treatment
    ),
    transition_data = transition_data,
    transition_models = transition_models,
    complete_case_data = complete_case_data,
    complete_case_models = complete_case_models,
    complete_case_effects = complete_case_effects,
    heteroscedastic_models = heteroscedastic_models,
    heteroscedastic_effects = heteroscedastic_effects,
    heteroscedastic_residual_multipliers = heteroscedastic_residual_multipliers,
    primary_effects = primary_effects,
    stage_effects = stage_effects,
    stage_interaction_tests = stage_interaction_tests,
    duration_data = duration_data,
    duration_models = duration_models,
    duration_effects = duration_effects,
    timing_index_data = timing_index_data,
    sem_data = sem_data,
    sem_fits = sem_fits,
    sem_effects = sem_effects,
    sem_paths = sem_paths,
    sem_heatmap_data = sem_heatmap_data,
    transition_missingness = transition_missingness,
    missingness_models = missingness_models,
    missingness_effects = missingness_effects,
    model_summary = model_summary
  ),
  file.path(output_dir, "phenology-overall-timing-analysis.rds")
)

cat("\nPrimary common-shift effects (raw days; negative = earlier):\n")
print(
  primary_effects %>%
    dplyr::select("species", "treatment_label", "estimate_days_raw", "lower_95_raw", "upper_95_raw", "p_value") %>%
    tibble::as_tibble(),
  n = Inf
)
cat("\nStage-interaction tests:\n")
print(
  stage_interaction_tests %>%
    dplyr::select("species", "test", "chi_square", "df", "p_value") %>%
    tibble::as_tibble(),
  n = Inf
)
cat("\nStage 2–4 duration effects (raw days; positive = longer):\n")
print(
  duration_effects %>%
    dplyr::select("species", "treatment_label", "estimate_days_raw", "lower_95_raw", "upper_95_raw", "p_value") %>%
    tibble::as_tibble(),
  n = Inf
)
cat("\nSEM decomposition (oriented standardized effects):\n")
print(
  sem_effects %>%
    dplyr::filter(.data$component %in% c("Direct", "Indirect via pre-leaf-out SWC", "Path-summed total")) %>%
    dplyr::select("species", "treatment_label", "component", "estimate", "lower_95", "upper_95", "p_boot", "n_boot_success") %>%
    tibble::as_tibble(),
  n = Inf
)
cat("\nModel status:\n")
print(model_summary, n = Inf)
message("Saved the complete phenology collection in: ", normalizePath(output_dir, winslash = "/", mustWork = TRUE))
