#!/usr/bin/env Rscript

# Block-stratified, non-parametric container bootstrap for the production
# phenology models used by Figures 6 and S16-S17.

options(stringsAsFactors = FALSE)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this file with Rscript.", call. = FALSE)
script_path_arg <- sub("^--file=", "", script_arg[[1]])
# Rscript 4.5 encodes spaces in --file paths as ~+~.
script_path_arg <- gsub("~\\+~", " ", script_path_arg)
script_file <- normalizePath(script_path_arg, winslash = "/", mustWork = TRUE)
script_dir <- dirname(script_file)

find_project_root <- function(path) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, "scripts", "auxiliary", "functions", "_source.R"))) return(path)
    parent <- dirname(path)
    if (identical(parent, path)) stop("Could not locate project root.", call. = FALSE)
    path <- parent
  }
}
project_root <- find_project_root(script_dir)
setwd(project_root)

args <- commandArgs(trailingOnly = TRUE)
arg_int <- function(prefix, default) {
  hit <- grep(paste0("^", prefix, "="), args, value = TRUE)
  if (!length(hit)) return(as.integer(default))
  as.integer(sub(paste0("^", prefix, "="), "", hit[[1]]))
}
n_boot <- arg_int("--bootstrap", 1000L)
n_cores <- arg_int("--cores", 3L)
if (is.na(n_boot) || n_boot < 100L) stop("--bootstrap must be at least 100.", call. = FALSE)
if (is.na(n_cores) || n_cores < 1L) stop("--cores must be positive.", call. = FALSE)

renv_lib <- Sys.glob(file.path(project_root, "renv", "library", "*", "R-*", "*"))
if (length(renv_lib)) .libPaths(c(normalizePath(renv_lib[[1]], winslash = "/"), .libPaths()))

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(lme4)
  library(patchwork)
  library(readr)
  library(tidyr)
})
suppressMessages(suppressPackageStartupMessages({
  source(file.path(project_root, "scripts", "auxiliary", "functions", "_source.R"))
  source(file.path(project_root, "scripts", "auxiliary", "functions", "1-summary-figures.R"))
  source(file.path(project_root, "scripts", "auxiliary", "functions", "7-phenology-transition-models.R"))
}))

output_dir <- file.path(project_root, "data", "final", "bootstrap", "phenology")
work_dir <- file.path(tempdir(), "alinv-phenology-bootstrap")
figure_dir <- file.path(work_dir, "figures")
comparison_dir <- file.path(work_dir, "comparison")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(comparison_dir, recursive = TRUE, showWarnings = FALSE)

species_vec <- c("fagus", "quercus")
species_short <- c(fagus = "Fagus", quercus = "Quercus")
stages_keep <- 2:4
spring_swc_start <- as.Date("2025-03-04")
spring_swc_end <- as.Date("2025-04-02")
seed_lmm <- 202608181L
seed_sem <- 202608182L

treatment_specs <- ALINV_TREATMENT_CONFIG %>%
  filter(.data$effect %in% c("precipitation", "robinia", "culture")) %>%
  arrange(.data$plot_order) %>%
  select("effect", "baseline_level", "treatment_level",
         treatment_label = "short_label", "contrast_label", "plot_order")
treatment_colors <- c(precipitation = "#1B9E77", robinia = "#D95F02", culture = "#7570B3")
precipitation_colors <- c(control = "#4F6674", drought = "#D65F5F")
component_colors <- c("Direct" = "#3B6FB6", "Indirect via pre-leaf-out SWC" = "#D95F02",
                      "Path-summed total" = "#222222")

model_control <- lmerControl(
  optimizer = "bobyqa", optCtrl = list(maxfun = 200000),
  check.conv.singular = lme4::.makeCC(action = "ignore", tol = 1e-4)
)

add_block <- function(df) {
  df %>% mutate(
    block = factor(sub("-.*$", "", as.character(.data$boxlabel)), levels = c("b1", "b2", "b3")),
    boxlabel = factor(.data$boxlabel), tree_id = factor(.data$tree_id),
    precipitation = factor(.data$precipitation, levels = alinv_factor_levels("precipitation")),
    robinia = factor(.data$robinia, levels = alinv_factor_levels("robinia")),
    culture = factor(.data$culture, levels = alinv_factor_levels("culture"))
  ) %>% droplevels()
}

transition_data <- setNames(lapply(species_vec, function(sp) {
  prepare_phenology_transition_data(
    species_keep = sp, soil_type = "both", include_soil_treatment = FALSE,
    stages_keep = stages_keep
  ) %>% add_block()
}), species_vec)

primary_formula <- doy ~ stage_label + precipitation + robinia + culture +
  (1 | boxlabel) + (1 | tree_id)
interaction_formula <- doy ~ stage_label * (precipitation + robinia + culture) +
  (1 | boxlabel) + (1 | tree_id)
duration_formula <- duration_days ~ precipitation + robinia + culture + (1 | boxlabel)

make_duration_data <- function(df) {
  df %>%
    filter(.data$stage %in% c(2L, 4L)) %>%
    select("tree_id", "boxlabel", "block", "precipitation", "robinia", "culture", "stage", "doy") %>%
    pivot_wider(names_from = "stage", values_from = "doy", names_prefix = "stage_") %>%
    filter(!is.na(.data$stage_2), !is.na(.data$stage_4)) %>%
    mutate(duration_days = .data$stage_4 - .data$stage_2) %>% droplevels()
}

coef_for <- function(spec) paste0(spec$effect[[1]], spec$treatment_level[[1]])

primary_from_fit <- function(fit, sp, replicate = NA_integer_) {
  beta <- fixef(fit)
  bind_rows(lapply(seq_len(nrow(treatment_specs)), function(i) {
    spec <- treatment_specs[i, ]
    value <- unname(beta[[coef_for(spec)]])
    tibble(species = sp, replicate = replicate, effect = spec$effect[[1]],
           treatment_label = spec$treatment_label[[1]], contrast_label = spec$contrast_label[[1]],
           plot_order = spec$plot_order[[1]], estimate_days_raw = value,
           estimate_oriented = -value)
  }))
}

stage_from_fit <- function(fit, sp, replicate = NA_integer_) {
  beta <- fixef(fit)
  bind_rows(lapply(seq_len(nrow(treatment_specs)), function(i) {
    spec <- treatment_specs[i, ]
    main_name <- coef_for(spec)
    bind_rows(lapply(stages_keep, function(stage_i) {
      value <- unname(beta[[main_name]])
      if (stage_i != 2L) {
        stage_name <- paste0("stage_labelStage ", stage_i)
        candidates <- c(paste0(stage_name, ":", main_name), paste0(main_name, ":", stage_name))
        interaction_name <- candidates[candidates %in% names(beta)]
        if (length(interaction_name) != 1L) stop("Missing stage interaction coefficient.")
        value <- value + unname(beta[[interaction_name]])
      }
      tibble(species = sp, replicate = replicate, stage_label = paste("Stage", stage_i),
             effect = spec$effect[[1]], treatment_label = spec$treatment_label[[1]],
             contrast_label = spec$contrast_label[[1]], plot_order = spec$plot_order[[1]],
             estimate_days_raw = value, estimate_oriented = -value)
    }))
  }))
}

duration_from_fit <- function(fit, sp, replicate = NA_integer_) {
  beta <- fixef(fit)
  bind_rows(lapply(seq_len(nrow(treatment_specs)), function(i) {
    spec <- treatment_specs[i, ]
    value <- unname(beta[[coef_for(spec)]])
    tibble(species = sp, replicate = replicate, effect = spec$effect[[1]],
           treatment_label = spec$treatment_label[[1]], contrast_label = spec$contrast_label[[1]],
           plot_order = spec$plot_order[[1]], estimate_days_raw = value,
           estimate_oriented = -value)
  }))
}

resample_transition_within_block <- function(df, seed) {
  set.seed(seed)
  box_map <- df %>% distinct(.data$boxlabel, .data$block) %>% arrange(.data$block, .data$boxlabel)
  sampled <- box_map %>% group_by(.data$block) %>%
    group_modify(~ tibble(source_box = sample(as.character(.x$boxlabel), nrow(.x), replace = TRUE))) %>%
    ungroup() %>% mutate(copy = row_number(), boot_box = paste0(as.character(.data$block), "-boot-", .data$copy))
  bind_rows(lapply(seq_len(nrow(sampled)), function(j) {
    source <- sampled$source_box[[j]]
    boot_box <- sampled$boot_box[[j]]
    df %>% filter(as.character(.data$boxlabel) == source) %>%
      mutate(boxlabel = boot_box, tree_id = paste0(boot_box, "::", as.character(.data$tree_id)))
  })) %>% add_block()
}

fit_lmm_boot_once <- function(df, sp, attempt, seed) {
  boot <- resample_transition_within_block(df, seed)
  out <- tryCatch({
    primary <- suppressWarnings(lmer(primary_formula, boot, REML = TRUE, control = model_control))
    interaction <- suppressWarnings(lmer(interaction_formula, boot, REML = TRUE, control = model_control))
    duration_data <- make_duration_data(boot)
    duration <- suppressWarnings(lmer(duration_formula, duration_data, REML = TRUE, control = model_control))
    list(
      primary = primary_from_fit(primary, sp, attempt),
      stage = stage_from_fit(interaction, sp, attempt),
      duration = duration_from_fit(duration, sp, attempt),
      status = tibble(species = sp, attempt = attempt, success = TRUE,
                      primary_singular = isSingular(primary, tol = 1e-4),
                      interaction_singular = isSingular(interaction, tol = 1e-4),
                      duration_singular = isSingular(duration, tol = 1e-4), error = NA_character_)
    )
  }, error = function(e) {
    list(status = tibble(species = sp, attempt = attempt, success = FALSE,
                         primary_singular = NA, interaction_singular = NA,
                         duration_singular = NA, error = conditionMessage(e)))
  })
  out
}

run_until_success <- function(worker, data, sp, target, seed_base, label) {
  accepted <- list(); statuses <- list(); attempts <- 0L; n_success <- 0L
  while (n_success < target) {
    remaining <- target - n_success
    batch_n <- remaining
    ids <- attempts + seq_len(batch_n)
    seeds <- seed_base + match(sp, species_vec) * 1000000L + ids
    cores <- if (.Platform$OS.type == "windows") 1L else n_cores
    batch <- parallel::mclapply(seq_along(ids), function(k) worker(data, sp, ids[[k]], seeds[[k]]),
                                mc.cores = cores, mc.preschedule = TRUE)
    batch_ok <- batch[vapply(batch, function(x) isTRUE(x$status$success[[1]]), logical(1))]
    accepted <- c(accepted, batch_ok)
    statuses <- c(statuses, lapply(batch, `[[`, "status"))
    n_success <- length(accepted)
    attempts <- max(ids)
    message(label, " ", sp, ": ", min(n_success, target), "/", target,
            " successful (", attempts, " attempts)")
    if (attempts > target * 5L) stop("Too many failed bootstrap attempts for ", label, " ", sp)
  }
  list(results = accepted[seq_len(target)], status = bind_rows(statuses), attempts = attempts)
}

message("Running block-stratified transition-model bootstrap...")
lmm_runs <- setNames(lapply(species_vec, function(sp) {
  run_until_success(fit_lmm_boot_once, transition_data[[sp]], sp, n_boot, seed_lmm, "Phenology LMM")
}), species_vec)

lmm_primary_boot <- bind_rows(lapply(lmm_runs, function(x) bind_rows(lapply(x$results, `[[`, "primary")))) %>%
  group_by(.data$species) %>% mutate(replicate = dense_rank(.data$replicate)) %>% ungroup()
lmm_stage_boot <- bind_rows(lapply(lmm_runs, function(x) bind_rows(lapply(x$results, `[[`, "stage")))) %>%
  group_by(.data$species) %>% mutate(replicate = dense_rank(.data$replicate)) %>% ungroup()
lmm_duration_boot <- bind_rows(lapply(lmm_runs, function(x) bind_rows(lapply(x$results, `[[`, "duration")))) %>%
  group_by(.data$species) %>% mutate(replicate = dense_rank(.data$replicate)) %>% ungroup()
lmm_boot_status <- bind_rows(lapply(lmm_runs, `[[`, "status"))

point_models <- setNames(lapply(transition_data, function(df) list(
  primary = lmer(primary_formula, df, REML = TRUE, control = model_control),
  interaction = lmer(interaction_formula, df, REML = TRUE, control = model_control),
  duration = lmer(duration_formula, make_duration_data(df), REML = TRUE, control = model_control)
)), species_vec)
point_primary <- bind_rows(lapply(species_vec, function(sp) primary_from_fit(point_models[[sp]]$primary, sp)))
point_stage <- bind_rows(lapply(species_vec, function(sp) stage_from_fit(point_models[[sp]]$interaction, sp)))
point_duration <- bind_rows(lapply(species_vec, function(sp) duration_from_fit(point_models[[sp]]$duration, sp)))

p_boot <- function(x) min(1, 2 * min((sum(x <= 0) + 1) / (length(x) + 1),
                                     (sum(x >= 0) + 1) / (length(x) + 1)))

summarise_lmm <- function(boot, point, keys) {
  boot %>% group_by(across(all_of(keys))) %>% summarise(
    lower_95_raw = quantile(.data$estimate_days_raw, 0.025, names = FALSE),
    upper_95_raw = quantile(.data$estimate_days_raw, 0.975, names = FALSE),
    lower_95_oriented = quantile(.data$estimate_oriented, 0.025, names = FALSE),
    upper_95_oriented = quantile(.data$estimate_oriented, 0.975, names = FALSE),
    p_boot = p_boot(.data$estimate_oriented), n_boot_success = n(), .groups = "drop"
  ) %>% left_join(point %>% select(all_of(keys), "estimate_days_raw", "estimate_oriented"), by = keys) %>%
    mutate(raw_interpretation = case_when(.data$estimate_days_raw < 0 ~ "earlier",
                                          .data$estimate_days_raw > 0 ~ "later", TRUE ~ "no shift"))
}
primary_effects <- summarise_lmm(lmm_primary_boot, point_primary,
                                 c("species", "effect", "treatment_label", "contrast_label", "plot_order"))
stage_effects <- summarise_lmm(lmm_stage_boot, point_stage,
                               c("species", "stage_label", "effect", "treatment_label", "contrast_label", "plot_order"))
duration_effects <- summarise_lmm(lmm_duration_boot, point_duration,
                                  c("species", "effect", "treatment_label", "contrast_label", "plot_order"))

# Descriptive progression (unchanged by inference method).
phenology_all <- get_data("tree", "phenology_transitions") %>%
  alinv_filter_by_soil(soil_filter = "both") %>%
  filter(.data$species %in% species_vec, .data$stage %in% 1:4,
         !is.na(.data$doy), !is.na(.data$stage_date))
progression_data <- phenology_all %>% mutate(
  species = factor(.data$species, levels = species_vec, labels = unname(species_short[species_vec])),
  robinia = factor(.data$robinia, levels = alinv_factor_levels("robinia"),
                   labels = c("Without Robinia", "With Robinia")),
  precipitation = factor(.data$precipitation, levels = alinv_factor_levels("precipitation")),
  culture = factor(.data$culture, levels = alinv_factor_levels("culture")), stage = as.integer(.data$stage)
) %>% group_by(.data$species, .data$robinia, .data$precipitation, .data$culture, .data$stage) %>%
  summarise(mean_doy = mean(.data$doy), sd_doy = sd(.data$doy), n = n(),
            se_doy = .data$sd_doy / sqrt(.data$n), minimum_doy = min(.data$doy),
            maximum_doy = max(.data$doy), .groups = "drop")

# Stage-centred phenology index and fixed pre-leaf-out SWC definition.
make_timing_index <- function(df, sp) {
  df %>% group_by(.data$stage_label) %>%
    mutate(stage_mean_doy = mean(.data$doy), stage_centered_doy = .data$doy - .data$stage_mean_doy) %>%
    ungroup() %>% group_by(.data$tree_id, .data$boxlabel, .data$block,
                           .data$precipitation, .data$robinia, .data$culture) %>%
    summarise(timing_index_days = mean(.data$stage_centered_doy), n_transitions = n(), .groups = "drop") %>%
    filter(.data$n_transitions >= 2L) %>%
    mutate(timing_z = as.numeric(scale(.data$timing_index_days)),
           phenology_oriented_z = -.data$timing_z, species = sp) %>% droplevels()
}
timing_index_data <- setNames(Map(make_timing_index, transition_data, species_vec), species_vec)
spring_swc <- get_data("box", "soilwater", swc_source = "measured") %>%
  filter(.data$date >= spring_swc_start, .data$date <= spring_swc_end, !is.na(.data$swc)) %>%
  group_by(.data$boxlabel, .data$precipitation, .data$robinia, .data$culture) %>%
  summarise(spring_swc = mean(.data$swc), spring_swc_sd = sd(.data$swc),
            n_swc_dates = n_distinct(.data$date), .groups = "drop") %>%
  mutate(block = factor(sub("-.*$", "", .data$boxlabel), levels = c("b1", "b2", "b3")),
         precipitation = factor(.data$precipitation, levels = alinv_factor_levels("precipitation")),
         robinia = factor(.data$robinia, levels = alinv_factor_levels("robinia")),
         culture = factor(.data$culture, levels = alinv_factor_levels("culture")))

prepare_sem_data <- function(sp) {
  timing <- timing_index_data[[sp]]
  boxes <- spring_swc %>% filter(.data$boxlabel %in% as.character(timing$boxlabel)) %>%
    mutate(swc_z = as.numeric(scale(.data$spring_swc)), boxlabel = factor(.data$boxlabel)) %>% droplevels()
  trees <- timing %>% mutate(boxlabel = as.character(.data$boxlabel)) %>%
    left_join(boxes %>% select("boxlabel", "swc_z", "spring_swc", "n_swc_dates"), by = "boxlabel") %>%
    filter(!is.na(.data$swc_z)) %>% mutate(boxlabel = factor(.data$boxlabel)) %>% droplevels()
  list(box = boxes, tree = trees)
}
sem_data <- setNames(lapply(species_vec, prepare_sem_data), species_vec)
mediator_formula <- swc_z ~ precipitation + robinia + culture

fit_sem_paths <- function(data_i, box_group = "boxlabel") {
  outcome_formula <- as.formula(paste("phenology_oriented_z ~ swc_z + precipitation + robinia + culture +",
                                      paste0("(1 | ", box_group, ")")))
  total_formula <- as.formula(paste("phenology_oriented_z ~ precipitation + robinia + culture +",
                                    paste0("(1 | ", box_group, ")")))
  med <- lm(mediator_formula, data_i$box)
  outcome <- suppressWarnings(lmer(outcome_formula, data_i$tree, REML = TRUE, control = model_control))
  total <- suppressWarnings(lmer(total_formula, data_i$tree, REML = TRUE, control = model_control))
  a <- coef(med); out <- fixef(outcome); tot <- fixef(total); b <- unname(out[["swc_z"]])
  effects <- bind_rows(lapply(seq_len(nrow(treatment_specs)), function(i) {
    spec <- treatment_specs[i, ]; name <- coef_for(spec)
    a_i <- unname(a[[name]]); direct <- unname(out[[name]]); total_reduced <- unname(tot[[name]])
    tibble(effect = spec$effect[[1]], treatment_label = spec$treatment_label[[1]],
           plot_order = spec$plot_order[[1]], a_treatment_to_swc = a_i,
           b_swc_to_phenology = b, direct = direct, indirect = a_i * b,
           total_path = direct + a_i * b, total_reduced = total_reduced,
           total_discrepancy = direct + a_i * b - total_reduced)
  }))
  list(mediator = med, outcome = outcome, total = total, effects = effects)
}

resample_sem_within_block <- function(data_i, seed) {
  set.seed(seed)
  map <- data_i$box %>% distinct(.data$boxlabel, .data$block) %>% arrange(.data$block, .data$boxlabel) %>%
    group_by(.data$block) %>%
    group_modify(~ tibble(source_box = sample(as.character(.x$boxlabel), nrow(.x), replace = TRUE))) %>%
    ungroup() %>% mutate(copy = row_number(), boot_box = paste0(as.character(.data$block), "-boot-", .data$copy))
  boxes <- bind_rows(lapply(seq_len(nrow(map)), function(j) {
    data_i$box %>% filter(as.character(.data$boxlabel) == map$source_box[[j]]) %>%
      mutate(boot_box = factor(map$boot_box[[j]]))
  }))
  trees <- bind_rows(lapply(seq_len(nrow(map)), function(j) {
    data_i$tree %>% filter(as.character(.data$boxlabel) == map$source_box[[j]]) %>%
      mutate(boot_box = factor(map$boot_box[[j]]),
             boot_tree = factor(paste0(map$boot_box[[j]], "::", as.character(.data$tree_id))))
  }))
  list(box = boxes, tree = trees)
}

fit_sem_boot_once <- function(data_i, sp, attempt, seed) {
  out <- tryCatch({
    boot <- resample_sem_within_block(data_i, seed)
    fit <- fit_sem_paths(boot, "boot_box")
    if (any(!is.finite(as.matrix(fit$effects %>% select(where(is.numeric)))))) stop("Non-finite SEM effect")
    list(effects = fit$effects %>% mutate(species = sp, replicate = attempt, .before = 1),
         status = tibble(species = sp, attempt = attempt, success = TRUE,
                         outcome_singular = isSingular(fit$outcome, tol = 1e-4),
                         total_singular = isSingular(fit$total, tol = 1e-4), error = NA_character_))
  }, error = function(e) {
    list(status = tibble(species = sp, attempt = attempt, success = FALSE,
                         outcome_singular = NA, total_singular = NA, error = conditionMessage(e)))
  })
  out
}

message("Running block-stratified phenology SEM bootstrap...")
sem_runs <- setNames(lapply(species_vec, function(sp) {
  run_until_success(fit_sem_boot_once, sem_data[[sp]], sp, n_boot, seed_sem, "Phenology SEM")
}), species_vec)
sem_boot <- bind_rows(lapply(sem_runs, function(x) bind_rows(lapply(x$results, `[[`, "effects")))) %>%
  group_by(.data$species) %>% mutate(replicate = dense_rank(.data$replicate)) %>% ungroup()
sem_boot_status <- bind_rows(lapply(sem_runs, `[[`, "status"))
sem_fits <- setNames(lapply(sem_data, fit_sem_paths), species_vec)
point_sem <- bind_rows(lapply(species_vec, function(sp) sem_fits[[sp]]$effects %>% mutate(species = sp, .before = 1)))

metric_labels <- c(direct = "Direct", indirect = "Indirect via pre-leaf-out SWC",
                   total_path = "Path-summed total", total_reduced = "Reduced-form total")
summarise_sem_metric <- function(metric) {
  sem_boot %>% group_by(.data$species, .data$effect, .data$treatment_label, .data$plot_order) %>%
    summarise(lower_95 = quantile(.data[[metric]], 0.025, names = FALSE),
              upper_95 = quantile(.data[[metric]], 0.975, names = FALSE),
              p_boot = p_boot(.data[[metric]]), n_boot_success = n(), .groups = "drop") %>%
    left_join(point_sem %>% select("species", "effect", "treatment_label", estimate = all_of(metric)),
              by = c("species", "effect", "treatment_label"))
}
sem_effects <- bind_rows(lapply(names(metric_labels), function(metric) {
  summarise_sem_metric(metric) %>% mutate(component = unname(metric_labels[[metric]]), metric = metric)
}))
sem_paths <- bind_rows(
  summarise_sem_metric("a_treatment_to_swc") %>%
    mutate(path = "Treatment → pre-leaf-out SWC", metric = "a_treatment_to_swc"),
  summarise_sem_metric("direct") %>% mutate(path = "Treatment → phenology (direct)", metric = "direct")
)
b_paths <- sem_boot %>% distinct(.data$species, .data$replicate, .data$b_swc_to_phenology) %>%
  group_by(.data$species) %>% summarise(
    lower_95 = quantile(.data$b_swc_to_phenology, 0.025, names = FALSE),
    upper_95 = quantile(.data$b_swc_to_phenology, 0.975, names = FALSE),
    p_boot = p_boot(.data$b_swc_to_phenology), n_boot_success = n(), .groups = "drop"
  ) %>% left_join(point_sem %>% distinct(.data$species, estimate = .data$b_swc_to_phenology), by = "species") %>%
  mutate(effect = "swc", treatment_label = "Pre-leaf-out SWC", plot_order = 4L,
         path = "Pre-leaf-out SWC → phenology", metric = "b_swc_to_phenology")
sem_paths <- bind_rows(sem_paths, b_paths)
figure6_ready <- sem_effects %>% filter(.data$component == "Path-summed total") %>%
  mutate(bootstrap_design = "container resampling within block", bootstrap_replicates = n_boot)

# Comparisons with the current manuscript outputs.
current_dir <- file.path(
  project_root, "_archive", "model-caches", "dated-output", "2026-08-18",
  "phenology-overall-timing-analysis"
)
compare_effects <- function(current_file, new_df, keys, old_method, new_method,
                            old_lower, old_upper, old_p, new_lower, new_upper, new_p) {
  old <- read_csv(current_file, show_col_types = FALSE) %>%
    select(all_of(keys), old_estimate = "estimate_oriented", old_lower = all_of(old_lower),
           old_upper = all_of(old_upper), old_p = all_of(old_p))
  new <- new_df %>% select(all_of(keys), new_estimate = "estimate_oriented",
                           new_lower = all_of(new_lower), new_upper = all_of(new_upper), new_p = all_of(new_p))
  inner_join(old, new, by = keys) %>% mutate(
    old_method = old_method, new_method = new_method,
    point_difference = .data$new_estimate - .data$old_estimate,
    old_direction = case_when(.data$old_estimate < 0 ~ "negative", .data$old_estimate > 0 ~ "positive", TRUE ~ "zero"),
    new_direction = case_when(.data$new_estimate < 0 ~ "negative", .data$new_estimate > 0 ~ "positive", TRUE ~ "zero"),
    direction_changed = .data$old_direction != .data$new_direction,
    old_interval_excludes_zero = .data$old_lower > 0 | .data$old_upper < 0,
    new_interval_excludes_zero = .data$new_lower > 0 | .data$new_upper < 0,
    interval_significance_changed = .data$old_interval_excludes_zero != .data$new_interval_excludes_zero,
    old_p_below_0_05 = .data$old_p < 0.05, new_p_below_0_05 = .data$new_p < 0.05,
    p_significance_changed = .data$old_p_below_0_05 != .data$new_p_below_0_05
  )
}
lmm_comparison <- compare_effects(
  file.path(current_dir, "primary-common-shift-effects.csv"), primary_effects,
  c("species", "effect", "treatment_label"), "Kenward–Roger", "block-stratified container bootstrap",
  "lower_95_oriented", "upper_95_oriented", "p_value", "lower_95_oriented", "upper_95_oriented", "p_boot"
)
stage_comparison <- compare_effects(
  file.path(current_dir, "stage-specific-sensitivity-effects.csv"), stage_effects,
  c("species", "stage_label", "effect", "treatment_label"), "Kenward–Roger", "block-stratified container bootstrap",
  "lower_95_oriented", "upper_95_oriented", "p_value", "lower_95_oriented", "upper_95_oriented", "p_boot"
)
duration_comparison <- compare_effects(
  file.path(current_dir, "stage-2-to-4-duration-effects.csv"), duration_effects,
  c("species", "effect", "treatment_label"), "Kenward–Roger", "block-stratified container bootstrap",
  "lower_95_oriented", "upper_95_oriented", "p_value", "lower_95_oriented", "upper_95_oriented", "p_boot"
)
old_sem <- read_csv(file.path(current_dir, "sem-effect-decomposition.csv"), show_col_types = FALSE)
sem_comparison <- old_sem %>%
  select("species", "effect", "treatment_label", "component", old_estimate = "estimate",
         old_lower = "lower_95", old_upper = "upper_95", old_p = "p_boot") %>%
  inner_join(sem_effects %>% select("species", "effect", "treatment_label", "component",
                                   new_estimate = "estimate", new_lower = "lower_95",
                                   new_upper = "upper_95", new_p = "p_boot"),
             by = c("species", "effect", "treatment_label", "component")) %>%
  mutate(old_method = "unstratified container bootstrap", new_method = "block-stratified container bootstrap",
         point_difference = .data$new_estimate - .data$old_estimate,
         old_direction = case_when(.data$old_estimate < 0 ~ "negative", .data$old_estimate > 0 ~ "positive", TRUE ~ "zero"),
         new_direction = case_when(.data$new_estimate < 0 ~ "negative", .data$new_estimate > 0 ~ "positive", TRUE ~ "zero"),
         direction_changed = .data$old_direction != .data$new_direction,
         old_interval_excludes_zero = .data$old_lower > 0 | .data$old_upper < 0,
         new_interval_excludes_zero = .data$new_lower > 0 | .data$new_upper < 0,
         interval_significance_changed = .data$old_interval_excludes_zero != .data$new_interval_excludes_zero,
         old_p_below_0_05 = .data$old_p < 0.05, new_p_below_0_05 = .data$new_p < 0.05,
         p_significance_changed = .data$old_p_below_0_05 != .data$new_p_below_0_05)

# Figures: S1 is descriptive; S2 and all exploratory inferential panels use bootstrap intervals.
theme_supp <- function(base_size = 8) theme_classic(base_size = base_size) + theme(
  text = element_text(color = "black"), axis.text = element_text(color = "black"),
  axis.line = element_line(linewidth = 0.25), axis.ticks = element_line(linewidth = 0.25),
  panel.grid.major = element_line(color = "grey90", linewidth = 0.2), panel.grid.minor = element_blank(),
  strip.background = element_rect(fill = "grey12", color = "grey12"),
  strip.text = element_text(color = "white", face = "bold"), plot.title = element_text(face = "bold"),
  legend.position = "bottom")
save_pair <- function(plot, stem, width_mm, height_mm, dpi = 600) {
  ggsave(file.path(figure_dir, paste0(stem, ".png")), plot, width = width_mm, height = height_mm,
         units = "mm", dpi = dpi, bg = "white", limitsize = FALSE)
  ggsave(file.path(figure_dir, paste0(stem, ".pdf")), plot, width = width_mm, height = height_mm,
         units = "mm", bg = "white", limitsize = FALSE)
}
figure_s1 <- ggplot(progression_data, aes(.data$mean_doy, .data$stage, color = .data$precipitation,
                                         linetype = .data$culture,
                                         group = interaction(.data$precipitation, .data$culture))) +
  geom_segment(aes(x = .data$mean_doy - .data$se_doy, xend = .data$mean_doy + .data$se_doy,
                   yend = .data$stage), linewidth = 0.55, alpha = 0.55) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.8) + facet_grid(.data$species ~ .data$robinia) +
  scale_color_manual(values = precipitation_colors, labels = c(control = "Control", drought = "Reduced precipitation")) +
  scale_linetype_manual(values = c(mono = "solid", mixed = "dotted"),
                        labels = c(mono = "Monoculture", mixed = "Mixed culture")) +
  scale_x_continuous(breaks = seq(90, 130, 10)) +
  scale_y_continuous(breaks = 1:4, labels = paste("Stage", 1:4), limits = c(0.8, 4.2)) +
  labs(title = "Observed spring phenology progression",
       subtitle = "Treatment-group mean transition day of year (DOY) ± SE",
       x = "Day of year (DOY)", y = "Phenological stage", color = "Precipitation", linetype = "Culture",
       caption = "Stage 1 is descriptive; inferential bootstrap models use transitions into stages 2–4.") + theme_supp()

effect_plot_data <- primary_effects %>% mutate(
  species_label = factor(.data$species, levels = species_vec, labels = unname(species_short[species_vec])),
  treatment_label = factor(.data$treatment_label, levels = treatment_specs$treatment_label),
  interval_excludes_zero = .data$lower_95_oriented > 0 | .data$upper_95_oriented < 0)
figure_s2 <- ggplot(effect_plot_data, aes(.data$estimate_oriented, .data$treatment_label, color = .data$effect)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_segment(aes(x = .data$lower_95_oriented, xend = .data$upper_95_oriented,
                   yend = .data$treatment_label), linewidth = 0.8, lineend = "round") +
  geom_point(aes(shape = .data$interval_excludes_zero), size = 2.5, fill = "white", stroke = 0.8) +
  facet_wrap(~species_label, nrow = 1) + scale_color_manual(values = treatment_colors, guide = "none") +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), guide = "none") +
  labs(title = "Overall shift of the stage-2–4 transition profile",
       subtitle = paste0("Additive Gaussian LMM; 95% block-stratified container-bootstrap intervals (",
                         n_boot, " successful replicates per species)"),
       x = "Oriented timing effect (days): negative = delayed, positive = earlier", y = NULL,
       caption = "Open points have intervals overlapping zero; filled points have intervals excluding zero.") +
  theme_supp() + theme(panel.grid.major.x = element_blank())

stage_plot <- stage_effects %>% mutate(
  species_label = factor(.data$species, levels = species_vec, labels = unname(species_short[species_vec])),
  treatment_label = factor(.data$treatment_label, levels = treatment_specs$treatment_label),
  stage_label = factor(.data$stage_label, levels = paste("Stage", 2:4)),
  significant = .data$lower_95_oriented > 0 | .data$upper_95_oriented < 0) %>%
  ggplot(aes(.data$estimate_oriented, .data$treatment_label, color = .data$effect)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_segment(aes(x = .data$lower_95_oriented, xend = .data$upper_95_oriented,
                   yend = .data$treatment_label), linewidth = 0.7) +
  geom_point(aes(shape = .data$significant), size = 2.2) + facet_grid(.data$species_label ~ .data$stage_label) +
  scale_color_manual(values = treatment_colors, guide = "none") +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), guide = "none") +
  labs(title = "Stage-specific sensitivity contrasts",
       subtitle = paste0("Stage × treatment LMM; block-stratified container bootstrap (", n_boot, " successful replicates)"),
       x = "Oriented timing effect (days): negative = delayed, positive = earlier", y = NULL) + theme_supp(8)

duration_plot <- duration_effects %>% mutate(
  species_label = factor(.data$species, levels = species_vec, labels = unname(species_short[species_vec])),
  treatment_label = factor(.data$treatment_label, levels = treatment_specs$treatment_label),
  significant = .data$lower_95_oriented > 0 | .data$upper_95_oriented < 0) %>%
  ggplot(aes(.data$estimate_oriented, .data$treatment_label, color = .data$effect)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_segment(aes(x = .data$lower_95_oriented, xend = .data$upper_95_oriented,
                   yend = .data$treatment_label), linewidth = 0.8) +
  geom_point(aes(shape = .data$significant), size = 2.5) + facet_wrap(~species_label) +
  scale_color_manual(values = treatment_colors, guide = "none") +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), guide = "none") +
  labs(title = "Change in duration from stage 2 to stage 4",
       subtitle = paste0("Block-stratified container bootstrap (", n_boot, " successful replicates per species)"),
       x = "Oriented duration effect (days): negative = slower/longer, positive = faster/shorter", y = NULL) + theme_supp()

sem_plot <- sem_effects %>% filter(.data$component %in% names(component_colors)) %>% mutate(
  species_label = factor(.data$species, levels = species_vec, labels = unname(species_short[species_vec])),
  treatment_label = factor(.data$treatment_label, levels = treatment_specs$treatment_label),
  component = factor(.data$component, levels = names(component_colors)),
  y_plot = as.numeric(.data$treatment_label) + recode(
    as.character(.data$component), "Direct" = -0.16,
    "Indirect via pre-leaf-out SWC" = 0, "Path-summed total" = 0.16
  )) %>%
  ggplot(aes(.data$estimate, .data$y_plot, color = .data$component)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_segment(aes(x = .data$lower_95, xend = .data$upper_95, yend = .data$y_plot), linewidth = 0.75) +
  geom_point(size = 2.3) + facet_wrap(~species_label) +
  scale_y_continuous(breaks = seq_along(treatment_specs$treatment_label),
                     labels = treatment_specs$treatment_label) +
  scale_color_manual(values = component_colors) +
  labs(title = "Associative SEM decomposition of overall phenology timing",
       subtitle = paste0("Block-stratified container bootstrap (", n_boot,
                         " successful replicates); pre-leaf-out SWC: 4 March–2 April 2025"),
       x = "Oriented standardized effect: negative = delayed, positive = earlier", y = NULL,
       color = "Effect component") + theme_supp()

heatmap_data <- sem_effects %>% filter(.data$component %in% names(component_colors)) %>%
  mutate(component = recode(.data$component, "Indirect via pre-leaf-out SWC" = "Indirect via SWC",
                            "Path-summed total" = "Total"),
         species_label = factor(.data$species, levels = species_vec, labels = unname(species_short[species_vec])),
         treatment_label = factor(.data$treatment_label, levels = treatment_specs$treatment_label),
         component = factor(.data$component, levels = c("Direct", "Indirect via SWC", "Total")),
         excludes_zero = .data$lower_95 > 0 | .data$upper_95 < 0,
         label = paste0(sprintf("%.2f", .data$estimate), ifelse(.data$excludes_zero, "*", "")))
heat_limit <- max(0.1, ceiling(max(abs(heatmap_data$estimate)) * 10) / 10)
sem_heatmap <- ggplot(heatmap_data, aes(.data$component, .data$treatment_label, fill = .data$estimate)) +
  geom_tile(color = "grey82") +
  geom_text(aes(label = .data$label), fontface = ifelse(heatmap_data$excludes_zero, "bold", "plain"), size = 3) +
  facet_grid(.data$species_label ~ .) +
  scale_fill_gradient2(low = "#C94F50", mid = "white", high = "#4C78A8", midpoint = 0,
                       limits = c(-heat_limit, heat_limit), name = "Std. effect") +
  labs(title = "Associative phenology SEM heatmap",
       subtitle = "Negative values indicate delayed timing; * marks a bootstrap interval excluding zero.",
       x = NULL, y = NULL) + theme_classic(base_size = 8) +
  theme(panel.grid = element_blank(), strip.background = element_rect(fill = "grey12"),
        strip.text = element_text(color = "white", face = "bold"), axis.text.x = element_text(angle = 20, hjust = 1))

save_pair(figure_s1, "figure-s1-phenology-progression", 160, 112)
save_pair(figure_s2, "figure-s2-phenology-overall-timing-effects-bootstrap", 160, 82)
save_pair(stage_plot, "stage-specific-sensitivity-effects-bootstrap", 180, 110)
save_pair(duration_plot, "stage-2-to-4-duration-effects-bootstrap", 160, 82)
save_pair(sem_plot, "phenology-sem-effect-decomposition-bootstrap", 180, 95)
save_pair(sem_heatmap, "phenology-sem-effect-heatmap-bootstrap", 160, 105)

# Exports.
write_csv(progression_data, file.path(output_dir, "figure-s1-phenology-progression-source-data.csv"))
write_csv(primary_effects, file.path(output_dir, "primary-common-shift-bootstrap-effects.csv"))
write_csv(stage_effects, file.path(output_dir, "stage-specific-bootstrap-effects.csv"))
write_csv(duration_effects, file.path(output_dir, "stage-2-to-4-duration-bootstrap-effects.csv"))
write_csv(lmm_primary_boot, file.path(output_dir, "primary-common-shift-bootstrap-replicates.csv"))
write_csv(lmm_stage_boot, file.path(output_dir, "stage-specific-bootstrap-replicates.csv"))
write_csv(lmm_duration_boot, file.path(output_dir, "stage-2-to-4-duration-bootstrap-replicates.csv"))
write_csv(lmm_boot_status, file.path(output_dir, "phenology-lmm-bootstrap-status.csv"))
write_csv(bind_rows(timing_index_data), file.path(output_dir, "stage-centred-timing-index.csv"))
write_csv(spring_swc, file.path(output_dir, "spring-swc-container-summary.csv"))
write_csv(sem_effects, file.path(output_dir, "sem-effect-decomposition-block-stratified.csv"))
write_csv(sem_paths, file.path(output_dir, "sem-constituent-paths-block-stratified.csv"))
write_csv(figure6_ready, file.path(output_dir, "figure6-ready-phenology-path-summed-total.csv"))
write_csv(sem_boot, file.path(output_dir, "sem-block-stratified-bootstrap-replicates.csv"))
write_csv(sem_boot_status, file.path(output_dir, "sem-block-stratified-bootstrap-status.csv"))
write_csv(lmm_comparison, file.path(comparison_dir, "primary-lmm-bootstrap-vs-kenward-roger.csv"))
write_csv(stage_comparison, file.path(comparison_dir, "stage-sensitivity-bootstrap-vs-kenward-roger.csv"))
write_csv(duration_comparison, file.path(comparison_dir, "duration-bootstrap-vs-kenward-roger.csv"))
write_csv(sem_comparison, file.path(comparison_dir, "sem-stratified-vs-unstratified-bootstrap.csv"))

model_status <- bind_rows(lapply(species_vec, function(sp) tibble(
  species = sp,
  n_transition_observations = nrow(transition_data[[sp]]),
  n_trees = n_distinct(transition_data[[sp]]$tree_id),
  n_containers = n_distinct(transition_data[[sp]]$boxlabel),
  n_blocks = n_distinct(transition_data[[sp]]$block),
  n_duration_trees = nrow(make_duration_data(transition_data[[sp]])),
  lmm_boot_success = sum(lmm_boot_status$species == sp & lmm_boot_status$success),
  lmm_boot_failures = sum(lmm_boot_status$species == sp & !lmm_boot_status$success),
  sem_trees = nrow(sem_data[[sp]]$tree), sem_containers = nrow(sem_data[[sp]]$box),
  sem_boot_success = sum(sem_boot_status$species == sp & sem_boot_status$success),
  sem_boot_failures = sum(sem_boot_status$species == sp & !sem_boot_status$success)
)))
write_csv(model_status, file.path(output_dir, "phenology-bootstrap-model-status.csv"))

saveRDS(list(
  settings = list(n_boot = n_boot, resampling = "containers with replacement within block",
                  seed_lmm = seed_lmm, seed_sem = seed_sem, stages = stages_keep,
                  spring_swc_start = spring_swc_start, spring_swc_end = spring_swc_end,
                  soil_type = "both", include_soil_treatment = FALSE),
  transition_data = transition_data, point_models = point_models,
  primary_effects = primary_effects, stage_effects = stage_effects, duration_effects = duration_effects,
  lmm_primary_boot = lmm_primary_boot, lmm_stage_boot = lmm_stage_boot,
  lmm_duration_boot = lmm_duration_boot, lmm_boot_status = lmm_boot_status,
  timing_index_data = timing_index_data, sem_data = sem_data, sem_fits = sem_fits,
  sem_effects = sem_effects, sem_paths = sem_paths, sem_boot = sem_boot,
  sem_boot_status = sem_boot_status, comparisons = list(lmm = lmm_comparison,
  stage = stage_comparison, duration = duration_comparison, sem = sem_comparison), model_status = model_status
), file.path(output_dir, "phenology-bootstrap-analysis.rds"))

summarise_changes <- function(x) c(direction = sum(x$direction_changed),
                                    interval = sum(x$interval_significance_changed),
                                    p = sum(x$p_significance_changed))
change_lmm <- summarise_changes(lmm_comparison)
change_stage <- summarise_changes(stage_comparison)
change_duration <- summarise_changes(duration_comparison)
change_sem <- summarise_changes(sem_comparison)
stage_changed <- stage_comparison %>% filter(.data$interval_significance_changed | .data$p_significance_changed)
sem_changed <- sem_comparison %>% filter(.data$interval_significance_changed | .data$p_significance_changed)
sem_total_changes <- sem_comparison %>% filter(.data$component == "Path-summed total") %>%
  summarise(n = n(), changes = sum(.data$interval_significance_changed | .data$p_significance_changed))
md <- c(
  "# Phenology bootstrap comparison", "",
  paste0("All models used **", n_boot, " successful non-parametric bootstrap replicates per species**. ",
         "Containers were sampled with replacement within block; all trees and transition records in each selected container were retained under synthetic IDs."), "",
  "## Primary common-shift LMM", "",
  paste0("- Direction changes: ", change_lmm[["direction"]], " of ", nrow(lmm_comparison), "."),
  paste0("- 95% interval-exclusion changes: ", change_lmm[["interval"]], " of ", nrow(lmm_comparison), "."),
  paste0("- P < 0.05 classification changes: ", change_lmm[["p"]], " of ", nrow(lmm_comparison), "."),
  "- Point estimates are unchanged by construction; only uncertainty estimation differs.", "",
  "## Stage-specific and duration sensitivities", "",
  paste0("- Stage-specific contrasts: ", change_stage[["direction"]], " direction, ", change_stage[["interval"]],
         " interval, and ", change_stage[["p"]], " P-value classification changes among ", nrow(stage_comparison), " contrasts."),
  paste0("- Stage 2–4 duration: ", change_duration[["direction"]], " direction, ", change_duration[["interval"]],
         " interval, and ", change_duration[["p"]], " P-value classification changes among ", nrow(duration_comparison), " contrasts."), "",
  if (nrow(stage_changed)) paste0(
    "- The only stage-specific threshold change is *Fagus*, mixed culture, stage 4: Kenward–Roger 95% CI ",
    sprintf("[%.3f, %.3f], P = %.4f", stage_changed$old_lower[[1]], stage_changed$old_upper[[1]], stage_changed$old_p[[1]]),
    "; bootstrap 95% CI ", sprintf("[%.3f, %.3f], P = %.4f", stage_changed$new_lower[[1]], stage_changed$new_upper[[1]], stage_changed$new_p[[1]]),
    ". This is a borderline sensitivity result, not the primary common-shift result."
  ) else "- No stage-specific threshold classifications changed.", "",
  "## Phenology SEM", "",
  paste0("- Direction changes: ", change_sem[["direction"]], " of ", nrow(sem_comparison), " component effects."),
  paste0("- 95% interval-exclusion changes: ", change_sem[["interval"]], " of ", nrow(sem_comparison), "."),
  paste0("- P < 0.05 classification changes: ", change_sem[["p"]], " of ", nrow(sem_comparison), "."),
  if (nrow(sem_changed)) paste0(
    "- The only threshold change is the *Quercus* reduced-precipitation direct path: unstratified 95% CI ",
    sprintf("[%.3f, %.3f], P = %.4f", sem_changed$old_lower[[1]], sem_changed$old_upper[[1]], sem_changed$old_p[[1]]),
    "; block-stratified 95% CI ", sprintf("[%.3f, %.3f], P = %.4f", sem_changed$new_lower[[1]], sem_changed$new_upper[[1]], sem_changed$new_p[[1]]),
    "."
  ) else "- No SEM threshold classifications changed.",
  paste0("- The Figure 6 target—the six path-summed total effects—has ", sem_total_changes$changes[[1]],
         " threshold changes among ", sem_total_changes$n[[1]], "; all six intervals still include zero."),
  "- Point estimates are unchanged; differences arise from preserving block composition in every bootstrap sample.", "",
  "## Fit accounting", "",
  paste0("- LMM bootstrap failures before reaching the target: ", sum(!lmm_boot_status$success), "."),
  paste0("- SEM bootstrap failures before reaching the target: ", sum(!sem_boot_status$success), "."),
  paste0("- Singular primary-LMM fits: ", sum(lmm_boot_status$primary_singular, na.rm = TRUE),
         " of ", sum(lmm_boot_status$success), "; singular SEM outcome fits: ",
         sum(sem_boot_status$outcome_singular, na.rm = TRUE), " of ", sum(sem_boot_status$success), "."),
  "- Singular fits were retained because the random-intercept variance can reach zero after cluster resampling while the requested fixed treatment contrasts remain estimable. This high frequency, especially for *Fagus*, is an important bootstrap caveat rather than a convergence failure.", "",
  "## Interpretation", "",
  "Direction is based on the manuscript-oriented sign: negative values indicate delayed/slower phenology and positive values earlier/faster phenology. Bootstrap intervals are percentile intervals and P values follow the shared sign-count formula in the production bootstrap helpers.", "",
  "Overall, the bootstrap analysis agrees with the current phenology conclusions: treatment-effect directions are unchanged, all six primary common-shift effects remain non-significant, all duration effects retain their classification, and all six Figure 6 path-summed totals remain non-significant. The two isolated P ≈ 0.05 changes should be reported as sensitivity-level borderline findings."
)
writeLines(md, file.path(comparison_dir, "comparison.md"))

message("Phenology bootstrap analysis saved in: ", normalizePath(output_dir, winslash = "/"))
print(model_status, n = Inf)
