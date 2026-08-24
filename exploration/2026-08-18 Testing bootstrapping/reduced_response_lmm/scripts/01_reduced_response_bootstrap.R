#!/usr/bin/env Rscript

# Reduced-form treatment models for the standardized response summary in
# Figure 6. These models use the response-model analysis samples and frozen
# covariate/random-effect structures from the repeated-response analysis, but
# omit SWC. Consequently, the treatment coefficients are estimated directly
# rather than reconstructed as SEM path sums.

suppressPackageStartupMessages({
  library(dplyr)
  library(lme4)
  library(purrr)
  library(readr)
  library(tibble)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default) {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^", flag, "="), "", hit[[1]])
}

B_target <- as.integer(arg_value("--bootstrap", "1000"))
n_cores <- as.integer(arg_value(
  "--cores",
  as.character(min(4L, parallel::detectCores(logical = FALSE)))
))
base_seed <- as.integer(arg_value("--seed", "20260824"))
if (!is.finite(B_target) || B_target < 1L) stop("--bootstrap must be positive.")
if (!is.finite(n_cores) || n_cores < 1L) stop("--cores must be positive.")

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- gsub("~[+]~", " ", sub("^--file=", "", script_arg[[1]]))
script_file <- normalizePath(script_path, winslash = "/", mustWork = TRUE)
analysis_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/", mustWork = TRUE)
project_root <- normalizePath(file.path(analysis_dir, "..", "..", ".."), winslash = "/", mustWork = TRUE)
output_dir <- file.path(analysis_dir, "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(project_root, "functions", "_source.R"))
source(file.path(project_root, "functions", "11-bootstrap-inference.R"))

paths <- alinv_ensure_bootstrap_family("repeated_sem", project_root, target = 1000L)
sem_status <- readr::read_csv(paths$repeated_sem_status, show_col_types = FALSE)

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
treatment_order <- c("precipitation", "robinia", "culture", "extreme_event")

model_grid <- tidyr::crossing(
  species = c("fagus", "quercus"),
  response_specs
) |>
  dplyr::mutate(
    model_id = paste(.data$species, .data$resp_var, sep = "__"),
    model_seed = .env$base_seed + dplyr::row_number() * 100003L
  ) |>
  dplyr::left_join(
    sem_status |>
      dplyr::select(
        .data$model_id,
        sem_status = .data$status,
        source_rds_file = .data$source_rds_file,
        source_note = .data$note
      ),
    by = "model_id"
  )

response_display_sign <- function(bundle) {
  matrix_data <- bundle$matrix_data
  if (is.null(matrix_data) || !nrow(matrix_data) ||
      !all(c("estimate", "estimate_raw") %in% names(matrix_data))) return(1)
  ratios <- matrix_data$estimate / matrix_data$estimate_raw
  ratios <- ratios[is.finite(ratios) & abs(matrix_data$estimate_raw) > 1e-12]
  if (!length(ratios)) 1 else sign(stats::median(ratios))
}

extract_treatment_effects <- function(model, treatments, display_sign) {
  purrr::map_dfr(treatments, function(treatment) {
    coef_name <- alinv_factor_main_effect_coef(model, treatment)
    if (is.na(coef_name)) return(tibble::tibble())
    estimate_raw <- unname(lme4::fixef(model)[coef_name])
    tibble::tibble(
      treatment = treatment,
      estimate_raw = estimate_raw,
      estimate = estimate_raw * display_sign
    )
  })
}

fit_one_model <- function(row) {
  species <- row$species[[1]]
  resp_var <- row$resp_var[[1]]
  response_label <- row$response_label[[1]]
  model_id <- row$model_id[[1]]
  model_seed <- row$model_seed[[1]]
  source_rds_file <- row$source_rds_file[[1]]

  unavailable <- is.na(row$sem_status[[1]]) || row$sem_status[[1]] != "complete" ||
    is.na(source_rds_file) || !file.exists(source_rds_file)
  if (unavailable) {
    return(list(
      point = tibble(), replicates = tibble(), failures = tibble(),
      status = tibble(
        species, resp_var, response_label, model_id,
        status = "unavailable", n_rows = 0L, n_containers = 0L, n_trees = 0L,
        n_boot_target = B_target, n_boot_success = 0L, n_attempts = 0L,
        n_failures = 0L, n_singular = 0L, formula_reduced = NA_character_,
        modeled_treatments = NA_character_, seed = model_seed,
        source_rds_file = source_rds_file,
        note = row$source_note[[1]]
      )
    ))
  }

  bundle <- readRDS(source_rds_file)
  data <- bundle$data
  reduced_formula <- update(stats::formula(bundle$mod_resp), . ~ . - swc)
  treatments <- intersect(treatment_order, as.character(bundle$effects$factor))
  display_sign <- response_display_sign(bundle)
  point_model <- suppressMessages(suppressWarnings(
    lme4::lmer(reduced_formula, data = data, REML = lme4::isREML(bundle$mod_resp))
  ))
  point <- extract_treatment_effects(point_model, treatments, display_sign) |>
    dplyr::mutate(
      species, resp_var, response_label, model_id,
      point_model_singular = lme4::isSingular(point_model, tol = 1e-4),
      source_rds_file = source_rds_file,
      .before = 1
    )

  set.seed(model_seed)
  successes <- vector("list", B_target)
  failures <- list()
  success_count <- 0L
  attempt <- 0L
  singular_count <- 0L
  max_attempts <- max(B_target * 3L, B_target + 250L)

  message("Starting ", model_id, " (", nrow(data), " rows; ",
          dplyr::n_distinct(data$boxlabel), " containers).")
  while (success_count < B_target && attempt < max_attempts) {
    attempt <- attempt + 1L
    result <- tryCatch({
      boot_data <- alinv_resample_containers_within_block(data, attempt)
      model <- suppressMessages(suppressWarnings(
        lme4::lmer(
          reduced_formula,
          data = boot_data,
          REML = lme4::isREML(bundle$mod_resp)
        )
      ))
      effects <- extract_treatment_effects(model, treatments, display_sign)
      if (nrow(effects) != length(treatments) || any(!is.finite(effects$estimate))) {
        stop("Incomplete or non-finite treatment coefficients.")
      }
      list(effects = effects, singular = lme4::isSingular(model, tol = 1e-4))
    }, error = function(e) e)

    if (inherits(result, "error")) {
      failures[[length(failures) + 1L]] <- tibble(
        species, resp_var, response_label, model_id, attempt,
        message = conditionMessage(result)
      )
      next
    }

    success_count <- success_count + 1L
    singular_count <- singular_count + as.integer(result$singular)
    successes[[success_count]] <- result$effects |>
      dplyr::mutate(
        species, resp_var, response_label, model_id,
        replicate = success_count,
        attempt = attempt,
        singular = result$singular,
        .before = 1
      )
  }

  if (success_count < B_target) {
    stop(model_id, " reached only ", success_count, " successful refits.")
  }

  list(
    point = point,
    replicates = dplyr::bind_rows(successes),
    failures = dplyr::bind_rows(failures),
    status = tibble(
      species, resp_var, response_label, model_id,
      status = "complete", n_rows = nrow(data),
      n_containers = dplyr::n_distinct(data$boxlabel),
      n_trees = dplyr::n_distinct(data$tree_id),
      n_boot_target = B_target, n_boot_success = success_count,
      n_attempts = attempt, n_failures = length(failures),
      n_singular = singular_count,
      formula_reduced = paste(deparse(reduced_formula), collapse = " "),
      modeled_treatments = paste(treatments, collapse = ";"),
      seed = model_seed, source_rds_file = source_rds_file,
      note = NA_character_
    )
  )
}

rows <- split(model_grid, seq_len(nrow(model_grid)))
if (.Platform$OS.type == "unix" && n_cores > 1L) {
  results <- parallel::mclapply(rows, fit_one_model, mc.cores = n_cores,
                                mc.preschedule = FALSE)
} else {
  results <- lapply(rows, fit_one_model)
}

point <- purrr::map_dfr(results, "point")
replicates <- purrr::map_dfr(results, "replicates")
failures <- purrr::map_dfr(results, "failures")
status <- purrr::map_dfr(results, "status")

effects <- replicates |>
  dplyr::group_by(
    .data$species, .data$resp_var, .data$response_label,
    .data$model_id, .data$treatment
  ) |>
  dplyr::summarise(
    estimate_boot_mean = mean(.data$estimate),
    lower = alinv_bootstrap_percentile(.data$estimate)[[1]],
    upper = alinv_bootstrap_percentile(.data$estimate)[[2]],
    p_boot = alinv_bootstrap_p(.data$estimate),
    n_boot = dplyr::n(),
    .groups = "drop"
  ) |>
  dplyr::left_join(
    point |>
      dplyr::select(
        .data$species, .data$resp_var, .data$response_label,
        .data$model_id, .data$treatment, estimate = .data$estimate,
        estimate_raw = .data$estimate_raw, .data$point_model_singular,
        .data$source_rds_file
      ),
    by = c("species", "resp_var", "response_label", "model_id", "treatment")
  ) |>
  dplyr::mutate(
    significant = is.finite(.data$p_boot) & .data$p_boot < 0.05,
    uncertainty_method = "block-stratified container-cluster bootstrap percentile",
    analysis = "reduced-form standardized response model"
  )

sem_totals <- alinv_read_repeated_sem_bootstrap_totals(project_root) |>
  dplyr::select(
    .data$species, .data$resp_var, .data$treatment,
    sem_total_estimate = .data$estimate,
    sem_total_lower = .data$lower,
    sem_total_upper = .data$upper,
    sem_total_p = .data$p_value
  )
comparison <- effects |>
  dplyr::left_join(sem_totals, by = c("species", "resp_var", "treatment")) |>
  dplyr::mutate(
    estimate_difference = .data$estimate - .data$sem_total_estimate,
    direction_agrees = sign(.data$estimate) == sign(.data$sem_total_estimate),
    significance_agrees = .data$significant == (.data$sem_total_p < 0.05)
  )

readr::write_csv(effects, file.path(output_dir, "reduced-response-bootstrap-effects.csv"))
readr::write_csv(status, file.path(output_dir, "reduced-response-bootstrap-status.csv"))
readr::write_csv(comparison, file.path(output_dir, "reduced-response-vs-sem-total.csv"))
readr::write_csv(failures, file.path(output_dir, "reduced-response-bootstrap-failures.csv"))
readr::write_csv(replicates, file.path(output_dir, "reduced-response-bootstrap-replicates.csv"))
saveRDS(
  list(point = point, effects = effects, status = status, comparison = comparison,
       replicates = replicates, failures = failures,
       settings = list(B_target = B_target, cores = n_cores, seed = base_seed)),
  file.path(output_dir, "reduced-response-bootstrap-results.rds")
)

cat("Reduced-response bootstrap complete.\n")
print(status |>
  dplyr::select(.data$model_id, .data$status, .data$n_boot_success,
                .data$n_attempts, .data$n_failures, .data$n_singular))
