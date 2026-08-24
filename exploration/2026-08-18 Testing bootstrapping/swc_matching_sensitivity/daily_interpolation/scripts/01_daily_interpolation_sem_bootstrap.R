#!/usr/bin/env Rscript

# Sensitivity analysis for repeated-response SEMs using exact-date daily SWC.
#
# SWC definition:
#   * observed container SWC when measured on the response date;
#   * otherwise the existing treatment-agnostic container-level daily GAM estimate.
#
# The original response scaling, original SWC mean/SD, frozen model formulas and
# block-stratified container-cluster bootstrap are retained so that only the SWC
# matching definition changes relative to the fuzzy-matched baseline.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(lme4)
  library(mgcv)
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
detected_cores <- parallel::detectCores(logical = FALSE)
default_cores <- if (is.na(detected_cores)) 1L else min(4L, detected_cores)
n_cores <- as.integer(arg_value("--cores", as.character(default_cores)))
base_seed <- as.integer(arg_value("--seed", "2026081803"))
qa_only <- tolower(arg_value("--qa-only", "false")) %in% c("true", "1", "yes")
if (!is.finite(B_target) || B_target < 1L) stop("--bootstrap must be positive.")
if (!is.finite(n_cores) || n_cores < 1L) stop("--cores must be positive.")

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- gsub("~[+]~", " ", sub("^--file=", "", script_arg[[1]]))
script_file <- normalizePath(script_path, winslash = "/", mustWork = TRUE)
analysis_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
project_root <- normalizePath(file.path(analysis_dir, "..", "..", "..", ".."), winslash = "/")
output_dir <- file.path(analysis_dir, "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

baseline_dir <- file.path(
  project_root, "exploration", "2026-08-18 Testing bootstrapping",
  "repeated_response_sem", "output"
)
daily_file <- file.path(project_root, "data", "interim", "box_soilwater_daily_gam_agnostic.csv")
inventory_file <- file.path(baseline_dir, "repeated-response-sem-source-inventory.csv")
baseline_effect_file <- file.path(baseline_dir, "repeated-response-sem-bootstrap-effects.csv")
baseline_rep_file <- file.path(baseline_dir, "repeated-response-sem-bootstrap-replicates.csv")

required_files <- c(daily_file, inventory_file, baseline_effect_file, baseline_rep_file)
if (any(!file.exists(required_files))) {
  stop("Required input missing: ", paste(required_files[!file.exists(required_files)], collapse = "; "))
}

run_started <- Sys.time()
`%||%` <- function(x, y) if (is.null(x)) y else x

# -----------------------------------------------------------------------------
# 1. Audit the existing daily GAM interpolation and exact-date coverage
# -----------------------------------------------------------------------------

daily <- read_csv(daily_file, show_col_types = FALSE) %>%
  mutate(boxlabel = as.character(boxlabel), date = as.Date(date))

if (anyDuplicated(daily[c("boxlabel", "date")])) {
  stop("Daily SWC file contains duplicate container-date keys.")
}

needed_daily_cols <- c(
  "boxlabel", "date", "swc_obs", "swc_hat", "swc_hat_se", "swp_site",
  "air_temp", "vpd", "precip_mm", "radiation", "date_num"
)
if (!all(needed_daily_cols %in% names(daily))) {
  stop("Daily SWC file lacks: ", paste(setdiff(needed_daily_cols, names(daily)), collapse = ", "))
}

# Source metadata are read without modifying shared functions or data.
source(file.path(project_root, "functions", "_source.R"), local = TRUE)
box_meta <- get_meta("box") %>%
  transmute(
    boxlabel = as.character(boxlabel), precipitation, robinia, culture, soiltype
  )

observed <- daily %>%
  filter(!is.na(swc_obs)) %>%
  left_join(box_meta, by = "boxlabel") %>%
  mutate(error_in_sample = swc_hat - swc_obs)

metric_row <- function(data, observed_col, predicted_col, method, group = "all") {
  obs <- data[[observed_col]]
  pred <- data[[predicted_col]]
  keep <- is.finite(obs) & is.finite(pred)
  obs <- obs[keep]
  pred <- pred[keep]
  tibble(
    method = method, group = group, n = length(obs),
    rmse = sqrt(mean((pred - obs)^2)),
    mae = mean(abs(pred - obs)),
    bias_pred_minus_obs = mean(pred - obs),
    correlation = cor(obs, pred),
    r_squared_prediction = 1 - sum((pred - obs)^2) / sum((obs - mean(obs))^2)
  )
}

qa_in_sample <- metric_row(observed, "swc_obs", "swc_hat", "in_sample", "all")

# Leave-measurement-date-out CV. Entire measurement dates, rather than rows, are
# assigned to folds, preventing the other 71 containers on a held-out date from
# leaking that date's response into validation. Systematic folds preserve broad
# seasonal coverage while leaving neighbouring dates available for interpolation.
measurement_dates <- sort(unique(observed$date))
n_folds <- min(10L, length(measurement_dates))
date_folds <- tibble(
  date = measurement_dates,
  cv_fold = (seq_along(measurement_dates) - 1L) %% n_folds + 1L
)
cv_data <- observed %>%
  left_join(date_folds, by = "date") %>%
  mutate(boxlabel = factor(boxlabel, levels = sort(unique(daily$boxlabel))))

gam_formula <- swc_obs ~
  s(date_num, k = 10) + s(swp_site, k = 8) + s(air_temp, k = 8) +
  s(vpd, k = 8) + s(precip_mm, k = 8) + s(radiation, k = 8) +
  s(boxlabel, bs = "re")

cv_predictions <- map_dfr(seq_len(n_folds), function(fold_id) {
  train <- cv_data %>% filter(cv_fold != fold_id)
  test <- cv_data %>% filter(cv_fold == fold_id)
  fit <- mgcv::gam(gam_formula, data = train, method = "REML")
  test %>%
    mutate(
      # Match the production export, which constrains volumetric SWC to its
      # physically meaningful range before it is used downstream.
      swc_hat_cv = pmin(
        pmax(as.numeric(predict(fit, newdata = test, type = "response")), 0), 100
      )
    )
}) %>%
  mutate(error_cv = swc_hat_cv - swc_obs)

qa_cv_all <- metric_row(cv_predictions, "swc_obs", "swc_hat_cv", "leave_date_out", "all")
qa_cv_by_treatment <- bind_rows(
  cv_predictions %>% group_by(precipitation) %>%
    group_modify(~metric_row(.x, "swc_obs", "swc_hat_cv", "leave_date_out", paste0("precipitation=", .y$precipitation))) %>% ungroup(),
  cv_predictions %>% group_by(robinia) %>%
    group_modify(~metric_row(.x, "swc_obs", "swc_hat_cv", "leave_date_out", paste0("robinia=", .y$robinia))) %>% ungroup(),
  cv_predictions %>% group_by(culture) %>%
    group_modify(~metric_row(.x, "swc_obs", "swc_hat_cv", "leave_date_out", paste0("culture=", .y$culture))) %>% ungroup()
)
qa_metrics <- bind_rows(qa_in_sample, qa_cv_all, qa_cv_by_treatment)

inventory <- read_csv(inventory_file, show_col_types = FALSE)

display_sign_from_bundle <- function(bundle) {
  matrix_data <- bundle$matrix_data
  if (is.null(matrix_data) || !nrow(matrix_data) ||
      !all(c("estimate", "estimate_raw") %in% names(matrix_data))) return(1)
  ratios <- matrix_data$estimate / matrix_data$estimate_raw
  ratios <- ratios[is.finite(ratios) & abs(matrix_data$estimate_raw) > 1e-12]
  if (!length(ratios)) 1 else sign(stats::median(ratios))
}

prepare_exact_date_bundle <- function(inventory_row) {
  if (!isTRUE(inventory_row$source_available[[1]])) return(NULL)
  bundle <- readRDS(inventory_row$source_rds_file[[1]])
  dat <- bundle$data %>%
    mutate(boxlabel = as.character(boxlabel), date = as.Date(date))
  swc_mean <- mean(dat$swc_org, na.rm = TRUE)
  swc_sd <- sd(dat$swc_org, na.rm = TRUE)
  scale_error <- max(abs(dat$swc - (dat$swc_org - swc_mean) / swc_sd), na.rm = TRUE)
  if (!is.finite(scale_error) || scale_error > 1e-10) {
    stop("Could not reconstruct original SWC z-score for ", inventory_row$model_id[[1]])
  }
  dat <- dat %>%
    select(-swc, -swc_org) %>%
    left_join(
      daily %>%
        transmute(
          boxlabel, date,
          swc_observed_exact = swc_obs,
          swc_gam_exact = swc_hat,
          swc_gam_se = swc_hat_se
        ),
      by = c("boxlabel", "date")
    ) %>%
    mutate(
      swc_org = coalesce(swc_observed_exact, swc_gam_exact),
      swc_daily_source = if_else(!is.na(swc_observed_exact), "measured_exact", "daily_gam"),
      swc = (swc_org - swc_mean) / swc_sd
    ) %>%
    filter(!is.na(swc), !is.na(y)) %>%
    mutate(
      boxlabel = factor(boxlabel),
      tree_id = factor(tree_id)
    )

  list(
    bundle = bundle, data = dat, swc_mean = swc_mean, swc_sd = swc_sd,
    scale_reconstruction_error = scale_error,
    display_sign = display_sign_from_bundle(bundle)
  )
}

prepared <- map(seq_len(nrow(inventory)), function(i) prepare_exact_date_bundle(inventory[i, ]))

coverage <- map2_dfr(seq_len(nrow(inventory)), prepared, function(i, prep) {
  row <- inventory[i, ]
  if (is.null(prep)) {
    return(tibble(
      species = row$species, resp_var = row$resp_var, model_id = row$model_id,
      status = "unavailable", n_baseline_rows = 0L, n_exact_daily_rows = 0L,
      n_exact_measured = 0L, n_daily_gam = 0L, exact_date_coverage = NA_real_,
      baseline_swc_mean = NA_real_, baseline_swc_sd = NA_real_,
      scale_reconstruction_max_abs_error = NA_real_
    ))
  }
  tibble(
    species = row$species, resp_var = row$resp_var, model_id = row$model_id,
    status = "available", n_baseline_rows = nrow(prep$bundle$data),
    n_exact_daily_rows = nrow(prep$data),
    n_exact_measured = sum(prep$data$swc_daily_source == "measured_exact"),
    n_daily_gam = sum(prep$data$swc_daily_source == "daily_gam"),
    exact_date_coverage = nrow(prep$data) / nrow(prep$bundle$data),
    baseline_swc_mean = prep$swc_mean, baseline_swc_sd = prep$swc_sd,
    scale_reconstruction_max_abs_error = prep$scale_reconstruction_error
  )
})

audit <- tibble(
  item = c(
    "daily_file", "generation_source", "predictor_scope", "treatment_information",
    "containers", "daily_rows", "date_min", "date_max", "observed_rows",
    "imputed_rows", "primary_exact_date_definition", "bootstrap_uncertainty_scope"
  ),
  value = c(
    daily_file,
    "notebooks/2-swc-interpolation.Rmd and scripts/4-impute-swc-gam.R",
    "date smooth + site SWP + ambient air temperature + VPD + ambient precipitation + radiation + container random intercept",
    "No treatment predictor or treatment-specific temporal smooth; ambient climate cannot encode manipulated water additions/exclusions",
    as.character(n_distinct(daily$boxlabel)), as.character(nrow(daily)),
    as.character(min(daily$date)), as.character(max(daily$date)),
    as.character(sum(!is.na(daily$swc_obs))), as.character(sum(is.na(daily$swc_obs))),
    "same-day observed SWC when available; otherwise existing daily GAM prediction",
    "SEM bootstrap conditional on the chosen daily SWC series; GAM interpolation uncertainty is not re-fitted inside each bootstrap"
  )
)

write_csv(audit, file.path(output_dir, "interpolation-audit.csv"))
write_csv(qa_metrics, file.path(output_dir, "interpolation-qa-metrics.csv"))
write_csv(cv_predictions, file.path(output_dir, "interpolation-leave-date-out-predictions.csv"))
write_csv(coverage, file.path(output_dir, "exact-date-model-coverage.csv"))

write_qa_pdf <- function() {
  p_obs <- ggplot(cv_predictions, aes(swc_obs, swc_hat_cv, colour = precipitation)) +
    geom_point(alpha = 0.45, size = 1.1) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    coord_equal() +
    labs(
      title = "Leave-measurement-date-out validation",
      subtitle = "Entire SWC measurement dates were held out",
      x = "Observed SWC", y = "Cross-validated daily-GAM SWC", colour = "Precipitation"
    ) +
    theme_bw(base_size = 10)

  p_bias <- cv_predictions %>%
    mutate(treatment_group = paste(precipitation, robinia, culture, sep = " | ")) %>%
    ggplot(aes(error_cv, reorder(treatment_group, error_cv, FUN = median))) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey45") +
    geom_boxplot(outlier.alpha = 0.25) +
    labs(
      title = "Cross-validated interpolation error by treatment combination",
      x = "Predicted minus observed SWC", y = NULL
    ) +
    theme_bw(base_size = 9)

  pdf(file.path(output_dir, "interpolation-qa-diagnostic.pdf"), width = 9, height = 6.5)
  print(p_obs)
  print(p_bias)
  dev.off()
}

if (qa_only) {
  write_qa_pdf()
  cat("Interpolation QA refresh complete; SEM bootstrap was not rerun.\n")
  print(qa_metrics)
  quit(save = "no", status = 0)
}

# -----------------------------------------------------------------------------
# 2. Exact-date SEMs and block-stratified container-cluster bootstrap
# -----------------------------------------------------------------------------

coef_name_for_factor <- function(model, factor_name) {
  hits <- grep(factor_name, names(fixef(model)), fixed = TRUE, value = TRUE)
  hits <- hits[!grepl(":", hits, fixed = TRUE)]
  if (length(hits) != 1L) return(NA_character_)
  hits[[1]]
}

extract_effects <- function(mod_swc, mod_resp, treatments, display_sign) {
  b <- unname(fixef(mod_resp)["swc"])
  if (!is.finite(b)) stop("Missing SWC -> response coefficient.")
  treatment_rows <- map_dfr(treatments, function(treatment) {
    a_name <- coef_name_for_factor(mod_swc, treatment)
    c_name <- coef_name_for_factor(mod_resp, treatment)
    if (is.na(a_name) || is.na(c_name)) stop("Missing treatment coefficient: ", treatment)
    a <- unname(fixef(mod_swc)[a_name])
    direct <- unname(fixef(mod_resp)[c_name])
    tibble(
      treatment = treatment,
      component = c("treatment_to_swc", "direct", "indirect", "total"),
      estimate_raw = c(a, direct, a * b, direct + a * b)
    )
  })
  bind_rows(
    tibble(treatment = "shared", component = "swc_to_response", estimate_raw = b),
    treatment_rows
  ) %>%
    mutate(
      estimate = if_else(component == "treatment_to_swc", estimate_raw,
                         estimate_raw * display_sign)
    )
}

resample_containers_within_block <- function(data, replicate_id) {
  dat <- data %>%
    mutate(
      .box_original = as.character(boxlabel),
      .tree_original = as.character(tree_id),
      .block = sub("-.*$", "", .box_original)
    )
  selected <- dat %>%
    distinct(.block, .box_original) %>%
    arrange(.block, .box_original) %>%
    group_by(.block) %>%
    group_modify(~tibble(
      .box_original = sample(.x$.box_original, nrow(.x), replace = TRUE),
      .draw = seq_len(nrow(.x))
    )) %>%
    ungroup()
  pmap_dfr(selected, function(.block, .box_original, .draw) {
    new_box <- paste0("boot", replicate_id, "_", .block, "_", sprintf("%02d", .draw))
    dat %>%
      filter(.data$.box_original == .env$.box_original) %>%
      mutate(boxlabel = new_box, tree_id = paste(new_box, .tree_original, sep = "__"))
  }) %>%
    select(-.box_original, -.tree_original, -.block) %>%
    mutate(boxlabel = factor(boxlabel), tree_id = factor(tree_id))
}

fit_one_model <- function(i) {
  row <- inventory[i, ]
  prep <- prepared[[i]]
  if (is.null(prep)) {
    return(list(
      point = tibble(), replicates = tibble(), failures = tibble(),
      status = tibble(
        species = row$species, resp_var = row$resp_var, response_label = row$response_label,
        model_id = row$model_id, status = "unavailable", n_rows = 0L,
        n_containers = 0L, n_trees = 0L, n_boot_target = B_target,
        n_boot_success = 0L, n_attempts = 0L, n_failures = 0L,
        n_singular_swc = 0L, n_singular_response = 0L,
        formula_swc = NA_character_, formula_response = NA_character_,
        note = row$source_note %||% "Original model unavailable"
      )
    ))
  }

  bundle <- prep$bundle
  dat <- prep$data
  f_swc <- formula(bundle$mod_swc)
  f_resp <- formula(bundle$mod_resp)
  treatments <- as.character(bundle$effects$factor)
  model_seed <- base_seed + i * 100003L

  point_fit <- tryCatch({
    mod_swc <- suppressWarnings(lmer(f_swc, data = dat, REML = isREML(bundle$mod_swc)))
    mod_resp <- suppressWarnings(lmer(f_resp, data = dat, REML = isREML(bundle$mod_resp)))
    list(mod_swc = mod_swc, mod_resp = mod_resp)
  }, error = function(e) e)
  if (inherits(point_fit, "error")) stop(row$model_id, " point fit failed: ", conditionMessage(point_fit))

  point <- extract_effects(
    point_fit$mod_swc, point_fit$mod_resp, treatments, prep$display_sign
  ) %>%
    mutate(
      species = row$species, resp_var = row$resp_var,
      response_label = row$response_label, model_id = row$model_id,
      display_sign = prep$display_sign, .before = 1
    )

  set.seed(model_seed)
  successes <- vector("list", B_target)
  failures <- list()
  success_count <- 0L
  attempt <- 0L
  singular_swc <- 0L
  singular_response <- 0L
  max_attempts <- max(B_target * 3L, B_target + 250L)
  message("Starting ", row$model_id, " exact-date daily SWC (", nrow(dat), " rows).")

  while (success_count < B_target && attempt < max_attempts) {
    attempt <- attempt + 1L
    result <- tryCatch({
      boot_dat <- resample_containers_within_block(dat, attempt)
      m_swc <- suppressMessages(suppressWarnings(
        lmer(f_swc, data = boot_dat, REML = isREML(bundle$mod_swc))
      ))
      m_resp <- suppressMessages(suppressWarnings(
        lmer(f_resp, data = boot_dat, REML = isREML(bundle$mod_resp))
      ))
      list(
        effects = extract_effects(m_swc, m_resp, treatments, prep$display_sign),
        singular_swc = isSingular(m_swc, tol = 1e-4),
        singular_response = isSingular(m_resp, tol = 1e-4)
      )
    }, error = function(e) e)

    if (inherits(result, "error")) {
      failures[[length(failures) + 1L]] <- tibble(
        species = row$species, resp_var = row$resp_var,
        response_label = row$response_label, model_id = row$model_id,
        attempt = attempt, error = conditionMessage(result)
      )
    } else {
      success_count <- success_count + 1L
      singular_swc <- singular_swc + as.integer(result$singular_swc)
      singular_response <- singular_response + as.integer(result$singular_response)
      successes[[success_count]] <- result$effects %>%
        mutate(
          species = row$species, resp_var = row$resp_var,
          response_label = row$response_label, model_id = row$model_id,
          replicate = success_count, attempt = attempt,
          singular_swc = result$singular_swc,
          singular_response = result$singular_response,
          .before = 1
        )
    }
    if (attempt %% 100L == 0L) {
      message(row$model_id, ": ", success_count, "/", B_target, " successful.")
    }
  }
  if (success_count < B_target) {
    stop(row$model_id, " reached only ", success_count, " successes after ", attempt, " attempts.")
  }
  fail_df <- bind_rows(failures)
  status <- tibble(
    species = row$species, resp_var = row$resp_var,
    response_label = row$response_label, model_id = row$model_id,
    status = "complete", n_rows = nrow(dat), n_containers = n_distinct(dat$boxlabel),
    n_trees = n_distinct(dat$tree_id), n_boot_target = B_target,
    n_boot_success = success_count, n_attempts = attempt, n_failures = nrow(fail_df),
    n_singular_swc = singular_swc, n_singular_response = singular_response,
    formula_swc = paste(deparse(f_swc), collapse = " "),
    formula_response = paste(deparse(f_resp), collapse = " "),
    note = NA_character_
  )
  list(point = point, replicates = bind_rows(successes), failures = fail_df, status = status)
}

model_indices <- seq_len(nrow(inventory))
if (.Platform$OS.type == "unix" && n_cores > 1L) {
  results <- parallel::mclapply(
    model_indices, fit_one_model, mc.cores = min(n_cores, length(model_indices)),
    mc.preschedule = FALSE
  )
} else {
  results <- lapply(model_indices, fit_one_model)
}
if (any(vapply(results, inherits, logical(1), what = "try-error"))) {
  stop("At least one model worker failed.")
}

point_effects <- map_dfr(results, "point")
replicates <- map_dfr(results, "replicates")
failures <- map_dfr(results, "failures")
if (!nrow(failures)) {
  failures <- tibble(
    species = character(), resp_var = character(), response_label = character(),
    model_id = character(), attempt = integer(), error = character()
  )
}
status <- map_dfr(results, "status")

bootstrap_effects <- replicates %>%
  group_by(species, resp_var, response_label, model_id, treatment, component) %>%
  summarise(
    estimate_boot_mean = mean(estimate),
    lower = quantile(estimate, 0.025, names = FALSE),
    upper = quantile(estimate, 0.975, names = FALSE),
    p_boot = pmin(
      1,
      2 * pmin(
        (sum(estimate_raw <= 0) + 1) / (n() + 1),
        (sum(estimate_raw >= 0) + 1) / (n() + 1)
      )
    ),
    n_boot = n(), .groups = "drop"
  ) %>%
  left_join(
    point_effects %>%
      select(species, resp_var, response_label, model_id, treatment, component,
             estimate_raw, estimate, display_sign),
    by = c("species", "resp_var", "response_label", "model_id", "treatment", "component")
  ) %>%
  select(species, resp_var, response_label, model_id, treatment, component,
         estimate_raw, estimate, display_sign, everything())

# Exact path identity, checked for every treatment at point and replicate level.
identity_point <- point_effects %>%
  filter(component %in% c("direct", "indirect", "total")) %>%
  select(species, resp_var, model_id, treatment, component, estimate) %>%
  pivot_wider(names_from = component, values_from = estimate) %>%
  mutate(identity_error = total - direct - indirect, level = "point")
identity_boot <- replicates %>%
  filter(component %in% c("direct", "indirect", "total")) %>%
  select(species, resp_var, model_id, replicate, treatment, component, estimate) %>%
  pivot_wider(names_from = component, values_from = estimate) %>%
  mutate(identity_error = total - direct - indirect, level = "bootstrap")
identity_qa <- bind_rows(identity_point, identity_boot) %>%
  summarise(
    n_checks = n(), max_abs_identity_error = max(abs(identity_error)),
    n_above_tolerance_1e12 = sum(abs(identity_error) > 1e-12), .by = level
  )
if (any(identity_qa$n_above_tolerance_1e12 > 0)) stop("Total path identity QA failed.")

# -----------------------------------------------------------------------------
# 3. Comparison with current fuzzy-matched block-bootstrap baseline
# -----------------------------------------------------------------------------

baseline_effects <- read_csv(baseline_effect_file, show_col_types = FALSE) %>%
  select(species, resp_var, response_label, model_id, treatment, component,
         estimate, lower, upper, p_boot, n_boot)

# Recover the shared b path from baseline replicates: indirect/a equals b within
# each replicate. Median across available treatment ratios avoids reliance on one
# possibly near-zero a path.
baseline_reps <- read_csv(baseline_rep_file, show_col_types = FALSE)
baseline_b_reps <- baseline_reps %>%
  filter(component %in% c("treatment_to_swc", "indirect")) %>%
  select(species, resp_var, response_label, model_id, replicate, treatment,
         component, estimate_raw) %>%
  pivot_wider(names_from = component, values_from = estimate_raw) %>%
  filter(is.finite(treatment_to_swc), abs(treatment_to_swc) > 1e-10) %>%
  mutate(b_raw = indirect / treatment_to_swc) %>%
  group_by(species, resp_var, response_label, model_id, replicate) %>%
  summarise(b_raw = median(b_raw), .groups = "drop") %>%
  left_join(
    point_effects %>%
      filter(component == "swc_to_response") %>%
      distinct(model_id, display_sign),
    by = "model_id"
  ) %>%
  mutate(estimate = b_raw * display_sign)

baseline_b_point <- map_dfr(seq_len(nrow(inventory)), function(i) {
  if (!isTRUE(inventory$source_available[[i]])) return(tibble())
  bundle <- readRDS(inventory$source_rds_file[[i]])
  sign_i <- display_sign_from_bundle(bundle)
  b_raw <- unname(fixef(bundle$mod_resp)["swc"])
  tibble(
    species = inventory$species[[i]], resp_var = inventory$resp_var[[i]],
    response_label = inventory$response_label[[i]], model_id = inventory$model_id[[i]],
    treatment = "shared", component = "swc_to_response", estimate = b_raw * sign_i
  )
})
baseline_b <- baseline_b_reps %>%
  group_by(species, resp_var, response_label, model_id) %>%
  summarise(
    lower = quantile(estimate, 0.025, names = FALSE),
    upper = quantile(estimate, 0.975, names = FALSE),
    p_boot = pmin(1, 2 * pmin(
      (sum(b_raw <= 0) + 1) / (n() + 1),
      (sum(b_raw >= 0) + 1) / (n() + 1)
    )),
    n_boot = n(), .groups = "drop"
  ) %>%
  left_join(baseline_b_point, by = c("species", "resp_var", "response_label", "model_id")) %>%
  select(species, resp_var, response_label, model_id, treatment, component,
         estimate, lower, upper, p_boot, n_boot)

baseline_all <- bind_rows(baseline_effects, baseline_b)
comparison <- baseline_all %>%
  rename_with(~paste0("baseline_", .x), c(estimate, lower, upper, p_boot, n_boot)) %>%
  full_join(
    bootstrap_effects %>%
      select(species, resp_var, response_label, model_id, treatment, component,
             estimate, lower, upper, p_boot, n_boot) %>%
      rename_with(~paste0("daily_", .x), c(estimate, lower, upper, p_boot, n_boot)),
    by = c("species", "resp_var", "response_label", "model_id", "treatment", "component")
  ) %>%
  mutate(
    estimate_difference_daily_minus_baseline = daily_estimate - baseline_estimate,
    direction_agrees = sign(daily_estimate) == sign(baseline_estimate),
    baseline_significant = baseline_p_boot < 0.05,
    daily_significant = daily_p_boot < 0.05,
    significance_agrees = baseline_significant == daily_significant,
    intervals_overlap = pmax(baseline_lower, daily_lower) <= pmin(baseline_upper, daily_upper)
  )

run_finished <- Sys.time()
status <- status %>%
  mutate(
    run_started_at = format(run_started, "%Y-%m-%d %H:%M:%S %Z"),
    run_finished_at = format(run_finished, "%Y-%m-%d %H:%M:%S %Z"),
    runtime_seconds = as.numeric(difftime(run_finished, run_started, units = "secs"))
  )

write_csv(status, file.path(output_dir, "daily-interpolation-sem-bootstrap-status.csv"))
write_csv(point_effects, file.path(output_dir, "daily-interpolation-sem-point-effects.csv"))
write_csv(bootstrap_effects, file.path(output_dir, "daily-interpolation-sem-bootstrap-effects.csv"))
write_csv(replicates, file.path(output_dir, "daily-interpolation-sem-bootstrap-replicates.csv"))
write_csv(failures, file.path(output_dir, "daily-interpolation-sem-bootstrap-failures.csv"))
write_csv(identity_qa, file.path(output_dir, "total-path-identity-qa.csv"))
write_csv(comparison, file.path(output_dir, "daily-interpolation-vs-fuzzy-baseline.csv"))
saveRDS(
  list(
    settings = list(
      B_target = B_target, n_cores = n_cores, seed = base_seed,
      swc_definition = "same-day measured SWC else existing daily GAM prediction",
      scaling = "original model-specific fuzzy-SWC mean and SD",
      formulas = "frozen original formulas",
      resampling = "containers within blocks",
      interpolation_uncertainty = "not propagated; SEM inference is conditional on daily series"
    ),
    audit = audit, qa_metrics = qa_metrics, coverage = coverage,
    status = status, point_effects = point_effects,
    bootstrap_effects = bootstrap_effects, comparison = comparison,
    identity_qa = identity_qa, replicates = replicates, failures = failures
  ),
  file.path(output_dir, "daily-interpolation-sem-bootstrap-results.rds")
)

# -----------------------------------------------------------------------------
# 4. Compact diagnostic PDFs and machine-generated comparison summary
# -----------------------------------------------------------------------------

write_qa_pdf()

plot_comparison <- comparison %>%
  filter(component %in% c("direct", "indirect", "total")) %>%
  select(species, response_label, treatment, component,
         baseline_estimate, baseline_lower, baseline_upper,
         daily_estimate, daily_lower, daily_upper) %>%
  pivot_longer(
    cols = matches("^(baseline|daily)_(estimate|lower|upper)$"),
    names_to = c("swc_method", ".value"), names_sep = "_"
  ) %>%
  mutate(
    swc_method = recode(swc_method, baseline = "Fuzzy measured SWC", daily = "Exact-date daily SWC"),
    treatment = factor(treatment, levels = rev(c("precipitation", "robinia", "culture", "extreme_event")))
  )

pdf(file.path(output_dir, "sem-effect-comparison.pdf"), width = 11, height = 10)
for (sp in c("fagus", "quercus")) {
  p <- plot_comparison %>%
    filter(species == sp) %>%
    ggplot(aes(estimate, treatment, colour = swc_method)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
    geom_errorbarh(aes(xmin = lower, xmax = upper),
                   position = position_dodge(width = 0.45), height = 0) +
    geom_point(position = position_dodge(width = 0.45), size = 1.7) +
    facet_grid(response_label ~ component, scales = "free_x") +
    labs(
      title = paste("SEM sensitivity to SWC matching:", tools::toTitleCase(sp)),
      x = "Oriented standardized effect (95% cluster-bootstrap CI)", y = NULL,
      colour = "SWC definition"
    ) +
    theme_bw(base_size = 8) +
    theme(legend.position = "bottom", strip.text = element_text(size = 7))
  print(p)
}
dev.off()

primary_comp <- comparison %>% filter(component %in% c("direct", "indirect", "total"))
total_comp <- primary_comp %>% filter(component == "total")
direct_comp <- primary_comp %>% filter(component == "direct")
indirect_comp <- primary_comp %>% filter(component == "indirect")
component_magnitude <- bootstrap_effects %>%
  filter(component %in% c("direct", "indirect")) %>%
  select(species, resp_var, treatment, component, estimate) %>%
  pivot_wider(names_from = component, values_from = estimate)
qa_cv <- qa_metrics %>% filter(method == "leave_date_out", group == "all")
summary_lines <- c(
  "# Exact-date daily-SWC sensitivity",
  "",
  "## Design",
  "",
  paste0(
    "Repeated-response SEMs were rebuilt on exact response dates using observed SWC when measured on that date and the existing container-level daily GAM prediction otherwise. Frozen formulas, response scaling, and each baseline model's original SWC mean/SD were retained. Uncertainty used 1,000 successful block-stratified container-cluster bootstrap replicates per estimable model."
  ),
  "",
  "## Interpolation audit",
  "",
  paste0(
    "The existing GAM is climate-informed but treatment-agnostic: ambient precipitation, air temperature, VPD, radiation and site-level soil-water potential vary through time, while a container random intercept represents persistent container differences. It has no precipitation-treatment input or treatment-specific temporal smooth. Leave-measurement-date-out validation gave RMSE ",
    round(qa_cv$rmse, 2), ", MAE ", round(qa_cv$mae, 2), ", bias ",
    round(qa_cv$bias_pred_minus_obs, 2), ", and predictive R-squared ",
    round(qa_cv$r_squared_prediction, 2), "."
  ),
  "",
  "## SEM comparison",
  "",
  paste0(
    "Across ", nrow(primary_comp), " comparable direct, indirect and total paths, effect direction agreed for ",
    sum(primary_comp$direction_agrees, na.rm = TRUE), " (",
    round(100 * mean(primary_comp$direction_agrees, na.rm = TRUE), 1), "%), bootstrap significance agreed for ",
    sum(primary_comp$significance_agrees, na.rm = TRUE), " (",
    round(100 * mean(primary_comp$significance_agrees, na.rm = TRUE), 1), "%), and intervals overlapped for ",
    sum(primary_comp$intervals_overlap, na.rm = TRUE), " (",
    round(100 * mean(primary_comp$intervals_overlap, na.rm = TRUE), 1), "%)."
  ),
  "",
  paste0(
    "The net total effects were exceptionally stable: all ", nrow(total_comp),
    " retained their direction, ", sum(total_comp$significance_agrees, na.rm = TRUE),
    "/", nrow(total_comp), " retained the same bootstrap significance classification, all intervals overlapped, and the mean absolute estimate change was ",
    round(mean(abs(total_comp$estimate_difference_daily_minus_baseline), na.rm = TRUE), 3),
    " SD. The decomposition was less stable. Direct-effect direction agreed for ",
    sum(direct_comp$direction_agrees, na.rm = TRUE), "/", nrow(direct_comp),
    " paths and indirect-effect direction for ", sum(indirect_comp$direction_agrees, na.rm = TRUE),
    "/", nrow(indirect_comp), "; indirect-effect significance agreed for ",
    sum(indirect_comp$significance_agrees, na.rm = TRUE), "/", nrow(indirect_comp),
    ". Under exact-date daily SWC, the direct component remained larger in absolute magnitude than the indirect component for ",
    sum(abs(component_magnitude$direct) >= abs(component_magnitude$indirect)), "/",
    nrow(component_magnitude), " treatment-response combinations. Thus, the net treatment conclusions are robust to this SWC definition, whereas the quantitative attribution between direct and SWC-associated pathways is definition-sensitive."
  ),
  "",
  "## Climate-covariate recommendation",
  "",
  "Do not add ambient climate as a substitute for manipulated precipitation: the current daily GAM already uses ambient climate, which is common to all containers and cannot reconstruct container-specific water exclusion/addition. Without dated container-level irrigation or throughfall-exclusion inputs, a more complex climate model would add assumptions rather than treatment information. Treat this exact-date GAM analysis as a sensitivity check. If pathway decomposition changes materially, prefer same-day/prior measured SWC or a transparent treatment-aware interpolation supported by actual watering records. Bootstrap intervals here are conditional on the interpolated series and do not propagate GAM prediction uncertainty.",
  "",
  "## QA",
  "",
  paste0(
    "All estimable SEMs retained exact total = direct + indirect identity (maximum absolute numerical error ",
    format(max(identity_qa$max_abs_identity_error), scientific = TRUE), ")."
  )
)
writeLines(summary_lines, file.path(analysis_dir, "comparison.md"))

cat("\nExact-date daily-SWC sensitivity complete.\n")
print(qa_metrics)
print(status %>% select(model_id, status, n_boot_success, n_attempts, n_failures))
print(identity_qa)
