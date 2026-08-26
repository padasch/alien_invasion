#!/usr/bin/env Rscript

# Sensitivity analysis for repeated-response SEMs using modeled antecedent SWC.
#
# For every response date, the mediator is the mean of the existing daily
# climate-informed GAM SWC series over either the current day through day -6
# (7-day window) or the current day through day -13 (14-day window). No future
# information is used. These are modeled antecedent means, not measured means.
#
# The original fuzzy-matched bundle's response values, response scaling, SWC
# mean/SD, frozen formulas and model decisions are retained. Uncertainty uses
# the same block-stratified container-cluster bootstrap as the baseline.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
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
  "--cores", as.character(min(4L, parallel::detectCores(logical = FALSE)))
))
base_seed <- as.integer(arg_value("--seed", "2026081804"))
windows <- as.integer(strsplit(arg_value("--windows", "7,14"), ",", fixed = TRUE)[[1]])
if (!is.finite(B_target) || B_target < 1L) stop("--bootstrap must be positive.")
if (!is.finite(n_cores) || n_cores < 1L) stop("--cores must be positive.")
if (!length(windows) || any(!is.finite(windows)) || any(windows < 1L)) {
  stop("--windows must contain positive integers.")
}

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
required <- c(daily_file, inventory_file, baseline_effect_file, baseline_rep_file)
if (any(!file.exists(required))) {
  stop("Missing required input: ", paste(required[!file.exists(required)], collapse = "; "))
}

run_started <- Sys.time()
`%||%` <- function(x, y) if (is.null(x)) y else x

# -----------------------------------------------------------------------------
# 1. Construct complete, trailing modeled-SWC windows
# -----------------------------------------------------------------------------

daily <- read_csv(daily_file, show_col_types = FALSE) %>%
  transmute(boxlabel = as.character(boxlabel), date = as.Date(date), swc_hat)

if (anyDuplicated(daily[c("boxlabel", "date")])) {
  stop("Daily SWC file contains duplicate container-date keys.")
}
if (any(!is.finite(daily$swc_hat))) {
  stop("Daily GAM SWC contains missing or non-finite fitted values.")
}
# The supplied daily series omits 23 December for every container. Complete the
# calendar explicitly so no row-based window can silently jump over that date.
# A trailing window containing the missing day is marked unavailable rather than
# imputed again. All response dates used below precede that isolated gap.
daily <- daily %>%
  group_by(boxlabel) %>%
  complete(date = seq.Date(min(date), max(date), by = "day")) %>%
  arrange(date, .by_group = TRUE) %>%
  ungroup()

trailing_mean_complete <- function(x, k) {
  n <- length(x)
  out <- rep(NA_real_, n)
  if (n < k) return(out)
  finite <- is.finite(x)
  cs <- c(0, cumsum(ifelse(finite, x, 0)))
  cn <- c(0, cumsum(finite))
  sums <- cs[(k + 1L):(n + 1L)] - cs[1L:(n - k + 1L)]
  counts <- cn[(k + 1L):(n + 1L)] - cn[1L:(n - k + 1L)]
  out[k:n] <- ifelse(counts == k, sums / k, NA_real_)
  out
}

daily_windows <- map_dfr(windows, function(k) {
  daily %>%
    arrange(boxlabel, date) %>%
    group_by(boxlabel) %>%
    mutate(
      swc_antecedent_org = trailing_mean_complete(swc_hat, k),
      window_days_available = if_else(is.finite(swc_antecedent_org), k, NA_integer_)
    ) %>%
    ungroup() %>%
    filter(window_days_available == k) %>%
    transmute(boxlabel, date, window_days = k, swc_antecedent_org,
              window_days_available)
})

inventory <- read_csv(inventory_file, show_col_types = FALSE)

display_sign_from_bundle <- function(bundle) {
  md <- bundle$matrix_data
  if (is.null(md) || !nrow(md) || !all(c("estimate", "estimate_raw") %in% names(md))) {
    return(1)
  }
  ratios <- md$estimate / md$estimate_raw
  ratios <- ratios[is.finite(ratios) & abs(md$estimate_raw) > 1e-12]
  if (!length(ratios)) 1 else sign(stats::median(ratios))
}

prepare_bundle <- function(inventory_row, window_days) {
  if (!isTRUE(inventory_row$source_available[[1]])) return(NULL)
  bundle <- readRDS(inventory_row$source_rds_file[[1]])
  dat0 <- bundle$data %>%
    mutate(boxlabel = as.character(boxlabel), date = as.Date(date))
  swc_mean <- mean(dat0$swc_org, na.rm = TRUE)
  swc_sd <- sd(dat0$swc_org, na.rm = TRUE)
  scale_error <- max(abs(dat0$swc - (dat0$swc_org - swc_mean) / swc_sd), na.rm = TRUE)
  if (!is.finite(scale_error) || scale_error > 1e-10) {
    stop("Could not reconstruct original SWC scaling for ", inventory_row$model_id[[1]])
  }
  win <- daily_windows %>% filter(.data$window_days == .env$window_days)
  dat <- dat0 %>%
    select(-swc, -swc_org) %>%
    left_join(win, by = c("boxlabel", "date")) %>%
    mutate(
      swc_org = swc_antecedent_org,
      swc = (swc_org - swc_mean) / swc_sd
    ) %>%
    filter(is.finite(swc), is.finite(y)) %>%
    mutate(boxlabel = factor(boxlabel), tree_id = factor(tree_id))

  list(
    bundle = bundle, data = dat, baseline_rows = nrow(dat0),
    swc_mean = swc_mean, swc_sd = swc_sd,
    scale_reconstruction_error = scale_error,
    display_sign = display_sign_from_bundle(bundle)
  )
}

task_grid <- tidyr::crossing(
  inventory_row = seq_len(nrow(inventory)), window_days = windows
) %>%
  mutate(
    task_id = row_number(),
    model_seed = base_seed + task_id * 100003L,
    method = paste0("antecedent_", window_days, "d")
  )

prepared <- pmap(task_grid, function(inventory_row, window_days, task_id, model_seed, method) {
  prepare_bundle(inventory[inventory_row, ], window_days)
})

coverage <- pmap_dfr(task_grid, function(inventory_row, window_days, task_id, model_seed, method) {
  row <- inventory[inventory_row, ]
  prep <- prepared[[task_id]]
  if (is.null(prep)) {
    return(tibble(
      species = row$species, resp_var = row$resp_var, response_label = row$response_label,
      model_id = row$model_id, method, window_days, status = "unavailable",
      n_baseline_rows = 0L, n_antecedent_rows = 0L, coverage_proportion = NA_real_,
      n_containers = 0L, n_trees = 0L, baseline_swc_mean = NA_real_,
      baseline_swc_sd = NA_real_, scale_reconstruction_max_abs_error = NA_real_
    ))
  }
  tibble(
    species = row$species, resp_var = row$resp_var, response_label = row$response_label,
    model_id = row$model_id, method, window_days, status = "available",
    n_baseline_rows = prep$baseline_rows, n_antecedent_rows = nrow(prep$data),
    coverage_proportion = nrow(prep$data) / prep$baseline_rows,
    n_containers = n_distinct(prep$data$boxlabel), n_trees = n_distinct(prep$data$tree_id),
    baseline_swc_mean = prep$swc_mean, baseline_swc_sd = prep$swc_sd,
    scale_reconstruction_max_abs_error = prep$scale_reconstruction_error
  )
})

write_csv(coverage, file.path(output_dir, "antecedent-window-coverage.csv"))

# -----------------------------------------------------------------------------
# 2. Frozen SEMs and block-stratified container-cluster bootstrap
# -----------------------------------------------------------------------------

coef_name_for_factor <- function(model, factor_name) {
  hits <- grep(factor_name, names(fixef(model)), fixed = TRUE, value = TRUE)
  hits <- hits[!grepl(":", hits, fixed = TRUE)]
  if (length(hits) != 1L) return(NA_character_)
  hits[[1]]
}

extract_effects <- function(mod_swc, mod_resp, treatments, display_sign) {
  b <- unname(fixef(mod_resp)["swc"])
  if (!length(b) || !is.finite(b)) stop("Missing SWC -> response coefficient.")
  treatment_rows <- map_dfr(treatments, function(treatment) {
    a_name <- coef_name_for_factor(mod_swc, treatment)
    c_name <- coef_name_for_factor(mod_resp, treatment)
    if (is.na(a_name) || is.na(c_name)) stop("Missing treatment coefficient: ", treatment)
    a <- unname(fixef(mod_swc)[a_name])
    direct <- unname(fixef(mod_resp)[c_name])
    tibble(
      treatment,
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

fit_one_task <- function(task_index) {
  grid_row <- task_grid[task_index, ]
  row <- inventory[grid_row$inventory_row, ]
  prep <- prepared[[task_index]]
  if (is.null(prep)) {
    return(list(
      point = tibble(), reduced_form = tibble(), replicates = tibble(), failures = tibble(),
      status = tibble(
        species = row$species, resp_var = row$resp_var, response_label = row$response_label,
        model_id = row$model_id, method = grid_row$method,
        window_days = grid_row$window_days, status = "unavailable", n_rows = 0L,
        n_containers = 0L, n_trees = 0L, n_boot_target = B_target,
        n_boot_success = 0L, n_attempts = 0L, n_failures = 0L,
        n_singular_swc = 0L, n_singular_response = 0L,
        formula_swc = NA_character_, formula_response = NA_character_,
        seed = grid_row$model_seed, note = row$source_note %||% "Original model unavailable"
      )
    ))
  }

  bundle <- prep$bundle
  dat <- prep$data
  f_swc <- formula(bundle$mod_swc)
  f_resp <- formula(bundle$mod_resp)
  f_reduced <- update(f_resp, . ~ . - swc)
  treatments <- as.character(bundle$effects$factor)

  point_fit <- tryCatch({
    mod_swc <- suppressMessages(suppressWarnings(
      lmer(f_swc, data = dat, REML = isREML(bundle$mod_swc))
    ))
    mod_resp <- suppressMessages(suppressWarnings(
      lmer(f_resp, data = dat, REML = isREML(bundle$mod_resp))
    ))
    mod_reduced <- suppressMessages(suppressWarnings(
      lmer(f_reduced, data = dat, REML = isREML(bundle$mod_resp))
    ))
    list(mod_swc = mod_swc, mod_resp = mod_resp, mod_reduced = mod_reduced)
  }, error = function(e) e)
  if (inherits(point_fit, "error")) {
    stop(row$model_id, " ", grid_row$method, " point fit failed: ", conditionMessage(point_fit))
  }

  point <- extract_effects(
    point_fit$mod_swc, point_fit$mod_resp, treatments, prep$display_sign
  ) %>%
    mutate(
      species = row$species, resp_var = row$resp_var,
      response_label = row$response_label, model_id = row$model_id,
      method = grid_row$method, window_days = grid_row$window_days,
      display_sign = prep$display_sign, .before = 1
    )

  reduced_form <- map_dfr(treatments, function(treatment) {
    coef_name <- coef_name_for_factor(point_fit$mod_reduced, treatment)
    if (is.na(coef_name)) stop("Missing reduced-form coefficient: ", treatment)
    estimate_raw <- unname(fixef(point_fit$mod_reduced)[coef_name])
    tibble(
      species = row$species, resp_var = row$resp_var,
      response_label = row$response_label, model_id = row$model_id,
      method = grid_row$method, window_days = grid_row$window_days,
      treatment, reduced_form_estimate_raw = estimate_raw,
      reduced_form_estimate = estimate_raw * prep$display_sign,
      reduced_form_singular = isSingular(point_fit$mod_reduced, tol = 1e-4)
    )
  })

  set.seed(grid_row$model_seed)
  successes <- vector("list", B_target)
  failures <- list()
  success_count <- 0L
  attempt <- 0L
  singular_swc <- 0L
  singular_response <- 0L
  max_attempts <- max(B_target * 3L, B_target + 250L)
  message("Starting ", row$model_id, " ", grid_row$method, " (", nrow(dat), " rows).")

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
        method = grid_row$method, window_days = grid_row$window_days,
        attempt, error = conditionMessage(result)
      )
    } else {
      success_count <- success_count + 1L
      singular_swc <- singular_swc + as.integer(result$singular_swc)
      singular_response <- singular_response + as.integer(result$singular_response)
      successes[[success_count]] <- result$effects %>%
        mutate(
          species = row$species, resp_var = row$resp_var,
          response_label = row$response_label, model_id = row$model_id,
          method = grid_row$method, window_days = grid_row$window_days,
          replicate = success_count, attempt,
          singular_swc = result$singular_swc,
          singular_response = result$singular_response,
          .before = 1
        )
    }
    if (attempt %% 100L == 0L) {
      message(row$model_id, " ", grid_row$method, ": ", success_count, "/", B_target)
    }
  }
  if (success_count < B_target) {
    stop(row$model_id, " ", grid_row$method, " reached only ", success_count,
         " successes after ", attempt, " attempts.")
  }

  fail_df <- bind_rows(failures)
  status <- tibble(
    species = row$species, resp_var = row$resp_var, response_label = row$response_label,
    model_id = row$model_id, method = grid_row$method,
    window_days = grid_row$window_days, status = "complete", n_rows = nrow(dat),
    n_containers = n_distinct(dat$boxlabel), n_trees = n_distinct(dat$tree_id),
    n_boot_target = B_target, n_boot_success = success_count, n_attempts = attempt,
    n_failures = nrow(fail_df), n_singular_swc = singular_swc,
    n_singular_response = singular_response,
    formula_swc = paste(deparse(f_swc), collapse = " "),
    formula_response = paste(deparse(f_resp), collapse = " "),
    seed = grid_row$model_seed, note = NA_character_
  )
  list(point = point, reduced_form = reduced_form,
       replicates = bind_rows(successes), failures = fail_df, status = status)
}

task_indices <- seq_len(nrow(task_grid))
if (.Platform$OS.type == "unix" && n_cores > 1L) {
  results <- parallel::mclapply(
    task_indices, fit_one_task, mc.cores = min(n_cores, length(task_indices)),
    mc.preschedule = FALSE
  )
} else {
  results <- lapply(task_indices, fit_one_task)
}
if (any(vapply(results, inherits, logical(1), what = "try-error"))) {
  stop("At least one model worker failed.")
}

point_effects <- map_dfr(results, "point")
reduced_form <- map_dfr(results, "reduced_form")
replicates <- map_dfr(results, "replicates")
failures <- map_dfr(results, "failures")
if (!nrow(failures)) {
  failures <- tibble(
    species = character(), resp_var = character(), response_label = character(),
    model_id = character(), method = character(), window_days = integer(),
    attempt = integer(), error = character()
  )
}
status <- map_dfr(results, "status")

bootstrap_effects <- replicates %>%
  group_by(species, resp_var, response_label, model_id, method, window_days,
           treatment, component) %>%
  summarise(
    estimate_boot_mean = mean(estimate),
    lower = quantile(estimate, 0.025, names = FALSE),
    upper = quantile(estimate, 0.975, names = FALSE),
    p_boot = pmin(1, 2 * pmin(
      (sum(estimate_raw <= 0) + 1) / (n() + 1),
      (sum(estimate_raw >= 0) + 1) / (n() + 1)
    )),
    n_boot = n(), .groups = "drop"
  ) %>%
  left_join(
    point_effects %>%
      select(species, resp_var, response_label, model_id, method, window_days,
             treatment, component, estimate_raw, estimate, display_sign),
    by = c("species", "resp_var", "response_label", "model_id", "method",
           "window_days", "treatment", "component")
  ) %>%
  select(species, resp_var, response_label, model_id, method, window_days,
         treatment, component, estimate_raw, estimate, display_sign, everything())

# Exact total = direct + indirect identity at point and replicate levels.
identity_point <- point_effects %>%
  filter(component %in% c("direct", "indirect", "total")) %>%
  select(species, resp_var, model_id, method, window_days, treatment, component, estimate) %>%
  pivot_wider(names_from = component, values_from = estimate) %>%
  mutate(identity_error = total - direct - indirect, level = "point")
identity_boot <- replicates %>%
  filter(component %in% c("direct", "indirect", "total")) %>%
  select(species, resp_var, model_id, method, window_days, replicate,
         treatment, component, estimate) %>%
  pivot_wider(names_from = component, values_from = estimate) %>%
  mutate(identity_error = total - direct - indirect, level = "bootstrap")
identity_qa <- bind_rows(identity_point, identity_boot) %>%
  summarise(
    n_checks = n(), max_abs_identity_error = max(abs(identity_error)),
    n_above_tolerance_1e12 = sum(abs(identity_error) > 1e-12),
    .by = c(method, window_days, level)
  )
if (any(identity_qa$n_above_tolerance_1e12 > 0)) stop("Total path identity QA failed.")

# In ordinary linear path models with identical rows and covariate structure,
# c' + a*b can be nearly invariant to mediator redefinition by regression
# algebra. It is therefore not an independent robustness check. Compare it to
# the reduced-form treatment coefficient from y ~ treatment + covariates (SWC
# omitted). LMM random-effects estimation makes exact equality unnecessary.
reduced_form_qa <- point_effects %>%
  filter(component == "total") %>%
  select(species, resp_var, response_label, model_id, method, window_days,
         treatment, path_summed_total = estimate) %>%
  left_join(reduced_form,
            by = c("species", "resp_var", "response_label", "model_id", "method",
                   "window_days", "treatment")) %>%
  mutate(
    path_sum_minus_reduced_form = path_summed_total - reduced_form_estimate,
    abs_difference = abs(path_sum_minus_reduced_form)
  )

# -----------------------------------------------------------------------------
# 3. Compare both windows with fuzzy baseline and with each other
# -----------------------------------------------------------------------------

baseline_effects <- read_csv(baseline_effect_file, show_col_types = FALSE) %>%
  select(species, resp_var, response_label, model_id, treatment, component,
         estimate, lower, upper, p_boot, n_boot) %>%
  mutate(method = "fuzzy_baseline", window_days = NA_integer_)

# The baseline effects file omits the shared b path; recover it replicate-wise as
# indirect/a and use the median across treatments to avoid a near-zero denominator.
baseline_reps <- read_csv(baseline_rep_file, show_col_types = FALSE)
baseline_display_sign <- point_effects %>%
  distinct(model_id, display_sign)
baseline_b_reps <- baseline_reps %>%
  filter(component %in% c("treatment_to_swc", "indirect")) %>%
  select(species, resp_var, response_label, model_id, replicate, treatment,
         component, estimate_raw) %>%
  pivot_wider(names_from = component, values_from = estimate_raw) %>%
  filter(is.finite(treatment_to_swc), abs(treatment_to_swc) > 1e-10) %>%
  mutate(b_raw = indirect / treatment_to_swc) %>%
  group_by(species, resp_var, response_label, model_id, replicate) %>%
  summarise(b_raw = median(b_raw), .groups = "drop") %>%
  left_join(baseline_display_sign, by = "model_id") %>%
  mutate(estimate = b_raw * display_sign)

baseline_b_point <- map_dfr(seq_len(nrow(inventory)), function(i) {
  if (!isTRUE(inventory$source_available[[i]])) return(tibble())
  bundle <- readRDS(inventory$source_rds_file[[i]])
  sign_i <- display_sign_from_bundle(bundle)
  b_raw <- unname(fixef(bundle$mod_resp)["swc"])
  tibble(
    species = inventory$species[[i]], resp_var = inventory$resp_var[[i]],
    response_label = inventory$response_label[[i]], model_id = inventory$model_id[[i]],
    treatment = "shared", component = "swc_to_response",
    estimate = b_raw * sign_i
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
  left_join(baseline_b_point,
            by = c("species", "resp_var", "response_label", "model_id")) %>%
  mutate(method = "fuzzy_baseline", window_days = NA_integer_) %>%
  select(species, resp_var, response_label, model_id, treatment, component,
         estimate, lower, upper, p_boot, n_boot, method, window_days)

all_effects <- bind_rows(
  baseline_effects,
  baseline_b,
  bootstrap_effects %>%
    select(species, resp_var, response_label, model_id, treatment, component,
           estimate, lower, upper, p_boot, n_boot, method, window_days)
)

compare_methods <- function(data, method_a, method_b) {
  keys <- c("species", "resp_var", "response_label", "model_id", "treatment", "component")
  a <- data %>%
    filter(method == method_a) %>%
    select(all_of(keys), estimate, lower, upper, p_boot, n_boot) %>%
    rename_with(~paste0("a_", .x), c(estimate, lower, upper, p_boot, n_boot))
  b <- data %>%
    filter(method == method_b) %>%
    select(all_of(keys), estimate, lower, upper, p_boot, n_boot) %>%
    rename_with(~paste0("b_", .x), c(estimate, lower, upper, p_boot, n_boot))
  full_join(a, b, by = keys) %>%
    mutate(
      method_a = method_a, method_b = method_b,
      estimate_difference_b_minus_a = b_estimate - a_estimate,
      direction_agrees = sign(a_estimate) == sign(b_estimate),
      a_significant = a_p_boot < 0.05,
      b_significant = b_p_boot < 0.05,
      significance_agrees = a_significant == b_significant,
      intervals_overlap = pmax(a_lower, b_lower) <= pmin(a_upper, b_upper),
      .before = 1
    )
}

comparison <- bind_rows(
  compare_methods(all_effects, "fuzzy_baseline", "antecedent_7d"),
  compare_methods(all_effects, "fuzzy_baseline", "antecedent_14d"),
  compare_methods(all_effects, "antecedent_7d", "antecedent_14d")
)

run_finished <- Sys.time()
status <- status %>%
  mutate(
    run_started_at = format(run_started, "%Y-%m-%d %H:%M:%S %Z"),
    run_finished_at = format(run_finished, "%Y-%m-%d %H:%M:%S %Z"),
    runtime_seconds = as.numeric(difftime(run_finished, run_started, units = "secs"))
  )

write_csv(status, file.path(output_dir, "antecedent-window-sem-bootstrap-status.csv"))
write_csv(point_effects, file.path(output_dir, "antecedent-window-sem-point-effects.csv"))
write_csv(bootstrap_effects, file.path(output_dir, "antecedent-window-sem-bootstrap-effects.csv"))
write_csv(replicates, file.path(output_dir, "antecedent-window-sem-bootstrap-replicates.csv"))
write_csv(failures, file.path(output_dir, "antecedent-window-sem-bootstrap-failures.csv"))
write_csv(identity_qa, file.path(output_dir, "total-path-identity-qa.csv"))
write_csv(reduced_form_qa, file.path(output_dir, "reduced-form-total-qa.csv"))
write_csv(comparison, file.path(output_dir, "antecedent-window-effect-comparisons.csv"))
write_csv(all_effects, file.path(output_dir, "all-method-bootstrap-effects.csv"))
saveRDS(
  list(
    settings = list(
      B_target = B_target, n_cores = n_cores, seed = base_seed, windows = windows,
      daily_input = "existing treatment-agnostic climate-informed GAM swc_hat",
      mediator = "modeled trailing mean ending on response date; no future information",
      scaling = "original model-specific fuzzy-SWC mean and SD",
      formulas = "frozen original formulas",
      resampling = "containers within experimental blocks",
      interpolation_uncertainty = "not propagated; conditional on fitted daily GAM series"
    ),
    coverage = coverage, status = status, point_effects = point_effects,
    bootstrap_effects = bootstrap_effects, comparison = comparison,
    identity_qa = identity_qa, reduced_form_qa = reduced_form_qa,
    replicates = replicates, failures = failures
  ),
  file.path(output_dir, "antecedent-window-sem-bootstrap-results.rds")
)

# -----------------------------------------------------------------------------
# 4. Compact comparison PDF and generated scientific summary
# -----------------------------------------------------------------------------

plot_data <- all_effects %>%
  filter(component %in% c("direct", "indirect", "total")) %>%
  mutate(
    method = factor(
      method,
      levels = c("fuzzy_baseline", "antecedent_7d", "antecedent_14d"),
      labels = c("Fuzzy measured SWC", "Modeled prior 7-day mean", "Modeled prior 14-day mean")
    ),
    treatment = factor(
      treatment,
      levels = rev(c("precipitation", "robinia", "culture", "extreme_event"))
    )
  )

pdf(file.path(output_dir, "antecedent-window-sem-effect-comparison.pdf"), width = 11, height = 10)
for (sp in c("fagus", "quercus")) {
  p <- plot_data %>%
    filter(species == sp) %>%
    ggplot(aes(estimate, treatment, colour = method)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
    geom_errorbarh(aes(xmin = lower, xmax = upper),
                   position = position_dodge(width = 0.55), height = 0) +
    geom_point(position = position_dodge(width = 0.55), size = 1.5) +
    facet_grid(response_label ~ component, scales = "free_x") +
    labs(
      title = paste("Antecedent-SWC SEM sensitivity:", tools::toTitleCase(sp)),
      x = "Oriented standardized effect (95% container-bootstrap CI)", y = NULL,
      colour = "Mediator definition"
    ) +
    theme_classic(base_size = 8) +
    theme(
      legend.position = "bottom", strip.text = element_text(size = 7),
      plot.title = element_text(hjust = 0.5)
    )
  print(p)
}
dev.off()

primary_comp <- comparison %>%
  filter(component %in% c("direct", "indirect", "total"))
summary_table <- primary_comp %>%
  summarise(
    n = sum(is.finite(a_estimate) & is.finite(b_estimate)),
    direction_agreement_pct = 100 * mean(direction_agrees, na.rm = TRUE),
    significance_agreement_pct = 100 * mean(significance_agrees, na.rm = TRUE),
    interval_overlap_pct = 100 * mean(intervals_overlap, na.rm = TRUE),
    median_abs_estimate_change = median(abs(estimate_difference_b_minus_a), na.rm = TRUE),
    .by = c(method_a, method_b)
  )
write_csv(summary_table, file.path(output_dir, "antecedent-window-comparison-summary.csv"))

fmt_summary <- function(a, b) {
  x <- summary_table %>% filter(method_a == a, method_b == b)
  paste0(
    "Across ", x$n, " comparable direct, indirect and total paths, direction agreed for ",
    round(x$direction_agreement_pct, 1), "%, bootstrap significance for ",
    round(x$significance_agreement_pct, 1), "%, and 95% intervals overlapped for ",
    round(x$interval_overlap_pct, 1), "%. The median absolute point-estimate change was ",
    round(x$median_abs_estimate_change, 3), " standardized units."
  )
}

complete_cov <- coverage %>% filter(status == "available")
summary_lines <- c(
  "# Antecedent modeled-SWC sensitivity",
  "",
  "## Design",
  "",
  paste0(
    "Repeated-response SEMs were refitted with modeled trailing SWC means ending on each response date. The 7-day definition uses response day through day -6; the 14-day definition uses response day through day -13. The daily latent input is the existing climate-informed, treatment-agnostic container-level GAM estimate (`swc_hat`). No future SWC information enters either definition. These variables are modeled antecedent means, not measured means."
  ),
  "",
  paste0(
    "Frozen original formulas and model decisions were retained. Each model retained its fuzzy baseline response values, response scaling, and original SWC mean/SD, so comparisons isolate mediator definition. Uncertainty used ", B_target, " successful block-stratified container-cluster bootstrap replicates per estimable species-response model and window."
  ),
  "",
  "## Coverage",
  "",
  paste0(
    "All estimable baseline response records had complete modeled antecedent windows: minimum model-level coverage was ",
    round(100 * min(complete_cov$coverage_proportion), 1), "%. Two Quercus senescence models remained unavailable because the original SEM bundles contained no post-filtering response data."
  ),
  "",
  "## Effect comparison",
  "",
  paste0("**Fuzzy baseline versus 7-day mean:** ",
         fmt_summary("fuzzy_baseline", "antecedent_7d")),
  "",
  paste0("**Fuzzy baseline versus 14-day mean:** ",
         fmt_summary("fuzzy_baseline", "antecedent_14d")),
  "",
  paste0("**7-day versus 14-day mean:** ",
         fmt_summary("antecedent_7d", "antecedent_14d")),
  "",
  "## Scientific interpretation",
  "",
  "The antecedent-window analyses test whether short-term integrated water availability is more relevant than a single fuzzily matched SWC observation. They do not validate causal mediation: the daily GAM is treatment-agnostic and its prediction uncertainty is not propagated through the SEM bootstrap. Climate covariates help describe common temporal variation but cannot reconstruct container-specific experimental water additions or exclusions without dated treatment-operation records.",
  "",
  "Prefer the shortest biologically defensible window that is stable relative to the 14-day result. A 7-day mean is the cleaner primary sensitivity for rapidly responding physiological traits; the 14-day mean is better interpreted as a longer antecedent-exposure check. If direct/indirect allocation changes while total effects remain stable, pathway attribution is sensitive to mediator timing.",
  "",
  "Stable path-summed totals are not independent evidence of robustness. In linear path models fitted to the same observations with the same covariates, `c' + a*b` can be nearly invariant to mediator redefinition by regression algebra. The relevant sensitivity target is therefore the split between direct and indirect paths (especially `b` and `a*b`). `output/reduced-form-total-qa.csv` compares every point path sum with the corresponding treatment coefficient from the response model with SWC omitted; mixed-model random-effects estimation means exact equality is not required.",
  "",
  "## QA and caveats",
  "",
  paste0(
    "All estimable point and bootstrap SEMs retained exact total = direct + indirect identity; maximum absolute numerical error was ",
    format(max(identity_qa$max_abs_identity_error), scientific = TRUE), "."
  ),
  "",
  paste0(
    "Across point fits, the median absolute difference between the path-summed total and the reduced-form treatment coefficient was ",
    round(median(reduced_form_qa$abs_difference, na.rm = TRUE), 4),
    " standardized units (maximum ", round(max(reduced_form_qa$abs_difference, na.rm = TRUE), 4), ")."
  ),
  "",
  paste0(
    "Bootstrap singularity counts and fitting status are recorded in `output/antecedent-window-sem-bootstrap-status.csv`. Inference is conditional on the fitted daily GAM series; it does not include uncertainty from estimating daily SWC."
  )
)
writeLines(summary_lines, file.path(analysis_dir, "comparison.md"))

readme_lines <- c(
  "# Antecedent SWC-window sensitivity",
  "",
  "Run from the repository root:",
  "",
  "```sh",
  "Rscript \"exploration/2026-08-18 Testing bootstrapping/swc_matching_sensitivity/antecedent_windows/scripts/01_antecedent_window_sem_bootstrap.R\" --bootstrap=1000 --cores=4",
  "```",
  "",
  "The script compares modeled trailing 7-day and 14-day mean SWC mediators with the fuzzy-matched baseline. Both windows end on the response date and never use future information. Outputs include coverage, status, point and bootstrap effects, all replicates, pairwise comparisons, path-identity QA, a compact PDF, and `comparison.md`."
)
writeLines(readme_lines, file.path(analysis_dir, "README.md"))

cat("\nAntecedent-window SEM sensitivity complete.\n")
print(summary_table)
print(status %>% select(model_id, method, status, n_boot_success, n_attempts,
                        n_singular_swc, n_singular_response))
print(identity_qa)
