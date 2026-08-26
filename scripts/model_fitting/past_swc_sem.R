#!/usr/bin/env Rscript

# Sensitivity analysis for the 14 non-phenology species-response SEM slots.
# Response observations retain measured SWC only when it was measured on the
# same day or in the preceding 7 days. Future measurements and the original
# unbounded past fallback are excluded. Existing response/SWC scaling and the
# original frozen formulas are preserved from each baseline model bundle.

suppressPackageStartupMessages({
  library(dplyr)
  library(lme4)
  library(purrr)
  library(readr)
  library(ggplot2)
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
n_cores <- as.integer(arg_value("--cores", as.character(min(4L, parallel::detectCores(logical = FALSE)))))
base_seed <- as.integer(arg_value("--seed", "20260818"))
run_started_at <- Sys.time()
if (!is.finite(B_target) || B_target < 1L) stop("--bootstrap must be a positive integer.")
if (!is.finite(n_cores) || n_cores < 1L) stop("--cores must be a positive integer.")

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path_arg <- gsub("~[+]~", " ", sub("^--file=", "", script_arg[[1]]), fixed = FALSE)
script_file <- normalizePath(script_path_arg, winslash = "/", mustWork = TRUE)
project_root <- normalizePath(file.path(dirname(script_file), "..", ".."), winslash = "/", mustWork = TRUE)
output_dir <- file.path(project_root, "data", "final", "bootstrap", "past_only_7d")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

baseline_dir <- file.path(
  project_root, "data", "final", "bootstrap", "repeated_response_sem"
)
baseline_inventory_file <- file.path(baseline_dir, "repeated-response-sem-source-inventory.csv")
baseline_effects_file <- file.path(baseline_dir, "repeated-response-sem-bootstrap-effects.csv")
if (!file.exists(baseline_inventory_file) || !file.exists(baseline_effects_file)) {
  stop("Baseline repeated-response SEM bootstrap outputs are required.")
}

swc_measured <- read_csv(
  file.path(project_root, "data", "interim", "box_soilwater.csv"),
  show_col_types = FALSE
) %>%
  transmute(boxlabel = as.character(boxlabel), swc_date = as.Date(date), swc_measured = swc)

`%||%` <- function(x, y) if (is.null(x)) y else x

response_specs <- tribble(
  ~resp_var, ~response_label,
  "volume", "Volume (total)",
  "volume_inc_phase_rel", "Volume (incr.)",
  "chl", "Chlorophyll",
  "condition", "Vitality",
  "qy", "Quantum yield",
  "remaining_green", "Senescence (%)",
  "chlavg", "Senescence (Chl)"
)

model_grid <- tidyr::crossing(
  species = c("fagus", "quercus"),
  response_specs
) %>%
  mutate(
    model_id = paste(species, resp_var, sep = "__"),
    model_seed = base_seed + row_number() * 100003L
  )

baseline_inventory <- read_csv(baseline_inventory_file, show_col_types = FALSE)
source_inventory <- model_grid %>%
  left_join(
    baseline_inventory %>%
      select(species, resp_var, source_matrix_file, source_rds_file,
             source_available, source_note),
    by = c("species", "resp_var")
  ) %>%
  mutate(
    source_available = replace_na(source_available, FALSE),
    source_note = if_else(
      !source_available & is.na(source_note),
      "No baseline fitted SEM bundle was available.", source_note
    )
  )

match_past_swc_7d <- function(data) {
  keys <- data %>%
    transmute(boxlabel = as.character(boxlabel), response_date = as.Date(date)) %>%
    distinct()
  matches <- pmap_dfr(keys, function(boxlabel, response_date) {
    candidates <- swc_measured %>%
      filter(
        .data$boxlabel == .env$boxlabel,
        swc_date <= response_date,
        swc_date >= response_date - 7
      ) %>%
      arrange(desc(swc_date))
    if (!nrow(candidates)) {
      return(tibble(
        boxlabel = boxlabel, response_date = response_date,
        swc_date = as.Date(NA), swc_measured = NA_real_, swc_lag_days = NA_integer_
      ))
    }
    tibble(
      boxlabel = boxlabel, response_date = response_date,
      swc_date = candidates$swc_date[[1]],
      swc_measured = candidates$swc_measured[[1]],
      swc_lag_days = as.integer(response_date - candidates$swc_date[[1]])
    )
  })
  matched <- data %>%
    mutate(boxlabel_chr = as.character(boxlabel), response_date = as.Date(date)) %>%
    left_join(matches, by = c("boxlabel_chr" = "boxlabel", "response_date"))
  retained <- matched %>% filter(!is.na(swc_lag_days))
  mismatch <- retained %>%
    filter(!is.na(swc_org), abs(swc_org - swc_measured) > 1e-8)
  if (nrow(mismatch)) {
    stop("Past-only match disagreed with the SWC value retained in the baseline bundle.")
  }
  retained %>%
    select(-boxlabel_chr, -response_date, -swc_measured)
}

make_balance_diagnostics <- function(before, after, species, resp_var, response_label, model_id) {
  dimensions <- c("precipitation", "robinia", "culture", "extreme_event", "date")
  map_dfr(dimensions, function(dimension) {
    before_counts <- before %>%
      mutate(.level = as.character(.data[[dimension]])) %>%
      count(.level, name = "n_baseline")
    after_counts <- after %>%
      mutate(.level = as.character(.data[[dimension]])) %>%
      count(.level, name = "n_past_only")
    full_join(before_counts, after_counts, by = ".level") %>%
      mutate(
        n_baseline = replace_na(n_baseline, 0L),
        n_past_only = replace_na(n_past_only, 0L),
        n_lost = n_baseline - n_past_only,
        pct_lost = if_else(n_baseline > 0, 100 * n_lost / n_baseline, NA_real_),
        dimension = dimension,
        level = .level,
        species = species, resp_var = resp_var, response_label = response_label,
        model_id = model_id,
        .before = 1
      ) %>%
      select(-.level)
  })
}

coef_name_for_factor <- function(model, factor_name) {
  hits <- grep(factor_name, names(fixef(model)), fixed = TRUE, value = TRUE)
  hits <- hits[!grepl(":", hits, fixed = TRUE)]
  if (length(hits) != 1L) return(NA_character_)
  hits[[1]]
}

extract_point_effects <- function(mod_swc, mod_resp, treatments, display_sign,
                                  species, resp_var, response_label, model_id,
                                  source_rds_file) {
  b <- unname(fixef(mod_resp)["swc"])
  if (!length(b) || !is.finite(b)) stop("Missing SWC -> response coefficient.")
  eff <- map_dfr(treatments, function(treatment) {
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
  eff <- bind_rows(
    tibble(treatment = "SWC", component = "swc_to_response", estimate_raw = b),
    eff
  )
  eff %>%
    mutate(
      species = species,
      resp_var = resp_var,
      response_label = response_label,
      model_id = model_id,
      estimate = if_else(component == "treatment_to_swc", estimate_raw,
                         estimate_raw * display_sign),
      display_sign = display_sign,
      source_rds_file = source_rds_file,
      .before = 1
    )
}

resample_containers_within_block <- function(data, replicate_id) {
  dat <- data %>%
    mutate(
      .box_original = as.character(boxlabel),
      .tree_original = as.character(tree_id),
      .block = sub("-.*$", "", .box_original)
    )
  boxes <- dat %>% distinct(.block, .box_original) %>% arrange(.block, .box_original)
  selected <- boxes %>%
    group_by(.block) %>%
    group_modify(~tibble(
      .box_original = sample(.x$.box_original, size = nrow(.x), replace = TRUE),
      .draw = seq_len(nrow(.x))
    )) %>%
    ungroup()

  pieces <- pmap(selected, function(.block, .box_original, .draw) {
    new_box <- paste0("boot", replicate_id, "_", .block, "_", sprintf("%02d", .draw))
    dat %>%
      filter(.data$.box_original == .env$.box_original) %>%
      mutate(
        boxlabel = new_box,
        tree_id = paste(new_box, .tree_original, sep = "__")
      )
  })
  bind_rows(pieces) %>%
    select(-.box_original, -.tree_original, -.block) %>%
    mutate(
      boxlabel = factor(boxlabel),
      tree_id = factor(tree_id)
    )
}

extract_boot_effects <- function(mod_swc, mod_resp, treatments) {
  b <- unname(fixef(mod_resp)["swc"])
  if (!length(b) || !is.finite(b)) stop("Missing or non-finite SWC -> response coefficient.")

  treatment_effects <- map_dfr(treatments, function(treatment) {
    a_name <- coef_name_for_factor(mod_swc, treatment)
    c_name <- coef_name_for_factor(mod_resp, treatment)
    if (is.na(a_name) || is.na(c_name)) {
      stop("Missing treatment coefficient for ", treatment, ".")
    }
    a <- unname(fixef(mod_swc)[a_name])
    c_direct <- unname(fixef(mod_resp)[c_name])
    values <- c(a = a, direct = c_direct, indirect = a * b, total = c_direct + a * b)
    if (any(!is.finite(values))) stop("Non-finite path estimate for ", treatment, ".")
    tibble(treatment, component = names(values), estimate_raw = unname(values))
  })
  bind_rows(
    tibble(treatment = "SWC", component = "swc_to_response", estimate_raw = b),
    treatment_effects
  )
}

bootstrap_one_model <- function(inventory_row) {
  species <- inventory_row$species[[1]]
  resp_var <- inventory_row$resp_var[[1]]
  response_label <- inventory_row$response_label[[1]]
  model_id <- inventory_row$model_id[[1]]
  model_seed <- inventory_row$model_seed[[1]]
  source_rds_file <- inventory_row$source_rds_file[[1]]

  if (!isTRUE(inventory_row$source_available[[1]])) {
    return(list(
      replicates = tibble(), failures = tibble(), point = tibble(),
      status = tibble(
        species, resp_var, response_label, model_id,
        status = "unavailable", n_rows_baseline = 0L, n_rows = 0L,
        sample_loss_n = 0L, sample_loss_pct = NA_real_,
        n_containers = 0L, n_trees = 0L,
        n_boot_target = B_target, n_boot_success = 0L, n_attempts = 0L,
        n_failures = 0L, n_singular_swc = 0L, n_singular_response = 0L,
        formula_swc = NA_character_, formula_response = NA_character_,
        modeled_treatments = NA_character_, seed = model_seed,
        source_rds_file = source_rds_file,
        note = inventory_row$source_note[[1]]
      ),
      lag_data = tibble(), balance = tibble()
    ))
  }

  bundle <- readRDS(source_rds_file)
  data_baseline <- bundle$data
  data <- match_past_swc_7d(data_baseline)
  formula_swc <- formula(bundle$mod_swc)
  formula_response <- formula(bundle$mod_resp)
  treatments <- as.character(bundle$effects$factor)
  display_sign <- 1
  matrix_data <- bundle$matrix_data
  if (!is.null(matrix_data) && nrow(matrix_data) &&
      all(c("estimate", "estimate_raw") %in% names(matrix_data))) {
    ratios <- matrix_data$estimate / matrix_data$estimate_raw
    ratios <- ratios[is.finite(ratios) & abs(matrix_data$estimate_raw) > 1e-12]
    if (length(ratios)) display_sign <- sign(stats::median(ratios))
  }
  if (!nrow(data) || n_distinct(data$boxlabel) < 2L || n_distinct(data$tree_id) < 2L) {
    return(list(
      replicates = tibble(), failures = tibble(), point = tibble(),
      status = tibble(
        species, resp_var, response_label, model_id,
        status = "unavailable_after_matching", n_rows_baseline = nrow(data_baseline),
        n_rows = nrow(data), sample_loss_n = nrow(data_baseline) - nrow(data),
        sample_loss_pct = 100 * (nrow(data_baseline) - nrow(data)) / nrow(data_baseline),
        n_containers = n_distinct(data$boxlabel), n_trees = n_distinct(data$tree_id),
        n_boot_target = B_target, n_boot_success = 0L, n_attempts = 0L,
        n_failures = 0L, n_singular_swc = 0L, n_singular_response = 0L,
        formula_swc = paste(deparse(formula_swc), collapse = " "),
        formula_response = paste(deparse(formula_response), collapse = " "),
        modeled_treatments = paste(treatments, collapse = ";"), seed = model_seed,
        source_rds_file = source_rds_file,
        note = "Too few observations/clusters after past-only 7-day SWC matching."
      ),
      lag_data = data %>%
        transmute(species, resp_var, response_label, model_id,
                  boxlabel = as.character(boxlabel), response_date = as.Date(date),
                  swc_date, swc_lag_days),
      balance = make_balance_diagnostics(
        data_baseline, data, species, resp_var, response_label, model_id
      )
    ))
  }
  point_mod_swc <- suppressMessages(suppressWarnings(
    lmer(formula_swc, data = data, REML = isREML(bundle$mod_swc))
  ))
  point_mod_resp <- suppressMessages(suppressWarnings(
    lmer(formula_response, data = data, REML = isREML(bundle$mod_resp))
  ))
  point <- extract_point_effects(
    point_mod_swc, point_mod_resp, treatments, display_sign,
    species, resp_var, response_label, model_id, source_rds_file
  )

  set.seed(model_seed)
  successes <- vector("list", B_target)
  failures <- list()
  success_count <- 0L
  attempt <- 0L
  singular_swc <- 0L
  singular_response <- 0L
  max_attempts <- max(B_target * 3L, B_target + 250L)

  message("Starting ", model_id, " (", nrow(data), " rows; ",
          n_distinct(data$boxlabel), " containers).")
  while (success_count < B_target && attempt < max_attempts) {
    attempt <- attempt + 1L
    fit_result <- tryCatch({
      boot_data <- resample_containers_within_block(data, attempt)
      mod_swc <- suppressMessages(suppressWarnings(
        lmer(formula_swc, data = boot_data, REML = isREML(bundle$mod_swc))
      ))
      mod_response <- suppressMessages(suppressWarnings(
        lmer(formula_response, data = boot_data, REML = isREML(bundle$mod_resp))
      ))
      estimates <- extract_boot_effects(mod_swc, mod_response, treatments)
      list(
        effects = estimates,
        singular_swc = isSingular(mod_swc, tol = 1e-4),
        singular_response = isSingular(mod_response, tol = 1e-4)
      )
    }, error = function(e) e)

    if (inherits(fit_result, "error")) {
      failures[[length(failures) + 1L]] <- tibble(
        species, resp_var, response_label, model_id, attempt,
        error = conditionMessage(fit_result)
      )
    } else {
      success_count <- success_count + 1L
      singular_swc <- singular_swc + as.integer(fit_result$singular_swc)
      singular_response <- singular_response + as.integer(fit_result$singular_response)
      successes[[success_count]] <- fit_result$effects %>%
        mutate(
          species, resp_var, response_label, model_id,
          replicate = success_count,
          attempt = attempt,
          component = recode(component, a = "treatment_to_swc"),
          estimate = if_else(component == "treatment_to_swc", estimate_raw,
                             estimate_raw * display_sign),
          singular_swc = fit_result$singular_swc,
          singular_response = fit_result$singular_response,
          .before = 1
        )
    }
    if (attempt %% 100L == 0L) {
      message(model_id, ": ", success_count, "/", B_target,
              " successful after ", attempt, " attempts.")
    }
  }

  if (success_count < B_target) {
    stop(model_id, " reached only ", success_count, " successful replicates after ",
         attempt, " attempts.")
  }

  replicate_df <- bind_rows(successes)
  failure_df <- bind_rows(failures)
  status <- tibble(
    species, resp_var, response_label, model_id,
    status = "complete", n_rows_baseline = nrow(data_baseline), n_rows = nrow(data),
    sample_loss_n = nrow(data_baseline) - nrow(data),
    sample_loss_pct = 100 * (nrow(data_baseline) - nrow(data)) / nrow(data_baseline),
    n_containers = n_distinct(data$boxlabel), n_trees = n_distinct(data$tree_id),
    n_boot_target = B_target, n_boot_success = success_count, n_attempts = attempt,
    n_failures = nrow(failure_df), n_singular_swc = singular_swc,
    n_singular_response = singular_response,
    formula_swc = paste(deparse(formula_swc), collapse = " "),
    formula_response = paste(deparse(formula_response), collapse = " "),
    modeled_treatments = paste(treatments, collapse = ";"), seed = model_seed,
    source_rds_file = source_rds_file, note = NA_character_
  )
  message("Completed ", model_id, ": ", success_count, " successful replicates.")
  lag_data <- data %>%
    transmute(
      species, resp_var, response_label, model_id,
      boxlabel = as.character(boxlabel), response_date = as.Date(date),
      swc_date, swc_lag_days
    )
  balance <- make_balance_diagnostics(
    data_baseline, data, species, resp_var, response_label, model_id
  )
  list(
    replicates = replicate_df, failures = failure_df, point = point,
    status = status, lag_data = lag_data, balance = balance
  )
}

inventory_rows <- split(source_inventory, seq_len(nrow(source_inventory)))
if (.Platform$OS.type == "unix" && n_cores > 1L) {
  results <- parallel::mclapply(inventory_rows, bootstrap_one_model,
                                mc.cores = min(n_cores, length(inventory_rows)),
                                mc.preschedule = FALSE)
} else {
  results <- lapply(inventory_rows, bootstrap_one_model)
}

bad_results <- which(vapply(results, inherits, logical(1), what = "try-error"))
if (length(bad_results)) {
  stop("Model bootstrap worker failure(s): ", paste(bad_results, collapse = ", "))
}

replicates <- map_dfr(results, "replicates")
failures <- map_dfr(results, "failures")
if (!nrow(failures)) {
  failures <- tibble(
    species = character(), resp_var = character(), response_label = character(),
    model_id = character(), attempt = integer(), error = character()
  )
}
point_effects <- map_dfr(results, "point")
status <- map_dfr(results, "status")
lag_data <- map_dfr(results, "lag_data")
balance <- map_dfr(results, "balance")
run_finished_at <- Sys.time()
runtime_seconds <- as.numeric(difftime(run_finished_at, run_started_at, units = "secs"))
status <- status %>%
  mutate(
    run_started_at = format(run_started_at, "%Y-%m-%d %H:%M:%S %Z"),
    run_finished_at = format(run_finished_at, "%Y-%m-%d %H:%M:%S %Z"),
    runtime_seconds = runtime_seconds
  )

bootstrap_effects <- replicates %>%
  group_by(species, resp_var, response_label, model_id, treatment, component) %>%
  summarise(
    estimate_raw_boot_mean = mean(estimate_raw),
    estimate_boot_mean = mean(estimate),
    lower_raw = quantile(estimate_raw, 0.025, names = FALSE),
    upper_raw = quantile(estimate_raw, 0.975, names = FALSE),
    lower = quantile(estimate, 0.025, names = FALSE),
    upper = quantile(estimate, 0.975, names = FALSE),
    p_boot = pmin(
      1,
      2 * pmin(
        (sum(estimate_raw <= 0) + 1) / (n() + 1),
        (sum(estimate_raw >= 0) + 1) / (n() + 1)
      )
    ),
    n_boot = n(),
    .groups = "drop"
  ) %>%
  left_join(
    point_effects %>%
      select(species, resp_var, response_label, model_id, treatment, component,
             estimate_raw, estimate, display_sign),
    by = c("species", "resp_var", "response_label", "model_id", "treatment", "component")
  ) %>%
  select(species, resp_var, response_label, model_id, treatment, component,
         estimate_raw, estimate, display_sign, everything())

baseline_effects <- read_csv(baseline_effects_file, show_col_types = FALSE)
comparison <- baseline_effects %>%
  select(species, resp_var, response_label, model_id, treatment, component,
         estimate_baseline = estimate, lower_baseline = lower,
         upper_baseline = upper, p_boot_baseline = p_boot) %>%
  full_join(
    bootstrap_effects %>%
      select(species, resp_var, response_label, model_id, treatment, component,
             estimate_past_only = estimate, lower_past_only = lower,
             upper_past_only = upper, p_boot_past_only = p_boot),
    by = c("species", "resp_var", "response_label", "model_id", "treatment", "component")
  ) %>%
  mutate(
    estimate_change = estimate_past_only - estimate_baseline,
    direction_agrees = sign(estimate_baseline) == sign(estimate_past_only),
    significant_baseline = !is.na(p_boot_baseline) & p_boot_baseline < 0.05,
    significant_past_only = !is.na(p_boot_past_only) & p_boot_past_only < 0.05,
    significance_changed = significant_baseline != significant_past_only,
    ci_overlap = pmax(lower_baseline, lower_past_only) <=
      pmin(upper_baseline, upper_past_only)
  )

figure6_ready <- bootstrap_effects %>%
  filter(component == "total") %>%
  transmute(
    species, treatment, response_label, resp_var,
    estimate, lower, upper, p_boot, n_boot,
    lower_95 = lower, upper_95 = upper, n_boot_success = n_boot,
    estimate_raw, lower_raw, upper_raw,
    source = "past-only 7-day SWC; block-stratified container-cluster bootstrap"
  )

lag_diagnostics <- lag_data %>%
  group_by(species, resp_var, response_label, model_id, swc_lag_days) %>%
  summarise(
    n_rows = n(),
    n_container_dates = n_distinct(paste(boxlabel, response_date)),
    .groups = "drop"
  )

qa_total_identity <- replicates %>%
  filter(component %in% c("direct", "indirect", "total")) %>%
  select(species, resp_var, model_id, replicate, treatment, component, estimate_raw) %>%
  pivot_wider(names_from = component, values_from = estimate_raw) %>%
  mutate(
    total_minus_direct_indirect = total - direct - indirect,
    qa_pass = is.finite(total_minus_direct_indirect) &
      abs(total_minus_direct_indirect) < 1e-10
  )

qa_summary <- qa_total_identity %>%
  group_by(species, resp_var, model_id) %>%
  summarise(
    n_checks = n(),
    n_failures = sum(!qa_pass),
    max_abs_difference = max(abs(total_minus_direct_indirect), na.rm = TRUE),
    .groups = "drop"
  )

plot_data <- comparison %>%
  filter(component %in% c("direct", "indirect", "total"), !is.na(estimate_change)) %>%
  mutate(
    species = recode(species, fagus = "Fagus", quercus = "Quercus"),
    component = factor(component, levels = c("direct", "indirect", "total"),
                       labels = c("Direct", "Indirect via SWC", "Total")),
    treatment = recode(treatment, precipitation = "Drought", robinia = "With Robinia",
                       culture = "Mixed culture", extreme_event = "Extreme event"),
    response_label = factor(
      response_label,
      levels = rev(c("Volume (total)", "Volume (incr.)", "Chlorophyll", "Vitality",
                     "Quantum yield", "Senescence (%)", "Senescence (Chl)"))
    )
  )

if (nrow(plot_data)) {
  p <- ggplot(plot_data, aes(treatment, response_label, fill = estimate_change)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%+.2f", estimate_change)), size = 2.5) +
    facet_grid(species ~ component, scales = "free_y", space = "free_y") +
    scale_fill_gradient2(
      low = "#b2182b", mid = "white", high = "#2166ac", midpoint = 0,
      name = expression(Delta*" effect")
    ) +
    labs(
      x = NULL, y = NULL,
      caption = "Cell = past-only 7-day estimate minus current fuzzy-match estimate."
    ) +
    theme_classic(base_size = 9) +
    theme(
      panel.grid = element_blank(), strip.background = element_rect(fill = "grey92"),
      axis.text.x = element_text(angle = 35, hjust = 1),
      legend.position = "right"
    )
  ggsave(file.path(output_dir, "past-only-vs-fuzzy-effect-deltas.pdf"), p,
         width = 12, height = 7.5, units = "in")
}

write_csv(source_inventory, file.path(output_dir, "past-only-sem-source-inventory.csv"))
write_csv(status, file.path(output_dir, "past-only-sem-bootstrap-status.csv"))
write_csv(point_effects, file.path(output_dir, "past-only-sem-point-effects.csv"))
write_csv(bootstrap_effects, file.path(output_dir, "past-only-sem-bootstrap-effects.csv"))
write_csv(comparison, file.path(output_dir, "past-only-vs-fuzzy-bootstrap-comparison.csv"))
write_csv(figure6_ready, file.path(output_dir, "past-only-figure6-ready-total-effects.csv"))
write_csv(failures, file.path(output_dir, "past-only-sem-bootstrap-failures.csv"))
write_csv(lag_diagnostics, file.path(output_dir, "past-only-lag-diagnostics.csv"))
write_csv(balance, file.path(output_dir, "past-only-sample-balance-diagnostics.csv"))
write_csv(qa_summary, file.path(output_dir, "past-only-qa-total-identity.csv"))
saveRDS(
  list(
    settings = list(B_target = B_target, cores = n_cores, base_seed = base_seed,
                    swc_matching = "same day or most recent measured SWC in preceding 7 days",
                    scaling = "preserved baseline y and SWC scaling constants",
                    resampling = "containers within block", formulas = "frozen original rfeAIC2 formulas",
                    run_started_at = run_started_at, run_finished_at = run_finished_at,
                    runtime_seconds = runtime_seconds),
    inventory = source_inventory, status = status, point_effects = point_effects,
    bootstrap_effects = bootstrap_effects, comparison = comparison,
    lag_diagnostics = lag_diagnostics, balance = balance,
    qa_total_identity = qa_total_identity, qa_summary = qa_summary,
    replicates = replicates, failures = failures
  ),
  file.path(output_dir, "past-only-sem-bootstrap-results.rds")
)

# CSV is intentionally retained for transparent effect-level auditing.
write_csv(replicates, file.path(output_dir, "past-only-sem-bootstrap-replicates.csv"))

cat("\nPast-only 7-day repeated-response SEM bootstrap complete.\n")
cat(sprintf("Runtime: %.1f seconds.\n", runtime_seconds))
print(status %>% select(model_id, status, n_boot_success, n_attempts, n_failures,
                        n_singular_swc, n_singular_response))
cat("\nSignificance changes versus fuzzy baseline (direct/indirect/total):\n")
print(comparison %>% filter(component %in% c("direct", "indirect", "total"), significance_changed) %>%
        select(species, resp_var, treatment, component, estimate_baseline,
               estimate_past_only, p_boot_baseline, p_boot_past_only))
cat("\nQA total = direct + indirect failures: ", sum(!qa_total_identity$qa_pass), "\n", sep = "")
