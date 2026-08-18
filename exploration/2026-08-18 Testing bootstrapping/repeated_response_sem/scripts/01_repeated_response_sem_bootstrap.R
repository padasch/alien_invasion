#!/usr/bin/env Rscript

# Container-cluster bootstrap for the 14 non-phenology species-response SEM
# slots used by Figure 6. Two Quercus senescence slots currently contain no
# post-filtering data and are retained as explicit unavailable status rows.

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
n_cores <- as.integer(arg_value("--cores", as.character(min(4L, parallel::detectCores(logical = FALSE)))))
base_seed <- as.integer(arg_value("--seed", "20260818"))
run_started_at <- Sys.time()
if (!is.finite(B_target) || B_target < 1L) stop("--bootstrap must be a positive integer.")
if (!is.finite(n_cores) || n_cores < 1L) stop("--cores must be a positive integer.")

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path_arg <- gsub("~[+]~", " ", sub("^--file=", "", script_arg[[1]]), fixed = FALSE)
script_file <- normalizePath(script_path_arg, winslash = "/", mustWork = TRUE)
analysis_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/", mustWork = TRUE)
project_root <- normalizePath(file.path(analysis_dir, "..", "..", ".."), winslash = "/", mustWork = TRUE)
output_dir <- file.path(analysis_dir, "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

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

find_current_sem_source <- function(species, resp_var) {
  root <- file.path(project_root, "output")
  pattern <- paste0(
    "sem-tree-.*-", resp_var, "-", species,
    "-soil-both_without_soil_treatment-noInt-scaled-.*-all-swcMeas.*-matrix_data[.]csv$"
  )
  files <- list.files(root, pattern = pattern, recursive = TRUE, full.names = TRUE)
  if (!length(files)) return(list(matrix_file = NA_character_, rds_file = NA_character_))

  date_part <- stringr::str_extract(files, "[0-9]{4}-[0-9]{2}-[0-9]{2}")
  prefer_rfe <- grepl("rfeAIC2", files, fixed = TRUE)
  ord <- order(as.Date(date_part), prefer_rfe, file.info(files)$mtime,
               decreasing = TRUE, na.last = TRUE)
  files <- files[ord]

  readable <- vapply(files, function(file) {
    cols <- tryCatch(names(readr::read_csv(file, show_col_types = FALSE, n_max = 0)),
                     error = function(e) character())
    "path_type" %in% cols
  }, logical(1))

  if (any(readable)) {
    matrix_file <- normalizePath(files[which(readable)[[1]]], winslash = "/", mustWork = TRUE)
    rds_file <- sub("-matrix_data[.]csv$", ".rds", matrix_file)
    return(list(matrix_file = matrix_file, rds_file = rds_file))
  }

  # Preserve the latest failed/empty slot in the audit trail.
  matrix_file <- normalizePath(files[[1]], winslash = "/", mustWork = TRUE)
  rds_file <- sub("-matrix_data[.]csv$", ".rds", matrix_file)
  list(matrix_file = matrix_file, rds_file = rds_file)
}

source_inventory <- pmap_dfr(model_grid, function(species, resp_var, response_label,
                                                   model_id, model_seed) {
  src <- find_current_sem_source(species, resp_var)
  bundle <- if (!is.na(src$rds_file) && file.exists(src$rds_file)) {
    tryCatch(readRDS(src$rds_file), error = function(e) NULL)
  } else NULL
  usable <- !is.null(bundle) && inherits(bundle$mod_swc, "merMod") &&
    inherits(bundle$mod_resp, "merMod") && !is.null(bundle$data) && nrow(bundle$data) > 0
  note <- if (usable) NA_character_ else {
    bundle$note %||% "No readable fitted SEM bundle was found."
  }
  tibble(
    species, resp_var, response_label, model_id, model_seed,
    source_matrix_file = src$matrix_file,
    source_rds_file = src$rds_file,
    source_available = usable,
    source_note = note
  )
})

coef_name_for_factor <- function(model, factor_name) {
  hits <- grep(factor_name, names(fixef(model)), fixed = TRUE, value = TRUE)
  hits <- hits[!grepl(":", hits, fixed = TRUE)]
  if (length(hits) != 1L) return(NA_character_)
  hits[[1]]
}

extract_point_effects <- function(bundle, species, resp_var, response_label, model_id, source_rds_file) {
  eff <- bundle$effects
  if (is.null(eff) || !nrow(eff)) return(tibble())
  display_sign <- 1
  matrix_data <- bundle$matrix_data
  if (!is.null(matrix_data) && nrow(matrix_data) &&
      all(c("estimate", "estimate_raw") %in% names(matrix_data))) {
    ratios <- matrix_data$estimate / matrix_data$estimate_raw
    ratios <- ratios[is.finite(ratios) & abs(matrix_data$estimate_raw) > 1e-12]
    if (length(ratios)) display_sign <- sign(stats::median(ratios))
  }

  bind_rows(
    eff %>% transmute(treatment = factor, component = "treatment_to_swc",
                      estimate_raw = a, se_old = se_a, p_old = p_a),
    eff %>% transmute(treatment = factor, component = "direct",
                      estimate_raw = c_direct, se_old = se_c, p_old = p_c),
    eff %>% transmute(treatment = factor, component = "indirect",
                      estimate_raw = indirect, se_old = se_ind, p_old = p_ind),
    eff %>% transmute(treatment = factor, component = "total",
                      estimate_raw = total, se_old = se_tot, p_old = p_tot)
  ) %>%
    mutate(
      species = species,
      resp_var = resp_var,
      response_label = response_label,
      model_id = model_id,
      estimate = if_else(component == "treatment_to_swc", estimate_raw,
                         estimate_raw * display_sign),
      display_sign = display_sign,
      lower_old_raw = estimate_raw - 1.96 * se_old,
      upper_old_raw = estimate_raw + 1.96 * se_old,
      lower_old = if_else(
        component == "treatment_to_swc", lower_old_raw,
        pmin(lower_old_raw * display_sign, upper_old_raw * display_sign)
      ),
      upper_old = if_else(
        component == "treatment_to_swc", upper_old_raw,
        pmax(lower_old_raw * display_sign, upper_old_raw * display_sign)
      ),
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

  map_dfr(treatments, function(treatment) {
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
        status = "unavailable", n_rows = 0L, n_containers = 0L, n_trees = 0L,
        n_boot_target = B_target, n_boot_success = 0L, n_attempts = 0L,
        n_failures = 0L, n_singular_swc = 0L, n_singular_response = 0L,
        formula_swc = NA_character_, formula_response = NA_character_,
        modeled_treatments = NA_character_, seed = model_seed,
        source_rds_file = source_rds_file,
        note = inventory_row$source_note[[1]]
      )
    ))
  }

  bundle <- readRDS(source_rds_file)
  data <- bundle$data
  formula_swc <- formula(bundle$mod_swc)
  formula_response <- formula(bundle$mod_resp)
  treatments <- as.character(bundle$effects$factor)
  point <- extract_point_effects(bundle, species, resp_var, response_label, model_id, source_rds_file)
  display_sign <- unique(point$display_sign)
  if (length(display_sign) != 1L) stop("Could not identify one response display sign for ", model_id)

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
          estimate = if_else(component == "a", estimate_raw,
                             estimate_raw * display_sign),
          component = recode(component, a = "treatment_to_swc"),
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
    status = "complete", n_rows = nrow(data),
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
  list(replicates = replicate_df, failures = failure_df, point = point, status = status)
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

comparison <- point_effects %>%
  left_join(
    bootstrap_effects %>%
      select(species, resp_var, response_label, model_id, treatment, component,
             estimate_raw_boot_mean, estimate_boot_mean, lower_raw, upper_raw,
             lower, upper, p_boot, n_boot),
    by = c("species", "resp_var", "response_label", "model_id", "treatment", "component")
  ) %>%
  mutate(
    old_ci_excludes_zero = lower_old_raw > 0 | upper_old_raw < 0,
    bootstrap_ci_excludes_zero = lower_raw > 0 | upper_raw < 0,
    old_significant = !is.na(p_old) & p_old < 0.05,
    bootstrap_significant = !is.na(p_boot) & p_boot < 0.05,
    direction_agrees = sign(estimate_raw) == sign(estimate_raw_boot_mean),
    significance_changed = old_significant != bootstrap_significant,
    ci_zero_exclusion_changed = old_ci_excludes_zero != bootstrap_ci_excludes_zero,
    old_interval_width = upper_old_raw - lower_old_raw,
    bootstrap_interval_width = upper_raw - lower_raw,
    interval_width_ratio = bootstrap_interval_width / old_interval_width
  )

figure6_ready <- bootstrap_effects %>%
  filter(component == "total") %>%
  transmute(
    species, treatment, response_label, resp_var,
    estimate, lower, upper, p_boot, n_boot,
    lower_95 = lower, upper_95 = upper, n_boot_success = n_boot,
    estimate_raw, lower_raw, upper_raw,
    source = "block-stratified container-cluster bootstrap"
  )

write_csv(source_inventory, file.path(output_dir, "repeated-response-sem-source-inventory.csv"))
write_csv(status, file.path(output_dir, "repeated-response-sem-bootstrap-status.csv"))
write_csv(point_effects, file.path(output_dir, "repeated-response-sem-point-delta-effects.csv"))
write_csv(bootstrap_effects, file.path(output_dir, "repeated-response-sem-bootstrap-effects.csv"))
write_csv(comparison, file.path(output_dir, "repeated-response-sem-delta-vs-bootstrap-comparison.csv"))
write_csv(figure6_ready, file.path(output_dir, "figure6-ready-repeated-response-sem-path-summed-total.csv"))
write_csv(failures, file.path(output_dir, "repeated-response-sem-bootstrap-failures.csv"))
saveRDS(
  list(
    settings = list(B_target = B_target, cores = n_cores, base_seed = base_seed,
                    resampling = "containers within block", formulas = "frozen original rfeAIC2 formulas",
                    run_started_at = run_started_at, run_finished_at = run_finished_at,
                    runtime_seconds = runtime_seconds),
    inventory = source_inventory, status = status, point_effects = point_effects,
    bootstrap_effects = bootstrap_effects, comparison = comparison,
    replicates = replicates, failures = failures
  ),
  file.path(output_dir, "repeated-response-sem-bootstrap-results.rds")
)

# CSV is intentionally retained for transparent effect-level auditing.
write_csv(replicates, file.path(output_dir, "repeated-response-sem-bootstrap-replicates.csv"))

cat("\nRepeated-response SEM bootstrap complete.\n")
cat(sprintf("Runtime: %.1f seconds.\n", runtime_seconds))
print(status %>% select(model_id, status, n_boot_success, n_attempts, n_failures,
                        n_singular_swc, n_singular_response))
cat("\nSignificance changes (direct/indirect/total):\n")
print(comparison %>% filter(component %in% c("direct", "indirect", "total"), significance_changed) %>%
        select(species, resp_var, treatment, component, estimate_raw, p_old, p_boot,
               lower_raw, upper_raw))
