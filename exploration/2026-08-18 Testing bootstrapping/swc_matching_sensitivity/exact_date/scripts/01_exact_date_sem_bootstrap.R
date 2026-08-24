#!/usr/bin/env Rscript

# Exact-date measured-SWC sensitivity analysis for repeated-response SEMs.
# Reuses the cached baseline data, scaling, factor coding, and selected formulas;
# only rows sharing an exact response/SWC container-date are retained.

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(lme4); library(purrr)
  library(readr); library(tibble); library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default) {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (!length(hit)) default else sub(paste0("^", flag, "="), "", hit[[1]])
}
B <- as.integer(arg_value("--bootstrap", "1000"))
cores <- as.integer(arg_value("--cores", "4"))
base_seed <- as.integer(arg_value("--seed", "2026081811"))
stopifnot(B > 0L, cores > 0L)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path_arg <- gsub("~[+]~", " ", sub("^--file=", "", script_arg[[1]]))
script_file <- normalizePath(script_path_arg, winslash = "/", mustWork = TRUE)
analysis_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
find_root <- function(x) {
  repeat {
    if (file.exists(file.path(x, "functions", "11-bootstrap-inference.R"))) return(x)
    parent <- dirname(x); if (identical(parent, x)) stop("Project root not found."); x <- parent
  }
}
project_root <- find_root(analysis_dir)
out <- file.path(analysis_dir, "output"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
source(file.path(project_root, "functions", "11-bootstrap-inference.R"), local = TRUE)
started <- Sys.time()

baseline_out <- file.path(project_root, "exploration", "2026-08-18 Testing bootstrapping",
                          "repeated_response_sem", "output")
inventory <- read_csv(file.path(baseline_out, "repeated-response-sem-source-inventory.csv"), show_col_types = FALSE) %>%
  mutate(model_seed = base_seed + row_number() * 100003L)
baseline_effects <- read_csv(file.path(baseline_out, "repeated-response-sem-bootstrap-effects.csv"), show_col_types = FALSE)
baseline_reps <- read_csv(file.path(baseline_out, "repeated-response-sem-bootstrap-replicates.csv"), show_col_types = FALSE)
swc_measured <- read_csv(file.path(project_root, "data", "interim", "box_soilwater.csv"), show_col_types = FALSE) %>%
  transmute(boxlabel = as.character(boxlabel), date = as.Date(date), swc_exact_org = swc)
if (anyDuplicated(swc_measured[c("boxlabel", "date")])) stop("Duplicate measured SWC container-date keys.")

display_sign <- function(bundle) {
  x <- bundle$matrix_data
  if (is.null(x) || !all(c("estimate", "estimate_raw") %in% names(x))) return(1)
  r <- x$estimate / x$estimate_raw
  r <- r[is.finite(r) & abs(x$estimate_raw) > 1e-12]
  if (length(r)) sign(median(r)) else 1
}
coef_for <- function(model, treatment) {
  x <- grep(treatment, names(fixef(model)), fixed = TRUE, value = TRUE)
  x <- x[!grepl(":", x, fixed = TRUE)]
  if (length(x) == 1L) x else NA_character_
}
extract_paths <- function(mod_swc, mod_resp, treatments) {
  b <- unname(fixef(mod_resp)["swc"])
  if (!is.finite(b)) stop("SWC-response coefficient unavailable.")
  tr <- map_dfr(treatments, function(treatment) {
    an <- coef_for(mod_swc, treatment); cn <- coef_for(mod_resp, treatment)
    if (is.na(an) || is.na(cn)) stop("Treatment coefficient unavailable: ", treatment)
    a <- unname(fixef(mod_swc)[an]); direct <- unname(fixef(mod_resp)[cn])
    tibble(treatment, component = c("treatment_to_swc", "direct", "indirect", "total"),
           estimate_raw = c(a, direct, a * b, direct + a * b))
  })
  bind_rows(tibble(treatment = "SWC", component = "swc_to_response", estimate_raw = b), tr)
}
scale_constants <- function(z, original, variable, model_id) {
  fit <- lm(z ~ original)
  slope <- unname(coef(fit)[[2]]); intercept <- unname(coef(fit)[[1]])
  center <- -intercept / slope; scale <- 1 / slope
  tibble(model_id, variable, center, scale,
         max_reconstruction_error = max(abs(z - (original - center) / scale)))
}
balance <- function(data, model_id, sample) {
  map_dfr(c("precipitation", "robinia", "culture", "extreme_event"), function(variable) {
    data %>% count(level = as.character(.data[[variable]]), name = "n") %>%
      mutate(model_id, sample, variable, proportion = n / sum(n), .before = 1)
  })
}
resample_boxes <- function(data, rep_id) {
  dat <- data %>% mutate(.old_box = as.character(boxlabel), .old_tree = as.character(tree_id),
                         .block = sub("-.*$", "", .old_box))
  selected <- dat %>% distinct(.block, .old_box) %>% group_by(.block) %>%
    group_modify(~tibble(.old_box = sample(.x$.old_box, nrow(.x), replace = TRUE), draw = seq_len(nrow(.x)))) %>%
    ungroup()
  pmap_dfr(selected, function(.block, .old_box, draw) {
    new_box <- paste0("boot", rep_id, "_", .block, "_", sprintf("%02d", draw))
    dat %>% filter(.data$.old_box == .env$.old_box) %>%
      mutate(boxlabel = new_box, tree_id = paste(new_box, .old_tree, sep = "__"))
  }) %>% select(-.old_box, -.old_tree, -.block) %>%
    mutate(boxlabel = factor(boxlabel), tree_id = factor(tree_id))
}

empty_result <- function(row, status, note, sample, bal = tibble(), scaling = tibble(),
                         fs = NA_character_, fr = NA_character_) {
  list(reps = tibble(), point = tibble(), failures = tibble(), sample = sample, balance = bal,
       scaling = scaling,
       status = tibble(species = row$species, resp_var = row$resp_var,
         response_label = row$response_label, model_id = row$model_id, status,
         n_rows_baseline = sample$n_rows_baseline, n_rows_exact = sample$n_rows_exact,
         retained_fraction = sample$retained_fraction, n_containers_exact = sample$n_containers_exact,
         n_trees_exact = sample$n_trees_exact, n_boot_target = B, n_boot_success = 0L,
         n_attempts = 0L, n_failures = 0L, n_singular_swc = 0L, n_singular_response = 0L,
         formula_swc = fs, formula_response = fr, seed = row$model_seed,
         source_rds_file = row$source_rds_file, note))
}

run_model <- function(row) {
  row <- as.list(row)
  blank_sample <- tibble(model_id = row$model_id, n_rows_baseline = 0L, n_rows_exact = 0L,
    retained_fraction = NA_real_, n_containers_baseline = 0L, n_containers_exact = 0L,
    n_trees_baseline = 0L, n_trees_exact = 0L, n_dates_baseline = 0L, n_dates_exact = 0L,
    max_abs_swc_exact_difference = NA_real_)
  if (!isTRUE(row$source_available) || !file.exists(row$source_rds_file))
    return(empty_result(row, "baseline_unavailable", row$source_note, blank_sample))
  bundle <- readRDS(row$source_rds_file); d0 <- bundle$data
  d <- d0 %>% mutate(boxlabel = as.character(boxlabel), date = as.Date(date)) %>%
    inner_join(swc_measured, by = c("boxlabel", "date")) %>%
    mutate(swc_exact_difference = swc_org - swc_exact_org,
           boxlabel = factor(boxlabel), tree_id = droplevels(factor(tree_id)))
  sample <- tibble(model_id = row$model_id, n_rows_baseline = nrow(d0), n_rows_exact = nrow(d),
    retained_fraction = nrow(d) / nrow(d0), n_containers_baseline = n_distinct(d0$boxlabel),
    n_containers_exact = n_distinct(d$boxlabel), n_trees_baseline = n_distinct(d0$tree_id),
    n_trees_exact = n_distinct(d$tree_id), n_dates_baseline = n_distinct(d0$date),
    n_dates_exact = n_distinct(d$date),
    max_abs_swc_exact_difference = if (nrow(d)) max(abs(d$swc_exact_difference)) else NA_real_)
  bal <- bind_rows(balance(d0, row$model_id, "fuzzy_baseline"), balance(d, row$model_id, "exact_date"))
  scaling <- bind_rows(scale_constants(d0$swc, d0$swc_org, "swc", row$model_id),
                       scale_constants(d0$y, d0$y_org, "response", row$model_id))
  fs <- formula(bundle$mod_swc); fr <- formula(bundle$mod_resp)
  fs_txt <- paste(deparse(fs), collapse = " "); fr_txt <- paste(deparse(fr), collapse = " ")
  if (!nrow(d)) return(empty_result(row, "no_exact_date_rows",
    "No response observation shared an exact container-date with measured SWC.", sample, bal, scaling, fs_txt, fr_txt))
  if (sample$max_abs_swc_exact_difference > 1e-10) return(empty_result(row, "swc_value_mismatch",
    "Raw exact-date SWC differs from cached swc_org.", sample, bal, scaling, fs_txt, fr_txt))
  treatments <- as.character(bundle$effects$factor); orient <- display_sign(bundle)
  point_fit <- tryCatch({
    ms <- suppressMessages(suppressWarnings(lmer(fs, d, REML = isREML(bundle$mod_swc))))
    mr <- suppressMessages(suppressWarnings(lmer(fr, d, REML = isREML(bundle$mod_resp))))
    list(ms = ms, mr = mr, paths = extract_paths(ms, mr, treatments))
  }, error = function(e) e)
  if (inherits(point_fit, "error")) return(empty_result(row, "formula_not_estimable",
    conditionMessage(point_fit), sample, bal, scaling, fs_txt, fr_txt))
  point <- point_fit$paths %>% mutate(species = row$species, resp_var = row$resp_var,
    response_label = row$response_label, model_id = row$model_id, display_sign = orient,
    estimate = if_else(component == "treatment_to_swc", estimate_raw, estimate_raw * orient), .before = 1)

  set.seed(row$model_seed); successes <- vector("list", B); failures <- list()
  success <- 0L; attempt <- 0L; sing_s <- 0L; sing_r <- 0L; max_attempts <- max(3L * B, B + 250L)
  message("Starting ", row$model_id, " exact-date (", nrow(d), " rows).")
  while (success < B && attempt < max_attempts) {
    attempt <- attempt + 1L
    fit <- tryCatch({
      db <- resample_boxes(d, attempt)
      ms <- suppressMessages(suppressWarnings(lmer(fs, db, REML = isREML(bundle$mod_swc))))
      mr <- suppressMessages(suppressWarnings(lmer(fr, db, REML = isREML(bundle$mod_resp))))
      list(paths = extract_paths(ms, mr, treatments), ss = isSingular(ms), sr = isSingular(mr))
    }, error = function(e) e)
    if (inherits(fit, "error")) {
      failures[[length(failures) + 1L]] <- tibble(species = row$species, resp_var = row$resp_var,
        model_id = row$model_id, attempt, error = conditionMessage(fit))
    } else {
      success <- success + 1L; sing_s <- sing_s + fit$ss; sing_r <- sing_r + fit$sr
      successes[[success]] <- fit$paths %>% mutate(species = row$species, resp_var = row$resp_var,
        response_label = row$response_label, model_id = row$model_id, replicate = success, attempt,
        estimate = if_else(component == "treatment_to_swc", estimate_raw, estimate_raw * orient),
        singular_swc = fit$ss, singular_response = fit$sr, .before = 1)
    }
    if (attempt %% 250L == 0L) message(row$model_id, ": ", success, "/", B)
  }
  if (success < B) return(empty_result(row, "bootstrap_incomplete",
    paste("Only", success, "successes after", attempt, "attempts."), sample, bal, scaling, fs_txt, fr_txt))
  fail <- bind_rows(failures)
  if (!nrow(fail)) fail <- tibble(species = character(), resp_var = character(), model_id = character(), attempt = integer(), error = character())
  stat <- tibble(species = row$species, resp_var = row$resp_var, response_label = row$response_label,
    model_id = row$model_id, status = "complete", n_rows_baseline = nrow(d0), n_rows_exact = nrow(d),
    retained_fraction = nrow(d) / nrow(d0), n_containers_exact = n_distinct(d$boxlabel),
    n_trees_exact = n_distinct(d$tree_id), n_boot_target = B, n_boot_success = success,
    n_attempts = attempt, n_failures = nrow(fail), n_singular_swc = sing_s,
    n_singular_response = sing_r, formula_swc = fs_txt, formula_response = fr_txt,
    seed = row$model_seed, source_rds_file = row$source_rds_file, note = NA_character_)
  list(reps = bind_rows(successes), point = point, failures = fail, sample = sample,
       balance = bal, scaling = scaling, status = stat)
}

rows <- split(inventory, seq_len(nrow(inventory)))
worker <- function(x) tryCatch(run_model(x), error = function(e) {
  row <- as.list(x); s <- tibble(model_id = row$model_id, n_rows_baseline = NA_integer_,
    n_rows_exact = NA_integer_, retained_fraction = NA_real_, n_containers_baseline = NA_integer_,
    n_containers_exact = NA_integer_, n_trees_baseline = NA_integer_, n_trees_exact = NA_integer_,
    n_dates_baseline = NA_integer_, n_dates_exact = NA_integer_, max_abs_swc_exact_difference = NA_real_)
  empty_result(row, "worker_error", conditionMessage(e), s)
})
if (.Platform$OS.type == "unix" && cores > 1L) {
  results <- parallel::mclapply(rows, worker, mc.cores = min(cores, length(rows)), mc.preschedule = FALSE)
} else results <- lapply(rows, worker)

reps <- map_dfr(results, "reps"); point <- map_dfr(results, "point")
failures <- map_dfr(results, "failures"); status <- map_dfr(results, "status")
samples <- map_dfr(results, "sample"); scales <- map_dfr(results, "scaling")
balance_long <- map_dfr(results, "balance")

effects <- reps %>% group_by(species, resp_var, response_label, model_id, treatment, component) %>%
  summarise(estimate_raw_boot_mean = mean(estimate_raw), estimate_boot_mean = mean(estimate),
    lower_raw = alinv_bootstrap_percentile(estimate_raw)[1], upper_raw = alinv_bootstrap_percentile(estimate_raw)[2],
    lower = alinv_bootstrap_percentile(estimate)[1], upper = alinv_bootstrap_percentile(estimate)[2],
    p_boot = alinv_bootstrap_p(estimate_raw), n_boot = n(), .groups = "drop") %>%
  left_join(point %>% select(species, resp_var, response_label, model_id, treatment, component,
                             estimate_raw, estimate, display_sign),
            by = c("species", "resp_var", "response_label", "model_id", "treatment", "component"))

# Existing baseline exports omit b; recover each replicate as indirect/a.
base_sign <- baseline_effects %>% distinct(model_id, display_sign)
base_b <- baseline_reps %>% select(model_id, replicate, treatment, component, estimate_raw) %>%
  pivot_wider(names_from = component, values_from = estimate_raw) %>%
  filter(abs(treatment_to_swc) > 1e-10) %>% mutate(b = indirect / treatment_to_swc) %>%
  group_by(model_id, replicate) %>% summarise(estimate_raw = median(b), .groups = "drop") %>%
  left_join(base_sign, by = "model_id") %>%
  mutate(estimate = estimate_raw * display_sign, raw_draw = estimate_raw, oriented_draw = estimate) %>%
  group_by(model_id) %>% summarise(treatment = "SWC", component = "swc_to_response",
    estimate_raw = median(raw_draw), estimate = median(oriented_draw),
    lower_raw = alinv_bootstrap_percentile(raw_draw)[1], upper_raw = alinv_bootstrap_percentile(raw_draw)[2],
    lower = alinv_bootstrap_percentile(oriented_draw)[1], upper = alinv_bootstrap_percentile(oriented_draw)[2],
    p_boot = alinv_bootstrap_p(raw_draw), n_boot = n(), .groups = "drop")
base <- bind_rows(baseline_effects %>% select(model_id, treatment, component, estimate_raw, estimate,
  lower_raw, upper_raw, lower, upper, p_boot, n_boot), base_b) %>%
  rename_with(~paste0("fuzzy_", .x), c(estimate_raw, estimate, lower_raw, upper_raw, lower, upper, p_boot, n_boot))
comparison <- effects %>% left_join(base, by = c("model_id", "treatment", "component")) %>%
  mutate(estimate_delta_exact_minus_fuzzy = estimate - fuzzy_estimate,
    raw_direction_agrees = sign(estimate_raw) == sign(fuzzy_estimate_raw),
    exact_significant = p_boot < .05, fuzzy_significant = fuzzy_p_boot < .05,
    significance_changed = exact_significant != fuzzy_significant)
balance_compare <- balance_long %>% select(model_id, sample, variable, level, n, proportion) %>%
  pivot_wider(names_from = sample, values_from = c(n, proportion), values_fill = 0) %>%
  mutate(proportion_change_exact_minus_fuzzy = proportion_exact_date - proportion_fuzzy_baseline)

qa <- bind_rows(
  point %>% filter(component %in% c("direct", "indirect", "total")) %>%
    select(model_id, treatment, component, estimate_raw) %>% pivot_wider(names_from = component, values_from = estimate_raw) %>%
    transmute(scope = "point", model_id, treatment, replicate = NA_integer_, error = total - direct - indirect),
  reps %>% filter(component %in% c("direct", "indirect", "total")) %>%
    select(model_id, treatment, replicate, component, estimate_raw) %>% pivot_wider(names_from = component, values_from = estimate_raw) %>%
    transmute(scope = "bootstrap", model_id, treatment, replicate, error = total - direct - indirect))
if (max(abs(qa$error), na.rm = TRUE) > 1e-12) stop("Path identity failed.")

plot_df <- comparison %>% filter(component %in% c("direct", "indirect", "total")) %>%
  select(species, response_label, treatment, component, exact_estimate = estimate,
    exact_lower = lower, exact_upper = upper, fuzzy_estimate, fuzzy_lower, fuzzy_upper) %>%
  pivot_longer(matches("^(exact|fuzzy)_"), names_to = c("matching", ".value"), names_pattern = "(exact|fuzzy)_(.*)") %>%
  mutate(matching = recode(matching, exact = "Exact date", fuzzy = "Current fuzzy match"),
    component = factor(component, c("direct", "indirect", "total"), c("Direct", "Indirect via SWC", "Total")))
p <- ggplot(plot_df, aes(estimate, treatment, colour = matching, shape = matching)) +
  geom_vline(xintercept = 0, colour = "grey65", linetype = 2) +
  geom_errorbar(aes(xmin = lower, xmax = upper), position = position_dodge(.45),
                width = 0, orientation = "y") +
  geom_point(position = position_dodge(.45), size = 1.7) +
  facet_grid(species + response_label ~ component, scales = "free_y", space = "free_y") +
  scale_colour_manual(values = c("Current fuzzy match" = "grey45", "Exact date" = "#0072B2")) +
  labs(x = "Oriented standardized path effect (95% percentile CI)", y = NULL, colour = NULL, shape = NULL) +
  theme_bw(9) + theme(legend.position = "bottom", strip.text.y = element_text(angle = 0), panel.grid.minor = element_blank())
ggsave(file.path(out, "exact-date-vs-fuzzy-sem-paths.pdf"), p, width = 11.5, height = 10.5)

status <- status %>% mutate(run_started_at = format(started), run_finished_at = format(Sys.time()),
  runtime_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")))
write_csv(inventory, file.path(out, "source-inventory.csv")); write_csv(status, file.path(out, "model-status.csv"))
write_csv(samples, file.path(out, "sample-loss.csv")); write_csv(scales, file.path(out, "baseline-standardization-constants.csv"))
write_csv(balance_compare, file.path(out, "treatment-balance.csv")); write_csv(point, file.path(out, "exact-date-point-effects.csv"))
write_csv(reps, file.path(out, "exact-date-bootstrap-replicates.csv")); write_csv(effects, file.path(out, "exact-date-bootstrap-effects.csv"))
write_csv(comparison, file.path(out, "exact-date-vs-fuzzy-comparison.csv")); write_csv(qa, file.path(out, "path-identity-qa.csv"))
write_csv(failures, file.path(out, "bootstrap-failures.csv"))
saveRDS(list(settings = list(B = B, cores = cores, matching = "exact measured container-date",
  scaling = "cached baseline z-values retained", formulas = "frozen; no simplification"), status = status,
  samples = samples, scales = scales, balance = balance_compare, point = point, effects = effects,
  comparison = comparison, qa = qa, reps = reps, failures = failures),
  file.path(out, "exact-date-sem-bootstrap-results.rds"))
print(status %>% select(model_id, status, n_rows_exact, retained_fraction, n_boot_success, note))
cat("Maximum identity error:", max(abs(qa$error), na.rm = TRUE), "\n")
