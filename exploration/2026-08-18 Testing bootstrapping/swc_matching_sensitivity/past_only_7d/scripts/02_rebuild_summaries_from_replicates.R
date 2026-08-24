#!/usr/bin/env Rscript

# Rebuild deterministic summaries from already completed bootstrap replicates.
# This avoids rerunning the expensive model fits when only summary/export code
# changes. It is safe to run after 01_past_only_7d_sem_bootstrap.R.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path_arg <- gsub("~[+]~", " ", sub("^--file=", "", script_arg[[1]]), fixed = FALSE)
script_file <- normalizePath(script_path_arg, winslash = "/", mustWork = TRUE)
analysis_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/", mustWork = TRUE)
project_root <- normalizePath(file.path(analysis_dir, "..", "..", "..", ".."), winslash = "/", mustWork = TRUE)
output_dir <- file.path(analysis_dir, "output")
baseline_file <- file.path(
  project_root, "exploration", "2026-08-18 Testing bootstrapping",
  "repeated_response_sem", "output", "repeated-response-sem-bootstrap-effects.csv"
)

replicate_file <- file.path(output_dir, "past-only-sem-bootstrap-replicates.csv")
point_file <- file.path(output_dir, "past-only-sem-point-effects.csv")
rds_file <- file.path(output_dir, "past-only-sem-bootstrap-results.rds")
if (!all(file.exists(c(replicate_file, point_file, baseline_file, rds_file)))) {
  stop("Required completed bootstrap outputs were not found.")
}

point <- read_csv(point_file, show_col_types = FALSE)
display_signs <- point %>% distinct(model_id, display_sign)
replicates <- read_csv(replicate_file, show_col_types = FALSE) %>%
  select(-any_of("estimate")) %>%
  mutate(component = recode(component, a = "treatment_to_swc")) %>%
  left_join(display_signs, by = "model_id") %>%
  mutate(
    estimate = if_else(component == "treatment_to_swc", estimate_raw,
                       estimate_raw * display_sign)
  )

effects <- replicates %>%
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
    point %>%
      select(species, resp_var, response_label, model_id, treatment, component,
             estimate_raw, estimate, display_sign),
    by = c("species", "resp_var", "response_label", "model_id", "treatment", "component")
  ) %>%
  select(species, resp_var, response_label, model_id, treatment, component,
         estimate_raw, estimate, display_sign, everything())

baseline <- read_csv(baseline_file, show_col_types = FALSE)
comparison <- baseline %>%
  select(species, resp_var, response_label, model_id, treatment, component,
         estimate_baseline = estimate, lower_baseline = lower,
         upper_baseline = upper, p_boot_baseline = p_boot) %>%
  full_join(
    effects %>%
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

figure6_ready <- effects %>%
  filter(component == "total") %>%
  transmute(
    species, treatment, response_label, resp_var,
    estimate, lower, upper, p_boot, n_boot,
    lower_95 = lower, upper_95 = upper, n_boot_success = n_boot,
    estimate_raw, lower_raw, upper_raw,
    source = "past-only 7-day SWC; block-stratified container-cluster bootstrap"
  )

qa_total <- replicates %>%
  filter(component %in% c("direct", "indirect", "total")) %>%
  select(species, resp_var, model_id, replicate, treatment, component, estimate_raw) %>%
  pivot_wider(names_from = component, values_from = estimate_raw) %>%
  mutate(
    total_minus_direct_indirect = total - direct - indirect,
    qa_pass = is.finite(total_minus_direct_indirect) &
      abs(total_minus_direct_indirect) < 1e-10
  )
qa_summary <- qa_total %>%
  group_by(species, resp_var, model_id) %>%
  summarise(
    n_checks = n(), n_failures = sum(!qa_pass),
    max_abs_difference = max(abs(total_minus_direct_indirect), na.rm = TRUE),
    .groups = "drop"
  )

write_csv(replicates, replicate_file)
write_csv(effects, file.path(output_dir, "past-only-sem-bootstrap-effects.csv"))
write_csv(comparison, file.path(output_dir, "past-only-vs-fuzzy-bootstrap-comparison.csv"))
write_csv(figure6_ready, file.path(output_dir, "past-only-figure6-ready-total-effects.csv"))
write_csv(qa_summary, file.path(output_dir, "past-only-qa-total-identity.csv"))

results <- readRDS(rds_file)
results$replicates <- replicates
results$bootstrap_effects <- effects
results$comparison <- comparison
results$qa_total_identity <- qa_total
results$qa_summary <- qa_summary
saveRDS(results, rds_file)

cat("Past-only summaries rebuilt from completed replicates.\n")
cat("QA failures: ", sum(!qa_total$qa_pass), "\n", sep = "")
