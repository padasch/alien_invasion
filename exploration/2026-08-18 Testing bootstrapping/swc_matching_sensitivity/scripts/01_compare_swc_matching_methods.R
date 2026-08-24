#!/usr/bin/env Rscript

# Cross-method synthesis for the repeated-response SWC matching sensitivity.
# This script does not fit models. It compares the completed B = 1,000
# bootstrap outputs from the independently reproducible method branches.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
})

arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- normalizePath(
  gsub("~[+]~", " ", sub("^--file=", "", arg[[1]])),
  winslash = "/", mustWork = TRUE
)
analysis_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
project_root <- normalizePath(file.path(analysis_dir, "..", "..", ".."), winslash = "/")
output_dir <- file.path(analysis_dir, "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_effects <- function(path, method, method_col = NULL) {
  if (!file.exists(path)) stop("Missing required input: ", path, call. = FALSE)
  x <- read_csv(path, show_col_types = FALSE)
  if (!is.null(method_col)) {
    x <- x %>% mutate(method = .data[[method_col]])
  } else {
    x <- x %>% mutate(method = method)
  }
  required <- c(
    "species", "resp_var", "response_label", "model_id", "treatment",
    "component", "estimate", "lower", "upper", "p_boot", "n_boot", "method"
  )
  missing <- setdiff(required, names(x))
  if (length(missing)) {
    stop("Input is missing columns [", paste(missing, collapse = ", "), "]: ", path)
  }
  x %>%
    select(all_of(required)) %>%
    mutate(method = as.character(method))
}

baseline <- read_effects(
  file.path(
    project_root, "exploration", "2026-08-18 Testing bootstrapping",
    "repeated_response_sem", "output", "repeated-response-sem-bootstrap-effects.csv"
  ),
  "fuzzy_baseline"
)
past <- read_effects(
  file.path(analysis_dir, "past_only_7d", "output", "past-only-sem-bootstrap-effects.csv"),
  "past_only_7d"
)
exact <- read_effects(
  file.path(analysis_dir, "exact_date", "output", "exact-date-bootstrap-effects.csv"),
  "exact_date"
)
daily <- read_effects(
  file.path(
    analysis_dir, "daily_interpolation", "output",
    "daily-interpolation-sem-bootstrap-effects.csv"
  ),
  "daily_interpolation"
)
antecedent <- read_effects(
  file.path(
    analysis_dir, "antecedent_windows", "output",
    "antecedent-window-sem-bootstrap-effects.csv"
  ),
  method = NA_character_, method_col = "method"
)

effects <- bind_rows(baseline, past, exact, daily, antecedent) %>%
  mutate(
    method = recode(
      method,
      fuzzy_baseline = "Fuzzy baseline",
      past_only_7d = "Past-only measured (0-7 d)",
      exact_date = "Exact-date measured",
      daily_interpolation = "Exact-date daily GAM",
      antecedent_7d = "Antecedent GAM mean (7 d)",
      antecedent_14d = "Antecedent GAM mean (14 d)"
    ),
    significant = is.finite(p_boot) & p_boot < 0.05
  )

bad_boot <- effects %>% filter(is.finite(n_boot), n_boot != 1000)
if (nrow(bad_boot)) {
  stop(
    "One or more completed effect rows do not contain 1,000 bootstrap replicates. ",
    "The sensitivity synthesis must not use smoke-test outputs.", call. = FALSE
  )
}

method_order <- c(
  "Fuzzy baseline", "Past-only measured (0-7 d)", "Exact-date measured",
  "Exact-date daily GAM", "Antecedent GAM mean (7 d)",
  "Antecedent GAM mean (14 d)"
)
effects <- effects %>% mutate(method = factor(method, levels = method_order))
write_csv(effects, file.path(output_dir, "all-swc-matching-bootstrap-effects.csv"))

key <- c("species", "resp_var", "response_label", "model_id", "treatment", "component")
baseline_ref <- effects %>%
  filter(method == "Fuzzy baseline") %>%
  select(all_of(key), estimate_baseline = estimate, lower_baseline = lower,
         upper_baseline = upper, significant_baseline = significant)

comparison <- effects %>%
  filter(method != "Fuzzy baseline") %>%
  left_join(baseline_ref, by = key) %>%
  mutate(
    estimate_delta = estimate - estimate_baseline,
    direction_agrees = sign(estimate) == sign(estimate_baseline),
    significance_agrees = significant == significant_baseline,
    intervals_overlap = pmax(lower, lower_baseline) <= pmin(upper, upper_baseline)
  )
write_csv(comparison, file.path(output_dir, "swc-matching-vs-fuzzy-comparison.csv"))

comparison_summary <- comparison %>%
  filter(component %in% c("direct", "indirect", "total")) %>%
  group_by(method, component) %>%
  summarise(
    n = n(),
    direction_agreement_pct = 100 * mean(direction_agrees, na.rm = TRUE),
    significance_agreement_pct = 100 * mean(significance_agrees, na.rm = TRUE),
    interval_overlap_pct = 100 * mean(intervals_overlap, na.rm = TRUE),
    median_absolute_estimate_change = median(abs(estimate_delta), na.rm = TRUE),
    mean_absolute_estimate_change = mean(abs(estimate_delta), na.rm = TRUE),
    .groups = "drop"
  )
write_csv(comparison_summary, file.path(output_dir, "swc-matching-comparison-summary.csv"))

dominance <- effects %>%
  filter(component %in% c("direct", "indirect")) %>%
  select(method, all_of(key), component, estimate) %>%
  pivot_wider(names_from = component, values_from = estimate) %>%
  filter(is.finite(direct), is.finite(indirect)) %>%
  group_by(method) %>%
  summarise(
    n = n(),
    direct_larger_n = sum(abs(direct) > abs(indirect)),
    indirect_larger_n = sum(abs(indirect) > abs(direct)),
    direct_larger_pct = 100 * mean(abs(direct) > abs(indirect)),
    indirect_magnitude_share =
      sum(abs(indirect)) / (sum(abs(direct)) + sum(abs(indirect))),
    opposing_sign_pct = 100 * mean(sign(direct) != sign(indirect)),
    .groups = "drop"
  )
write_csv(dominance, file.path(output_dir, "swc-matching-pathway-dominance.csv"))

plot_data <- comparison %>%
  filter(component %in% c("direct", "indirect", "total")) %>%
  mutate(
    component = factor(
      component, levels = c("direct", "indirect", "total"),
      labels = c("Direct", "Indirect via SWC", "Total")
    ),
    method = droplevels(method),
    response = factor(
      paste(if_else(species == "fagus", "Fagus", "Quercus"), response_label, sep = ": "),
      levels = rev(unique(paste(
        rep(c("Fagus", "Quercus"), each = 7),
        rep(c("Volume (total)", "Volume (incr.)", "Chlorophyll", "Vitality",
              "Quantum yield", "Senescence (%)", "Senescence (Chl)"), 2), sep = ": "
      )))
    ),
    treatment = factor(
      treatment,
      levels = c("precipitation", "robinia", "culture", "extreme_event"),
      labels = c("Precipitation", "Robinia", "Culture", "Extreme event")
    )
  )

fill_limit <- max(abs(plot_data$estimate_delta), na.rm = TRUE)
p <- ggplot(plot_data, aes(treatment, response, fill = estimate_delta)) +
  geom_tile(color = "white", linewidth = 0.2) +
  facet_grid(method ~ component, scales = "free_y", space = "free_y", drop = FALSE) +
  scale_fill_gradient2(
    low = "#D65F5F", mid = "white", high = "#3C6E8F", midpoint = 0,
    limits = c(-fill_limit, fill_limit),
    name = expression(Delta~"effect vs fuzzy")
  ) +
  labs(
    title = "SWC matching sensitivity of SEM effects",
    subtitle = paste(
      "Coefficient change relative to the fuzzy-match baseline;",
      "red is more negative and blue more positive. Missing cells were not estimable."
    ),
    x = NULL, y = NULL
  ) +
  theme_classic(base_size = 8) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8, color = "grey30", margin = margin(b = 6)),
    axis.text.x = element_text(angle = 35, hjust = 1),
    axis.ticks = element_blank(), axis.line = element_blank(),
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(
  file.path(output_dir, "swc-matching-sensitivity-summary.pdf"), p,
  width = 210, height = 285, units = "mm", device = cairo_pdf,
  bg = "white", limitsize = FALSE
)

message("Created cross-method SWC matching sensitivity synthesis in: ", output_dir)
