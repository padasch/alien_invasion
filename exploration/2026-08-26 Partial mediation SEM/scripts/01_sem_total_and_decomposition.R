#!/usr/bin/env Rscript

# Re-express the completed partial-mediation SEM bootstraps as one coherent
# analysis: total = SWC-adjusted path + SWC-associated indirect path.
#
# This script intentionally does not fit a second competing model. It rebuilds
# percentile intervals and bootstrap probabilities from the existing 1,000
# replicate-level draws, validates the path-sum identity, and produces the
# exploratory tables and comparison figure used to decide the final layout.

options(stringsAsFactors = FALSE)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this file with Rscript.", call. = FALSE)
script_file <- normalizePath(
  gsub("~\\+~", " ", sub("^--file=", "", script_arg[[1]])),
  winslash = "/",
  mustWork = TRUE
)
analysis_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/", mustWork = TRUE)
project_root <- normalizePath(file.path(analysis_dir, "..", ".."), winslash = "/", mustWork = TRUE)
output_dir <- file.path(analysis_dir, "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

renv_lib <- Sys.glob(file.path(project_root, "renv", "library", "*", "R-*", "*"))
if (length(renv_lib)) {
  .libPaths(c(normalizePath(renv_lib[[1]], winslash = "/", mustWork = TRUE), .libPaths()))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
})

n_boot_target <- 1000L
bootstrap_p <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  min(
    1,
    2 * min(
      (sum(x <= 0) + 1) / (length(x) + 1),
      (sum(x >= 0) + 1) / (length(x) + 1)
    )
  )
}

bootstrap_root <- file.path(project_root, "exploration", "2026-08-18 Testing bootstrapping")
repeated_dir <- file.path(
  bootstrap_root,
  "swc_matching_sensitivity", "past_only_7d", "output"
)
phenology_dir <- file.path(bootstrap_root, "phenology", "output")

source_paths <- c(
  repeated_replicates = file.path(repeated_dir, "past-only-sem-bootstrap-replicates.csv"),
  repeated_effects = file.path(repeated_dir, "past-only-sem-bootstrap-effects.csv"),
  repeated_status = file.path(repeated_dir, "past-only-sem-bootstrap-status.csv"),
  phenology_replicates = file.path(phenology_dir, "sem-block-stratified-bootstrap-replicates.csv"),
  phenology_effects = file.path(phenology_dir, "sem-effect-decomposition-block-stratified.csv"),
  phenology_status = file.path(phenology_dir, "phenology-bootstrap-model-status.csv")
)
missing_sources <- source_paths[!file.exists(source_paths)]
if (length(missing_sources)) {
  stop("Missing bootstrap source file(s): ", paste(missing_sources, collapse = ", "), call. = FALSE)
}

repeated_replicates <- read_csv(source_paths[["repeated_replicates"]], show_col_types = FALSE)
repeated_effects_source <- read_csv(source_paths[["repeated_effects"]], show_col_types = FALSE)
repeated_status <- read_csv(source_paths[["repeated_status"]], show_col_types = FALSE)
phenology_replicates <- read_csv(source_paths[["phenology_replicates"]], show_col_types = FALSE)
phenology_effects_source <- read_csv(source_paths[["phenology_effects"]], show_col_types = FALSE)
phenology_status <- read_csv(source_paths[["phenology_status"]], show_col_types = FALSE)

if (any(repeated_status$n_boot_success[repeated_status$status == "complete"] < n_boot_target)) {
  stop("At least one repeated-response SEM has fewer than 1,000 successful bootstrap refits.", call. = FALSE)
}
if (any(phenology_status$sem_boot_success < n_boot_target)) {
  stop("At least one phenology SEM has fewer than 1,000 successful bootstrap refits.", call. = FALSE)
}

# Recalculate the reported intervals and probabilities from the replicate-level
# draws. Point estimates remain the estimates from the original-data fits.
repeated_boot_summary <- repeated_replicates %>%
  filter(.data$component %in% c("direct", "indirect", "total")) %>%
  group_by(
    .data$species, .data$resp_var, .data$response_label,
    .data$treatment, .data$component
  ) %>%
  summarise(
    lower = quantile(.data$estimate, 0.025, names = FALSE, type = 7),
    upper = quantile(.data$estimate, 0.975, names = FALSE, type = 7),
    p_boot = bootstrap_p(.data$estimate),
    n_boot_success = n_distinct(.data$replicate),
    .groups = "drop"
  ) %>%
  left_join(
    repeated_effects_source %>%
      filter(.data$component %in% c("direct", "indirect", "total")) %>%
      select(
        "species", "resp_var", "response_label",
        "treatment", "component", "estimate"
      ),
    by = c("species", "resp_var", "response_label", "treatment", "component")
  ) %>%
  mutate(
    response_label = recode(.data$response_label, `Volume (incr.)` = "Volume (phase incr.)"),
    swc_definition = "Same-day or latest preceding measured SWC within 7 days",
    bootstrap_design = "Containers sampled with replacement within experimental block"
  )

phenology_boot_long <- phenology_replicates %>%
  select(
    "species", "replicate", "effect",
    "direct", "indirect", "total_path"
  ) %>%
  pivot_longer(
    cols = c("direct", "indirect", "total_path"),
    names_to = "component",
    values_to = "estimate_boot"
  ) %>%
  mutate(component = recode(.data$component, total_path = "total"))

phenology_point <- phenology_effects_source %>%
  filter(.data$metric %in% c("direct", "indirect", "total_path")) %>%
  transmute(
    species = .data$species,
    treatment = .data$effect,
    component = recode(.data$metric, total_path = "total"),
    estimate = .data$estimate
  )

phenology_boot_summary <- phenology_boot_long %>%
  group_by(.data$species, .data$effect, .data$component) %>%
  summarise(
    lower = quantile(.data$estimate_boot, 0.025, names = FALSE, type = 7),
    upper = quantile(.data$estimate_boot, 0.975, names = FALSE, type = 7),
    p_boot = bootstrap_p(.data$estimate_boot),
    n_boot_success = n_distinct(.data$replicate),
    .groups = "drop"
  ) %>%
  rename(treatment = "effect") %>%
  left_join(phenology_point, by = c("species", "treatment", "component")) %>%
  mutate(
    resp_var = "phenology_timing",
    response_label = "Leaf-out",
    swc_definition = "Container mean measured SWC from 4 March to 2 April 2025",
    bootstrap_design = "Containers sampled with replacement within experimental block"
  )

effects <- bind_rows(repeated_boot_summary, phenology_boot_summary) %>%
  mutate(
    component_label = recode(
      .data$component,
      total = "Total (c' + ab)",
      direct = "SWC-adjusted (c')",
      indirect = "Via SWC (ab)"
    ),
    significant = is.finite(.data$p_boot) & .data$p_boot < 0.05,
    interval_excludes_zero = .data$lower > 0 | .data$upper < 0,
    source_replicates = if_else(
      .data$resp_var == "phenology_timing",
      source_paths[["phenology_replicates"]],
      source_paths[["repeated_replicates"]]
    )
  ) %>%
  arrange(.data$species, .data$resp_var, .data$treatment, .data$component)

if (any(effects$n_boot_success < n_boot_target)) {
  stop("At least one displayed effect has fewer than 1,000 bootstrap draws.", call. = FALSE)
}
if (any(effects$significant != effects$interval_excludes_zero)) {
  warning("A bootstrap P-value and percentile-interval decision differ for at least one effect.")
}

# Validate total = direct + indirect within every bootstrap replicate.
repeated_identity <- repeated_replicates %>%
  filter(.data$component %in% c("direct", "indirect", "total")) %>%
  select(
    "species", "resp_var", "response_label",
    "replicate", "treatment", "component", "estimate"
  ) %>%
  pivot_wider(names_from = "component", values_from = "estimate") %>%
  mutate(discrepancy = .data$total - (.data$direct + .data$indirect)) %>%
  group_by(.data$species, .data$resp_var, .data$response_label) %>%
  summarise(
    bootstrap_replicates = n_distinct(.data$replicate),
    path_sum_checks = n(),
    max_abs_discrepancy = max(abs(.data$discrepancy), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(analysis = "Repeated response", swc_definition = "Past-only measured SWC, maximum lag 7 days")

phenology_identity <- phenology_replicates %>%
  mutate(discrepancy = .data$total_path - (.data$direct + .data$indirect)) %>%
  group_by(.data$species) %>%
  summarise(
    bootstrap_replicates = n_distinct(.data$replicate),
    path_sum_checks = n(),
    max_abs_discrepancy = max(abs(.data$discrepancy), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    resp_var = "phenology_timing",
    response_label = "Leaf-out",
    analysis = "Phenology",
    swc_definition = "Mean measured pre-leaf-out SWC, 4 March to 2 April 2025"
  )

identity_qa <- bind_rows(repeated_identity, phenology_identity) %>%
  select(
    "analysis", "species", "resp_var", "response_label",
    "swc_definition", "bootstrap_replicates",
    "path_sum_checks", "max_abs_discrepancy"
  )
if (max(identity_qa$max_abs_discrepancy, na.rm = TRUE) > 1e-10) {
  stop("Bootstrap total is not equal to direct + indirect within numerical tolerance.", call. = FALSE)
}

# A compact diagnostic for interpreting the decomposition. A percentage
# mediated is deliberately not calculated when paths oppose each other.
decomposition <- effects %>%
  select(
    "species", "resp_var", "response_label", "treatment",
    "component", "estimate", "lower", "upper",
    "p_boot", "significant"
  ) %>%
  pivot_wider(
    names_from = "component",
    values_from = c("estimate", "lower", "upper", "p_boot", "significant")
  ) %>%
  mutate(
    path_relation = case_when(
      .data$estimate_direct * .data$estimate_indirect < 0 ~ "opposing",
      .data$estimate_direct * .data$estimate_indirect > 0 ~ "same direction",
      TRUE ~ "one component is zero"
    ),
    larger_absolute_component = case_when(
      abs(.data$estimate_direct) > abs(.data$estimate_indirect) ~ "SWC-adjusted",
      abs(.data$estimate_direct) < abs(.data$estimate_indirect) ~ "via SWC",
      TRUE ~ "equal"
    ),
    indirect_fraction_of_absolute_sum = if_else(
      .data$path_relation == "same direction",
      abs(.data$estimate_indirect) /
        (abs(.data$estimate_direct) + abs(.data$estimate_indirect)),
      NA_real_
    )
  )

write_csv(effects, file.path(output_dir, "sem-effects-total-direct-indirect.csv"))
write_csv(decomposition, file.path(output_dir, "sem-decomposition-diagnostics.csv"))
write_csv(identity_qa, file.path(output_dir, "bootstrap-path-sum-qa.csv"))

# Plot contract: matrix/decomposition heatmap; rows are responses, columns are
# treatments, facets compare species and the three additive SEM components.
species_order <- c("fagus", "quercus")
response_order <- c(
  "Leaf-out", "Volume (total)", "Volume (phase incr.)", "Chlorophyll",
  "Vitality", "Quantum yield", "Senescence (%)", "Senescence (Chl)"
)
treatment_order <- c("precipitation", "robinia", "culture", "extreme_event")
component_order <- c("Total (c' + ab)", "SWC-adjusted (c')", "Via SWC (ab)")

plot_data <- expand_grid(
  species = species_order,
  response_label = response_order,
  treatment = treatment_order,
  component_label = component_order
) %>%
  left_join(
    effects,
    by = c("species", "response_label", "treatment", "component_label")
  ) %>%
  mutate(
    species_label = factor(.data$species, levels = species_order, labels = c("Fagus", "Quercus")),
    response_label = factor(.data$response_label, levels = rev(response_order)),
    treatment_label = factor(
      recode(
        .data$treatment,
        precipitation = "Precipitation",
        robinia = "Robinia",
        culture = "Culture",
        extreme_event = "Extreme event"
      ),
      levels = c("Precipitation", "Robinia", "Culture", "Extreme event")
    ),
    component_label = factor(.data$component_label, levels = component_order),
    missing = !is.finite(.data$estimate),
    coefficient_label = if_else(.data$significant, sprintf("%.2f", .data$estimate), "")
  )

fill_limit <- ceiling(max(abs(plot_data$estimate), na.rm = TRUE) * 10) / 10
figure <- ggplot(plot_data, aes(.data$treatment_label, .data$response_label)) +
  geom_tile(
    data = filter(plot_data, .data$missing),
    fill = "grey85",
    color = NA
  ) +
  geom_raster(
    data = filter(plot_data, !.data$missing),
    aes(fill = .data$estimate)
  ) +
  geom_text(
    data = filter(plot_data, .data$significant),
    aes(label = .data$coefficient_label),
    size = 2.4,
    color = "black"
  ) +
  facet_grid(.data$species_label ~ .data$component_label) +
  scale_fill_gradient2(
    low = "#D05A5A",
    mid = "white",
    high = "#547DA1",
    midpoint = 0,
    limits = c(-fill_limit, fill_limit),
    name = "Effect size"
  ) +
  scale_y_discrete(limits = rev(response_order), drop = FALSE) +
  labs(
    title = "SEM total and pathway decomposition",
    subtitle = "Total = SWC-adjusted + via SWC; coefficients are shown only when the 95% bootstrap interval excludes zero",
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 9) +
  theme(
    text = element_text(color = "black"),
    axis.text.x = element_text(angle = 32, hjust = 1, vjust = 1),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey12", color = NA),
    strip.text = element_text(color = "white", face = "bold"),
    legend.position = "right",
    plot.title = element_text(face = "plain", size = 11),
    plot.subtitle = element_text(size = 8.2),
    plot.margin = margin(5, 7, 5, 7)
  ) +
  guides(fill = guide_colorbar(barheight = grid::unit(35, "mm")))

ggsave(
  filename = file.path(output_dir, "sem-total-and-path-decomposition.pdf"),
  plot = figure,
  width = 250,
  height = 145,
  units = "mm",
  bg = "white",
  limitsize = FALSE
)

summary_counts <- effects %>%
  group_by(.data$species, .data$component) %>%
  summarise(n_effects = n(), n_significant = sum(.data$significant), .groups = "drop")
opposing_counts <- decomposition %>%
  group_by(.data$species) %>%
  summarise(
    n_effects = n(),
    n_opposing = sum(.data$path_relation == "opposing"),
    n_via_swc_larger = sum(.data$larger_absolute_component == "via SWC"),
    .groups = "drop"
  )

count_value <- function(data, species, component, column) {
  data %>%
    filter(.data$species == .env$species, .data$component == .env$component) %>%
    pull(all_of(column)) %>%
    first()
}
opposing_value <- function(species, column) {
  opposing_counts %>%
    filter(.data$species == .env$species) %>%
    pull(all_of(column)) %>%
    first()
}

summary_lines <- c(
  "# Partial-mediation SEM exploration", "",
  "## Question", "",
  "For each imposed treatment, separate the overall standardized response into an SWC-adjusted component and an SWC-associated indirect component, while retaining their sum as the total treatment effect.", "",
  "## Model", "",
  "For treatment k, the two linked equations are:", "",
  "`SWC_z = time terms + sum(a_k X_k) + random effects + error`", "",
  "`Y_z = b SWC_z + time terms + sum(c'_k X_k) + random effects + error`", "",
  "The SWC-adjusted component is `c'_k`, the SWC-associated indirect component is `a_k b`, and the path-summed total is `c'_k + a_k b`.", "",
  "Repeated responses use same-day or the latest preceding measured SWC within seven days. Leaf-out uses mean pre-leaf-out SWC at container level (4 March to 2 April 2025). Negative effects denote deterioration or delayed leaf-out; positive effects denote improvement or earlier leaf-out.", "",
  "## Bootstrap uncertainty", "",
  "Yes, bootstrapping is part of this analysis. Containers were sampled with replacement within experimental block, retaining all associated trees and observations. Each estimable model contributed 1,000 successful refits. Direct, indirect, and total effects were recalculated within every replicate; the confidence interval for the total therefore comes from the empirical distribution of `c' + ab`.", "",
  "The present exploration reuses the completed replicate-level fits and independently recalculates their percentile intervals and two-sided bootstrap probabilities. It does not perform a second set of identical model refits.", "",
  "## Numerical QA", "",
  paste0("All bootstrap path-sum checks passed: the maximum absolute difference between the stored total and `c' + ab` was ", format(max(identity_qa$max_abs_discrepancy), scientific = TRUE), "."), "",
  "## Main findings", "",
  paste0("- Fagus: ", count_value(summary_counts, "fagus", "direct", "n_significant"), " of ", count_value(summary_counts, "fagus", "direct", "n_effects"), " SWC-adjusted effects, ", count_value(summary_counts, "fagus", "indirect", "n_significant"), " of ", count_value(summary_counts, "fagus", "indirect", "n_effects"), " indirect effects, and ", count_value(summary_counts, "fagus", "total", "n_significant"), " of ", count_value(summary_counts, "fagus", "total", "n_effects"), " total effects excluded zero."),
  paste0("- Quercus: ", count_value(summary_counts, "quercus", "direct", "n_significant"), " of ", count_value(summary_counts, "quercus", "direct", "n_effects"), " SWC-adjusted effects, ", count_value(summary_counts, "quercus", "indirect", "n_significant"), " of ", count_value(summary_counts, "quercus", "indirect", "n_effects"), " indirect effects, and ", count_value(summary_counts, "quercus", "total", "n_significant"), " of ", count_value(summary_counts, "quercus", "total", "n_effects"), " total effects excluded zero."),
  paste0("- Direct and indirect estimates opposed one another in ", opposing_value("fagus", "n_opposing"), " of ", opposing_value("fagus", "n_effects"), " Fagus combinations and ", opposing_value("quercus", "n_opposing"), " of ", opposing_value("quercus", "n_effects"), " Quercus combinations."),
  paste0("- The indirect estimate was larger in absolute magnitude than the SWC-adjusted estimate in ", opposing_value("fagus", "n_via_swc_larger"), " Fagus and ", opposing_value("quercus", "n_via_swc_larger"), " Quercus combinations."), "",
  "## Interpretation", "",
  "This is a clean additive decomposition for presentation: the total panel answers what each treatment did overall, while the other panels show how that estimate splits into an SWC-adjusted remainder and an SWC-associated pathway. However, it does not support a blanket claim that Fagus responses were mainly mediated by measured SWC. Many Fagus indirect paths oppose the corresponding SWC-adjusted paths, particularly for chlorophyll, vitality, and quantum yield. Those patterns are better described as statistical suppression or competing associations than as straightforward positive water mediation.", "",
  "Because the response-SWC matching and the SWC-to-response path are observational within the experiment, `ab` should be called an SWC-associated indirect component rather than a proven causal mediation effect. A percentage mediated should only be reported when direct and indirect paths point in the same direction; the diagnostics table therefore leaves that quantity missing for opposing paths.", "",
  "## Files", "",
  "- `output/sem-total-and-path-decomposition.pdf`: shared-scale heatmap of total, SWC-adjusted, and indirect effects.",
  "- `output/sem-effects-total-direct-indirect.csv`: complete effect table with bootstrap intervals and probabilities.",
  "- `output/sem-decomposition-diagnostics.csv`: component dominance and path-direction diagnostics.",
  "- `output/bootstrap-path-sum-qa.csv`: replicate counts and numerical identity checks."
)
writeLines(summary_lines, file.path(analysis_dir, "analysis-summary.md"), useBytes = TRUE)

message("Wrote partial-mediation SEM exploration to: ", analysis_dir)
