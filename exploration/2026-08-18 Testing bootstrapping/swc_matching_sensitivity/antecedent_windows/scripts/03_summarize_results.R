#!/usr/bin/env Rscript

# Build the final branch-level scientific comparison from saved B=1000 outputs.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- gsub("~[+]~", " ", sub("^--file=", "", script_arg[[1]]))
script_file <- normalizePath(script_path, winslash = "/", mustWork = TRUE)
analysis_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
output_dir <- file.path(analysis_dir, "output")

comparison <- read_csv(
  file.path(output_dir, "antecedent-window-effect-comparisons.csv"),
  show_col_types = FALSE
)
effects <- read_csv(
  file.path(output_dir, "all-method-bootstrap-effects.csv"),
  show_col_types = FALSE
)
coverage <- read_csv(
  file.path(output_dir, "antecedent-window-coverage.csv"),
  show_col_types = FALSE
)
status <- read_csv(
  file.path(output_dir, "antecedent-window-sem-bootstrap-status.csv"),
  show_col_types = FALSE
)
identity <- read_csv(
  file.path(output_dir, "total-path-identity-qa.csv"),
  show_col_types = FALSE
)
reduced <- read_csv(
  file.path(output_dir, "reduced-form-total-qa.csv"),
  show_col_types = FALSE
)

primary <- comparison %>%
  filter(component %in% c("direct", "indirect", "total"))
overall_summary <- primary %>%
  summarise(
    n = sum(is.finite(a_estimate) & is.finite(b_estimate)),
    direction_agreement_pct = 100 * mean(direction_agrees, na.rm = TRUE),
    significance_agreement_pct = 100 * mean(significance_agrees, na.rm = TRUE),
    interval_overlap_pct = 100 * mean(intervals_overlap, na.rm = TRUE),
    median_abs_estimate_change = median(abs(estimate_difference_b_minus_a), na.rm = TRUE),
    .by = c(method_a, method_b)
  )
component_summary <- comparison %>%
  filter(component %in% c("direct", "indirect", "total", "swc_to_response")) %>%
  summarise(
    n = sum(is.finite(a_estimate) & is.finite(b_estimate)),
    direction_agreement_pct = 100 * mean(direction_agrees, na.rm = TRUE),
    significance_agreement_pct = 100 * mean(significance_agrees, na.rm = TRUE),
    interval_overlap_pct = 100 * mean(intervals_overlap, na.rm = TRUE),
    median_abs_estimate_change = median(abs(estimate_difference_b_minus_a), na.rm = TRUE),
    maximum_abs_estimate_change = max(abs(estimate_difference_b_minus_a), na.rm = TRUE),
    .by = c(method_a, method_b, component)
  )
method_path_summary <- effects %>%
  filter(component %in% c("treatment_to_swc", "swc_to_response", "direct", "indirect", "total")) %>%
  summarise(
    n = n(), n_significant = sum(p_boot < 0.05),
    median_abs_estimate = median(abs(estimate)),
    .by = c(method, component)
  )
indirect_share <- effects %>%
  filter(component %in% c("direct", "indirect", "total")) %>%
  select(method, species, resp_var, response_label, treatment, component, estimate) %>%
  pivot_wider(names_from = component, values_from = estimate) %>%
  mutate(indirect_share_abs_components = abs(indirect) / (abs(direct) + abs(indirect)))
indirect_share_summary <- indirect_share %>%
  summarise(
    n = n(), median_indirect_share = median(indirect_share_abs_components),
    mean_indirect_share = mean(indirect_share_abs_components),
    n_indirect_larger_than_direct = sum(indirect_share_abs_components > 0.5),
    .by = method
  )
b_comparison <- comparison %>%
  filter(component == "swc_to_response")

write_csv(overall_summary, file.path(output_dir, "antecedent-window-comparison-summary.csv"))
write_csv(component_summary, file.path(output_dir, "component-comparison-summary.csv"))
write_csv(method_path_summary, file.path(output_dir, "method-path-summary.csv"))
write_csv(indirect_share_summary, file.path(output_dir, "indirect-share-summary.csv"))
write_csv(b_comparison, file.path(output_dir, "swc-to-response-comparison.csv"))

get_comp <- function(method_b, component) {
  component_summary %>%
    filter(method_a == "fuzzy_baseline", .data$method_b == .env$method_b,
           .data$component == .env$component)
}
get_count <- function(method, component, field) {
  row <- method_path_summary %>%
    filter(.data$method == .env$method, .data$component == .env$component)
  row[[field]][[1]]
}
get_share <- function(method, field) {
  row <- indirect_share_summary %>% filter(.data$method == .env$method)
  row[[field]][[1]]
}

d7 <- get_comp("antecedent_7d", "direct")
i7 <- get_comp("antecedent_7d", "indirect")
t7 <- get_comp("antecedent_7d", "total")
b7 <- get_comp("antecedent_7d", "swc_to_response")
d14 <- get_comp("antecedent_14d", "direct")
i14 <- get_comp("antecedent_14d", "indirect")
t14 <- get_comp("antecedent_14d", "total")
b14 <- get_comp("antecedent_14d", "swc_to_response")

fmt <- function(x, digits = 1) format(round(x, digits), nsmall = digits, trim = TRUE)

lines <- c(
  "# Antecedent modeled-SWC sensitivity",
  "",
  "## Design",
  "",
  "Repeated-response SEMs were refitted using modeled trailing SWC means ending on each response date. The 7-day definition uses response day through day -6; the 14-day definition uses response day through day -13. The daily latent input is the existing climate-informed, treatment-agnostic container-level GAM estimate (`swc_hat`). Neither definition uses future information. These mediators are modeled antecedent means, not measured means.",
  "",
  "Frozen original formulas and model decisions were retained. Each model retained its fuzzy-baseline response observations, response scaling, and original SWC mean/SD, so comparisons isolate the mediator definition. Uncertainty used 1,000 successful block-stratified container-cluster bootstrap replicates per estimable species-response model and window.",
  "",
  "## Coverage and fitting",
  "",
  paste0(
    "All estimable baseline response records had complete 7- and 14-day modeled windows (minimum model-level coverage ",
    fmt(100 * min(coverage$coverage_proportion, na.rm = TRUE)), "%). The two Quercus senescence models remained unavailable because their original SEM bundles contained no post-filtering response data. All 24 estimable model-window tasks reached 1,000 successful replicates in 1,000 attempts; there were no failed refits."
  ),
  "",
  "The SWC submodels were never singular. Response-model singularity was frequent for several endpoints: 69.1% of Quercus total-volume replicates under the 7-day definition, 64.7% under 14 days, 60.6% for Fagus incremental volume under 14 days, and approximately half of Fagus vitality replicates. These are retained, as in the baseline bootstrap, but they indicate weak or boundary random-effect variance and require cautious interval interpretation.",
  "",
  "## Pathway comparison",
  "",
  paste0(
    "**Path-summed totals were effectively unchanged, but this is largely algebraic.** Relative to fuzzy matching, total-effect direction and interval overlap were both 100% for both antecedent windows. Median absolute changes were only ",
    fmt(t7$median_abs_estimate_change, 4), " (7 days) and ",
    fmt(t14$median_abs_estimate_change, 4), " (14 days), with maxima of ",
    fmt(t7$maximum_abs_estimate_change, 4), " and ", fmt(t14$maximum_abs_estimate_change, 4),
    " standardized units."
  ),
  "",
  paste0(
    "**The direct/indirect allocation was sensitive.** Direct-effect direction agreed with the fuzzy baseline for ",
    fmt(d7$direction_agreement_pct), "% of paths under 7 days and ",
    fmt(d14$direction_agreement_pct), "% under 14 days; corresponding interval overlap was ",
    fmt(d7$interval_overlap_pct), "% and ", fmt(d14$interval_overlap_pct), "%. Indirect-effect direction agreement was only ",
    fmt(i7$direction_agreement_pct), "% and ", fmt(i14$direction_agreement_pct), "%, and interval overlap was ",
    fmt(i7$interval_overlap_pct), "% and ", fmt(i14$interval_overlap_pct), "%."
  ),
  "",
  paste0(
    "The shared SWC-to-response path (`b`) was especially unstable: its direction agreed with the fuzzy baseline for ",
    fmt(b7$direction_agreement_pct), "% of the 12 response models under 7 days and only ",
    fmt(b14$direction_agreement_pct), "% under 14 days. Significant `b` paths numbered ",
    get_count("fuzzy_baseline", "swc_to_response", "n_significant"), " for fuzzy matching, ",
    get_count("antecedent_7d", "swc_to_response", "n_significant"), " for 7 days, and ",
    get_count("antecedent_14d", "swc_to_response", "n_significant"), " for 14 days. Significant indirect paths numbered ",
    get_count("fuzzy_baseline", "indirect", "n_significant"), ", ",
    get_count("antecedent_7d", "indirect", "n_significant"), ", and ",
    get_count("antecedent_14d", "indirect", "n_significant"), ", respectively."
  ),
  "",
  "Examples show why a single mechanistic conclusion is unsafe: the oriented Fagus chlorophyll `b` path changed from -0.198 under fuzzy matching to -0.074 under 7 days and +0.055 under 14 days; Quercus quantum yield changed from +0.076 to -0.083 and -0.054. These sign changes are conditional associations after DOY and treatment adjustment, not necessarily biological reversals.",
  "",
  paste0(
    "The broad magnitude statement that most modeled treatment effects were not allocated to the SWC-mediated path was nevertheless retained: median absolute indirect shares were ",
    fmt(100 * get_share("fuzzy_baseline", "median_indirect_share")), "% (fuzzy), ",
    fmt(100 * get_share("antecedent_7d", "median_indirect_share")), "% (7 days), and ",
    fmt(100 * get_share("antecedent_14d", "median_indirect_share")), "% (14 days). Indirect magnitude exceeded direct magnitude in ",
    get_share("fuzzy_baseline", "n_indirect_larger_than_direct"), ", ",
    get_share("antecedent_7d", "n_indirect_larger_than_direct"), ", and ",
    get_share("antecedent_14d", "n_indirect_larger_than_direct"), " of 48 paths. This supports only a broad descriptive statement, not stable attribution for individual treatment-response combinations."
  ),
  "",
  "## Regression-algebra check",
  "",
  paste0(
    "Stable path-summed totals are not independent evidence of robustness. In linear path models using the same observations and covariates, `c' + a*b` can be nearly invariant to mediator redefinition by regression algebra. The point path sum differed from the corresponding reduced-form treatment coefficient (response model with SWC omitted) by a median of ",
    fmt(median(reduced$abs_difference), 4), " standardized units; the 95th percentile was ",
    fmt(quantile(reduced$abs_difference, 0.95), 4), " and the maximum was ",
    fmt(max(reduced$abs_difference), 4), ". Mixed-model random-effects estimation prevents exact equality, but the near-identity confirms that decomposition sensitivity, not total-effect stability, is the meaningful test."
  ),
  "",
  "## Scientific recommendation",
  "",
  "Do not select a window because it produces more significant indirect paths. For a simple, pre-specified sensitivity analysis, the prior 7-day mean is the more defensible common window for rapidly responsive physiological traits; retain the 14-day mean as a longer-exposure check. Growth and senescence could plausibly integrate water availability over longer or stage-specific periods, but adopting endpoint-specific windows post hoc would add researcher degrees of freedom and should require an a priori biological rationale.",
  "",
  "More importantly, neither modeled window should be treated as definitive causal mediation. The daily GAM is climate-informed but treatment-agnostic, and its prediction uncertainty was not propagated. Ambient climate describes common temporal forcing but cannot reconstruct container-specific experimental watering or exclusion without dated treatment-operation records. Until a treatment-aware daily SWC reconstruction is available, report the total treatment patterns separately and describe direct-versus-SWC-mediated allocation as exploratory and sensitive to mediator timing.",
  "",
  "## QA",
  "",
  paste0(
    "All estimable point and bootstrap SEMs retained exact `total = direct + indirect`; the maximum absolute numerical error was ",
    format(max(identity$max_abs_identity_error), scientific = TRUE), "."
  )
)

writeLines(lines, file.path(analysis_dir, "comparison.md"))
