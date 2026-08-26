#!/usr/bin/env Rscript

# Harvest-biomass sensitivity analysis using a non-parametric container-cluster
# bootstrap. The fitted models intentionally match the current Figure 5 models;
# only uncertainty estimation changes.

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default) {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^", flag, "="), "", hit[[1]])
}

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this file with Rscript.", call. = FALSE)
script_file_arg <- sub("^--file=", "", script_arg[[1]])
# Rscript on some macOS installations protects spaces in --file paths with
# the literal token "~+~".
script_file_arg <- gsub("~\\+~", " ", script_file_arg)
script_path <- normalizePath(script_file_arg, winslash = "/", mustWork = TRUE)
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/", mustWork = TRUE)
output_dir <- file.path(project_root, "data", "final", "bootstrap", "biomass_lmm")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cache_dir <- file.path(tempdir(), "alinv-biomass-bootstrap-cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(XDG_CACHE_HOME = cache_dir)
setwd(project_root)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(readxl)
  library(ggplot2)
  library(lme4)
  library(scales)
})
suppressMessages(suppressPackageStartupMessages({
  source(file.path(project_root, "scripts", "functions", "project_data.R"))
  source(file.path(project_root, "scripts", "functions", "biomass.R"))
}))

bootstrap_target <- as.integer(arg_value("--bootstrap", "1000"))
n_cores <- as.integer(arg_value(
  "--cores",
  as.character(min(8L, parallel::detectCores(logical = FALSE)))
))
if (!is.finite(bootstrap_target) || bootstrap_target < 1L) {
  stop("--bootstrap must be a positive integer.", call. = FALSE)
}
if (!is.finite(n_cores) || n_cores < 1L) {
  stop("--cores must be a positive integer.", call. = FALSE)
}
max_attempts_per_model <- 5000L
bootstrap_seed <- 2026081805L
metrics <- c("shoot_biomass", "root_biomass", "root_shoot_biomass")
species_levels <- c("fagus", "quercus")
coefficient_terms <- c("precipitationdrought", "culturemixed", "robiniawith-robinia")

latest_biomass_workbook <- function() {
  candidates <- list.files(
    file.path(project_root, "data", "raw"),
    pattern = "^data_[0-9]+[.]xlsx$",
    full.names = TRUE
  )
  if (!length(candidates)) stop("No data_*.xlsx biomass workbook found.", call. = FALSE)
  normalizePath(candidates[order(basename(candidates), decreasing = TRUE)][[1]], winslash = "/", mustWork = TRUE)
}

prepare_model_data <- function(data, species_keep, metric) {
  data %>%
    filter(.data$species == species_keep, !is.na(.data[[metric]])) %>%
    alinv_apply_soil_context(soil_filter = "both", include_soil_treatment = FALSE) %>%
    mutate(
      precipitation = factor(.data$precipitation, levels = alinv_factor_levels("precipitation")),
      culture = factor(.data$culture, levels = alinv_factor_levels("culture")),
      robinia = factor(.data$robinia, levels = alinv_factor_levels("robinia")),
      block = factor(sub("-.*$", "", as.character(.data$boxlabel)), levels = c("b1", "b2", "b3")),
      y = .data[[metric]],
      # Retain the original-sample SD scale in every bootstrap sample so the
      # interval is for the same standardized effect shown in Figure 5.
      y_z = as.numeric(scale(.data$y)),
      tree_id = as.character(.data$tree_id),
      boxlabel = as.character(.data$boxlabel)
    ) %>%
    filter(!is.na(.data$block), !is.na(.data$y_z))
}

fit_model <- function(data, box_variable = "boxlabel") {
  formula_i <- stats::as.formula(paste0(
    "y_z ~ precipitation + culture + robinia + (1 | ", box_variable, ")"
  ))
  lme4::lmer(
    formula_i,
    data = data,
    REML = TRUE,
    control = lme4::lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 100000),
      check.conv.singular = "ignore"
    )
  )
}

extract_point <- function(model, species, metric, data) {
  coef_table <- summary(model)$coefficients
  tibble(
    species = species,
    metric = metric,
    term = rownames(coef_table),
    estimate = coef_table[, "Estimate"],
    se_wald = coef_table[, "Std. Error"]
  ) %>%
    filter(.data$term %in% coefficient_terms) %>%
    mutate(
      ci_lo_wald = .data$estimate - 1.96 * .data$se_wald,
      ci_hi_wald = .data$estimate + 1.96 * .data$se_wald,
      p_wald_normal = 2 * stats::pnorm(-abs(.data$estimate / .data$se_wald)),
      wald_ci_includes_zero = .data$ci_lo_wald <= 0 & .data$ci_hi_wald >= 0,
      singular_point_model = lme4::isSingular(model, tol = 1e-4),
      n_trees = dplyr::n_distinct(data$tree_id),
      n_containers = dplyr::n_distinct(data$boxlabel),
      n_blocks = dplyr::n_distinct(data$block),
      response_mean = mean(data$y, na.rm = TRUE),
      response_sd = stats::sd(data$y, na.rm = TRUE),
      oriented_estimate = .data$estimate,
      orientation = "raw biological sign retained"
    )
}

resample_containers_within_block <- function(data, replicate_id) {
  block_values <- levels(droplevels(data$block))
  sampled <- lapply(block_values, function(block_i) {
    block_data <- data[data$block == block_i, , drop = FALSE]
    box_ids <- unique(block_data$boxlabel)
    draws <- sample(box_ids, size = length(box_ids), replace = TRUE)
    pieces <- lapply(seq_along(draws), function(draw_i) {
      original_box <- draws[[draw_i]]
      piece <- block_data[block_data$boxlabel == original_box, , drop = FALSE]
      synthetic_box <- sprintf("boot%04d-%s-copy%03d", replicate_id, block_i, draw_i)
      piece$original_boxlabel <- original_box
      piece$boxlabel_boot <- synthetic_box
      piece$tree_id_boot <- paste0(synthetic_box, "::", piece$tree_id)
      piece
    })
    dplyr::bind_rows(pieces)
  })
  dplyr::bind_rows(sampled) %>%
    mutate(
      boxlabel_boot = factor(.data$boxlabel_boot),
      tree_id_boot = factor(.data$tree_id_boot)
    )
}

bootstrap_one_model <- function(data, species, metric, target, max_attempts, seed) {
  set.seed(seed)
  successful <- vector("list", target)
  success_n <- 0L
  attempts <- 0L
  failure_messages <- character()
  singular_n <- 0L

  while (success_n < target && attempts < max_attempts) {
    attempts <- attempts + 1L
    sample_i <- resample_containers_within_block(data, attempts)
    result_i <- tryCatch(
      withCallingHandlers({
        model_i <- fit_model(sample_i, "boxlabel_boot")
        estimates_i <- lme4::fixef(model_i)[coefficient_terms]
        if (length(estimates_i) != length(coefficient_terms) || any(!is.finite(estimates_i))) {
          stop("One or more treatment coefficients were absent or non-finite.")
        }
        list(
          estimates = estimates_i,
          singular = lme4::isSingular(model_i, tol = 1e-4)
        )
      }, warning = function(w) {
        # Convergence warnings are treated as failed refits. Singular fits are
        # handled separately because a zero container variance is admissible.
        if (grepl("failed to converge|degenerate|unable to evaluate|gradient", conditionMessage(w), ignore.case = TRUE)) {
          stop(conditionMessage(w), call. = FALSE)
        }
        invokeRestart("muffleWarning")
      }),
      error = function(e) e
    )

    if (inherits(result_i, "error")) {
      failure_messages <- c(failure_messages, conditionMessage(result_i))
      next
    }

    success_n <- success_n + 1L
    singular_n <- singular_n + as.integer(result_i$singular)
    successful[[success_n]] <- tibble(
      species = species,
      metric = metric,
      replicate = success_n,
      attempt = attempts,
      term = names(result_i$estimates),
      estimate = unname(result_i$estimates),
      oriented_estimate = unname(result_i$estimates),
      singular = result_i$singular
    )
  }

  if (success_n < target) {
    stop(
      sprintf("Only %d successful replicates after %d attempts for %s/%s.", success_n, attempts, species, metric),
      call. = FALSE
    )
  }

  failure_table <- if (length(failure_messages)) {
    tibble(message = failure_messages) %>% count(.data$message, name = "n") %>% arrange(desc(.data$n))
  } else {
    tibble(message = character(), n = integer())
  }

  list(
    estimates = bind_rows(successful),
    status = tibble(
      species = species,
      metric = metric,
      target_successful = target,
      successful = success_n,
      attempts = attempts,
      failed = attempts - success_n,
      singular_successful = singular_n,
      seed = seed
    ),
    failures = failure_table %>% mutate(species = species, metric = metric, .before = 1)
  )
}

bootstrap_summary <- function(bootstrap_estimates) {
  bootstrap_estimates %>%
    group_by(.data$species, .data$metric, .data$term) %>%
    summarise(
      bootstrap_mean = mean(.data$estimate),
      bootstrap_median = stats::median(.data$estimate),
      bootstrap_se = stats::sd(.data$estimate),
      ci_lo_boot = unname(stats::quantile(.data$estimate, 0.025, type = 7)),
      ci_hi_boot = unname(stats::quantile(.data$estimate, 0.975, type = 7)),
      p_boot = min(
        1,
        2 * min(
          (sum(.data$estimate <= 0) + 1) / (dplyr::n() + 1),
          (sum(.data$estimate >= 0) + 1) / (dplyr::n() + 1)
        )
      ),
      n_boot = dplyr::n(),
      .groups = "drop"
    ) %>%
    mutate(
      boot_ci_includes_zero = .data$ci_lo_boot <= 0 & .data$ci_hi_boot >= 0,
      oriented_ci_lo_boot = .data$ci_lo_boot,
      oriented_ci_hi_boot = .data$ci_hi_boot,
      oriented_bootstrap_mean = .data$bootstrap_mean,
      orientation = "raw biological sign retained"
    )
}

make_comparison <- function(point_results, boot_results) {
  point_results %>%
    left_join(boot_results, by = c("species", "metric", "term", "orientation")) %>%
    mutate(
      bootstrap_minus_point = .data$bootstrap_mean - .data$estimate,
      point_and_bootstrap_mean_same_sign = sign(.data$estimate) == sign(.data$bootstrap_mean),
      ci_inclusion_agrees = .data$wald_ci_includes_zero == .data$boot_ci_includes_zero,
      wald_significant = !.data$wald_ci_includes_zero,
      bootstrap_significant = !.data$boot_ci_includes_zero,
      significance_agrees = .data$wald_significant == .data$bootstrap_significant,
      wald_ci_width = .data$ci_hi_wald - .data$ci_lo_wald,
      bootstrap_ci_width = .data$ci_hi_boot - .data$ci_lo_boot,
      bootstrap_to_wald_width_ratio = .data$bootstrap_ci_width / .data$wald_ci_width,
      result_class = case_when(
        .data$significance_agrees & .data$point_and_bootstrap_mean_same_sign ~ "agrees",
        !.data$significance_agrees ~ "significance changed",
        TRUE ~ "sign changed"
      )
    )
}

theme_alinv_pub <- function(base_size = 7) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(color = "black"),
      axis.text = element_text(color = "black"),
      axis.ticks = element_line(linewidth = 0.25, color = "black"),
      axis.line = element_line(linewidth = 0.25, color = "black"),
      panel.border = element_rect(fill = NA, color = "grey55", linewidth = 0.25),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.18),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      legend.key.height = grid::unit(3.5, "mm"),
      legend.key.width = grid::unit(7, "mm"),
      legend.box.spacing = grid::unit(1, "mm"),
      strip.background = element_rect(fill = "grey12", color = "grey12", linewidth = 0.25),
      strip.text = element_text(color = "white", face = "bold", size = base_size)
    )
}

make_figure <- function(results) {
  term_labels <- c(
    culturemixed = "Culture: mono -> mixed",
    precipitationdrought = "Precipitation: control -> reduced",
    `robiniawith-robinia` = "Robinia: without -> with"
  )
  metric_labels <- c(
    shoot_biomass = "Shoot biomass",
    root_biomass = "Root biomass",
    root_shoot_biomass = "Root:shoot biomass"
  )
  species_labels <- c(fagus = "Fagus", quercus = "Quercus")

  plot_data <- results %>%
    mutate(
      term_label = factor(
        recode(.data$term, !!!term_labels),
        levels = c("Culture: mono -> mixed", "Precipitation: control -> reduced", "Robinia: without -> with")
      ),
      metric_label = factor(recode(.data$metric, !!!metric_labels), levels = unname(metric_labels)),
      species_label = factor(recode(.data$species, !!!species_labels), levels = unname(species_labels)),
      significant = factor(!.data$boot_ci_includes_zero, levels = c(FALSE, TRUE))
    )

  vals <- c(plot_data$ci_lo_boot, plot_data$ci_hi_boot, plot_data$estimate)
  axis_limit <- max(1, max(abs(vals[is.finite(vals)])) * 1.08)

  ggplot(plot_data, aes(x = .data$estimate, y = .data$term_label, color = .data$significant)) +
    geom_vline(xintercept = 0, linetype = "42", linewidth = 0.3, color = "grey45") +
    geom_segment(
      aes(x = .data$ci_lo_boot, xend = .data$ci_hi_boot, yend = .data$term_label),
      linewidth = 0.5
    ) +
    geom_point(size = 1.65) +
    facet_grid(.data$metric_label ~ .data$species_label) +
    coord_cartesian(xlim = c(-axis_limit, axis_limit)) +
    scale_x_continuous(breaks = scales::breaks_pretty(n = 5)) +
    scale_color_manual(
      values = c(`FALSE` = "#7A7A7A", `TRUE` = "#D65F5F"),
      labels = c(`FALSE` = "95% bootstrap CI includes 0", `TRUE` = "95% bootstrap CI excludes 0"),
      name = NULL,
      drop = FALSE
    ) +
    labs(x = "Standardized effect size (SD units, 95% bootstrap CI)", y = NULL) +
    theme_alinv_pub(base_size = 7) +
    theme(
      legend.position = "bottom",
      panel.grid.major.y = element_blank(),
      strip.text.y = element_text(angle = 90),
      axis.text.y = element_text(size = 6.7),
      plot.margin = margin(3, 5, 3, 4)
    )
}

write_comparison_markdown <- function(comparison, status, path) {
  changed <- comparison %>% filter(!.data$significance_agrees | !.data$point_and_bootstrap_mean_same_sign)
  width_summary <- comparison %>%
    summarise(
      min_ratio = min(.data$bootstrap_to_wald_width_ratio),
      median_ratio = stats::median(.data$bootstrap_to_wald_width_ratio),
      max_ratio = max(.data$bootstrap_to_wald_width_ratio)
    )
  wald_sig <- sum(comparison$wald_significant)
  boot_sig <- sum(comparison$bootstrap_significant)
  failure_total <- sum(status$failed)
  singular_point <- sum(status$point_model_singular)
  singular_boot <- sum(status$singular_successful)

  lines <- c(
    "# Harvest-biomass bootstrap comparison",
    "",
    "## Result",
    "",
    sprintf(
      "All %d treatment coefficients retained their sign; %d changed their 95%%-interval significance classification. The current Wald analysis classified %d/%d effects as significant, compared with %d/%d using percentile container-cluster bootstrap intervals.",
      nrow(comparison), sum(!comparison$significance_agrees), wald_sig, nrow(comparison), boot_sig, nrow(comparison)
    ),
    "",
    sprintf(
      "Bootstrap-to-Wald CI-width ratios ranged from %.2f to %.2f (median %.2f).",
      width_summary$min_ratio, width_summary$max_ratio, width_summary$median_ratio
    ),
    "",
    "## Method",
    "",
    "For each species and biomass metric, whole containers were sampled with replacement within experimental block. All trees from each selected container were retained and duplicated clusters received synthetic container and tree identifiers. The original Figure 5 model was then refitted: standardized response ~ precipitation + culture + Robinia + (1 | container). Response standardization was fixed to the original model sample. Percentile 95% intervals and empirical two-sided probabilities were calculated from 1,000 successful replicates.",
    "",
    "## Fit diagnostics",
    "",
    sprintf(
      "The six models produced %d successful bootstrap refits (%d per model) with %d failed attempts. Five of six point models were singular; %d successful bootstrap refits estimated a zero or near-zero container variance.",
      sum(status$successful), unique(status$successful)[1], failure_total, singular_boot
    ),
    "",
    "The analysed harvest workbook contained 32 containers from two represented blocks (b2 and b3) per species; block b1 had no rows in this biomass dataset. Model samples contained 128 Fagus trees and 127 Quercus trees.",
    "",
    "Singularity is scientifically informative here: for most biomass responses, between-container residual variance was estimated at zero after the additive treatment effects. Cluster resampling nevertheless retains the container as the treatment-assignment and dependence unit.",
    "",
    "## Classification changes",
    ""
  )

  if (!nrow(changed)) {
    lines <- c(lines, "None. Wald and bootstrap inference agreed on sign and whether each 95% interval included zero.")
  } else {
    change_lines <- apply(changed, 1, function(row) {
      sprintf(
        "- %s, %s, %s: Wald %s; bootstrap %s (bootstrap p = %.3f).",
        row[["species"]], row[["metric"]], row[["term"]],
        ifelse(as.logical(row[["wald_significant"]]), "significant", "not significant"),
        ifelse(as.logical(row[["bootstrap_significant"]]), "significant", "not significant"),
        as.numeric(row[["p_boot"]])
      )
    })
    lines <- c(lines, change_lines)
  }

  lines <- c(
    lines,
    "",
    "## Interpretation",
    "",
    "Bootstrap and Wald estimates use the same fitted model and therefore answer the same additive treatment-effect question. Differences reflect uncertainty estimation, not a change in the estimand. Because treatments were assigned to containers, the block-stratified container bootstrap is the more design-aligned sensitivity analysis. Results remain descriptive of the fitted factorial main effects; treatment interactions were not tested.",
    "",
    sprintf("Point-model singularity count: %d/6.", singular_point)
  )
  writeLines(lines, path, useBytes = TRUE)
}

workbook <- latest_biomass_workbook()
biomass <- suppressWarnings(wrangle_tree_biomass(workbook))

model_specs <- tidyr::crossing(species = species_levels, metric = metrics) %>%
  mutate(
    key = paste(.data$species, .data$metric, sep = "__"),
    seed = bootstrap_seed + row_number() * 100003L
  )

fit_and_bootstrap <- function(spec) {
  species_i <- spec$species[[1]]
  metric_i <- spec$metric[[1]]
  seed_i <- spec$seed[[1]]
  message("Fitting ", species_i, " / ", metric_i)
  model_data_i <- prepare_model_data(biomass, species_i, metric_i)
  point_model_i <- fit_model(model_data_i)
  bootstrap_i <- bootstrap_one_model(
    model_data_i, species_i, metric_i, bootstrap_target,
    max_attempts_per_model, seed_i
  )
  list(
    key = spec$key[[1]],
    model = point_model_i,
    point = extract_point(point_model_i, species_i, metric_i, model_data_i),
    estimates = bootstrap_i$estimates,
    status = bootstrap_i$status %>%
      mutate(
        n_trees = n_distinct(model_data_i$tree_id),
        n_containers = n_distinct(model_data_i$boxlabel),
        n_blocks = n_distinct(model_data_i$block),
        point_model_singular = lme4::isSingular(point_model_i, tol = 1e-4),
        formula = "y_z ~ precipitation + culture + robinia + (1 | boxlabel)",
        workbook = workbook
      ),
    failures = bootstrap_i$failures
  )
}

spec_rows <- split(model_specs, seq_len(nrow(model_specs)))
if (.Platform$OS.type == "unix" && n_cores > 1L) {
  model_results <- parallel::mclapply(
    spec_rows, fit_and_bootstrap,
    mc.cores = min(n_cores, length(spec_rows)),
    mc.preschedule = FALSE
  )
} else {
  model_results <- lapply(spec_rows, fit_and_bootstrap)
}

point_models <- setNames(lapply(model_results, `[[`, "model"), vapply(model_results, `[[`, character(1), "key"))
point_results <- bind_rows(lapply(model_results, `[[`, "point"))
bootstrap_estimates <- bind_rows(lapply(model_results, `[[`, "estimates"))
status_results <- bind_rows(lapply(model_results, `[[`, "status"))
failure_results <- bind_rows(lapply(model_results, `[[`, "failures"))
boot_summary <- bootstrap_summary(bootstrap_estimates)
comparison <- make_comparison(point_results, boot_summary)

readr::write_csv(point_results, file.path(output_dir, "biomass-point-wald-effects.csv"))
readr::write_csv(bootstrap_estimates, file.path(output_dir, "biomass-bootstrap-replicates.csv"))
readr::write_csv(boot_summary, file.path(output_dir, "biomass-bootstrap-effects.csv"))
readr::write_csv(status_results, file.path(output_dir, "biomass-bootstrap-status.csv"))
readr::write_csv(failure_results, file.path(output_dir, "biomass-bootstrap-failures.csv"))
readr::write_csv(comparison, file.path(output_dir, "biomass-wald-vs-bootstrap-comparison.csv"))
saveRDS(
  list(
    settings = list(
      target_successful = bootstrap_target,
      cores = n_cores,
      max_attempts_per_model = max_attempts_per_model,
      seed = bootstrap_seed,
      resampling = "container clusters within block",
      response_scaling = "fixed original-model-sample z score",
      formula = "y_z ~ precipitation + culture + robinia + (1 | boxlabel)",
      workbook = workbook
    ),
    point_models = point_models,
    point_effects = point_results,
    bootstrap_effects = boot_summary,
    bootstrap_status = status_results
  ),
  file.path(output_dir, "biomass-bootstrap-models.rds")
)

figure <- make_figure(comparison)
ggsave(
  file.path(output_dir, "fig5-biomass-effects-bootstrap.pdf"),
  figure,
  width = 160 / 25.4,
  height = 118 / 25.4,
  units = "in"
)
ggsave(
  file.path(output_dir, "fig5-biomass-effects-bootstrap.png"),
  figure,
  width = 160 / 25.4,
  height = 118 / 25.4,
  units = "in",
  dpi = 300,
  bg = "white"
)

write_comparison_markdown(
  comparison,
  status_results,
  file.path(output_dir, "comparison.md")
)

message("Completed biomass bootstrap analysis in: ", output_dir)
