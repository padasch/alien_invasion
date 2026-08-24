#!/usr/bin/env Rscript

# Temporal LMM container-cluster bootstrap exploration.
#
# This script intentionally writes only below exploration/2026-08-18 Testing
# bootstrapping/temporal_lmm. It preserves the manuscript model formulas and
# the response z-scaling calculated from the original species dataset. Only
# the uncertainty calculation changes.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(lme4)
  library(lubridate)
  library(patchwork)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

current_script <- function() {
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(arg)) stop("Run this file with Rscript.", call. = FALSE)
  path <- sub("^--file=", "", arg[[1]])
  # Rscript encodes spaces in --file paths as ~+~ on some macOS builds.
  path <- gsub("~\\+~", " ", path)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

SCRIPT_DIR <- dirname(current_script())
TEMPORAL_DIR <- normalizePath(file.path(SCRIPT_DIR, ".."), winslash = "/", mustWork = TRUE)
PROJECT_ROOT <- normalizePath(file.path(TEMPORAL_DIR, "..", "..", ".."), winslash = "/", mustWork = TRUE)
OUTPUT_DIR <- file.path(TEMPORAL_DIR, "output")
FIGURE_DIR <- file.path(OUTPUT_DIR, "figures")
DATA_DIR <- file.path(OUTPUT_DIR, "data")
MODEL_DIR <- file.path(OUTPUT_DIR, "models")
STATUS_DIR <- file.path(OUTPUT_DIR, "status")
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(MODEL_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(STATUS_DIR, recursive = TRUE, showWarnings = FALSE)

setwd(PROJECT_ROOT)
source(file.path(PROJECT_ROOT, "scripts", "auxiliary", "functions", "_source.R"))
source(file.path(PROJECT_ROOT, "scripts", "auxiliary", "functions", "3-effect-size-factorial.R"))

B_TARGET <- as.integer(Sys.getenv("ALINV_TEMPORAL_BOOT_B", unset = "1000"))
if (!is.finite(B_TARGET) || B_TARGET < 1L) stop("ALINV_TEMPORAL_BOOT_B must be a positive integer.", call. = FALSE)
detected_cores <- suppressWarnings(parallel::detectCores(logical = TRUE))
if (!length(detected_cores) || !is.finite(detected_cores)) detected_cores <- 2L
requested_cores <- as.integer(Sys.getenv("ALINV_TEMPORAL_BOOT_CORES", unset = "2"))
if (!is.finite(requested_cores) || requested_cores < 1L) requested_cores <- 1L
N_CORES <- as.integer(max(1L, min(requested_cores, detected_cores)))
BASE_SEED <- 2026081801L

MODEL_SPECS <- tibble::tribble(
  ~data_name,       ~resp_var,          ~species,  ~seed,
  "growth",        "diameter_inc_t0", "fagus",   BASE_SEED + 101L,
  "growth",        "diameter_inc_t0", "quercus", BASE_SEED + 102L,
  "quantum_yield", "qy",              "fagus",   BASE_SEED + 201L,
  "quantum_yield", "qy",              "quercus", BASE_SEED + 202L
)

TREATMENTS <- c("culture", "precipitation", "robinia")
BASELINES <- c(culture = "mono", precipitation = "control", robinia = "without-robinia")
TREATED <- c(culture = "mixed", precipitation = "drought", robinia = "with-robinia")
CONTRASTS <- c(
  culture = "mixed - mono",
  precipitation = "drought - control",
  robinia = "with-robinia - without-robinia"
)

temporal_formula <- y ~ date + date:robinia + date:precipitation + date:culture +
  (1 | boxlabel) + (1 | tree_id)

derive_block <- function(boxlabel) {
  block <- sub("^b([0-9]+)-.*$", "\\1", as.character(boxlabel))
  if (any(block == boxlabel) || any(!nzchar(block))) {
    stop("Could not derive block from every boxlabel.", call. = FALSE)
  }
  factor(block)
}

prepare_model_data <- function(data_name, resp_var, species) {
  out <- prepare_df_generic(
    type = "tree",
    data_name = data_name,
    resp_var = resp_var,
    species_keep = species,
    standardize_response = TRUE,
    add_covars = FALSE,
    soil_type = "both",
    include_soil_treatment = FALSE,
    swc_source = "measured"
  ) %>%
    mutate(
      block = derive_block(.data$boxlabel),
      date = factor(as.character(.data$date), levels = sort(unique(as.character(.data$date)))),
      culture = factor(.data$culture, levels = c("mono", "mixed")),
      precipitation = factor(.data$precipitation, levels = c("control", "drought")),
      robinia = factor(.data$robinia, levels = c("without-robinia", "with-robinia")),
      boxlabel = factor(.data$boxlabel),
      tree_id = factor(.data$tree_id)
    )

  if (any(!is.finite(out$y))) stop("Non-finite standardized response values.", call. = FALSE)
  out
}

fit_temporal <- function(df) {
  suppressMessages(suppressWarnings(
    lme4::lmer(
      temporal_formula,
      data = df,
      REML = TRUE,
      control = lme4::lmerControl(
        check.rankX = "stop.deficient",
        check.conv.singular = "ignore"
      )
    )
  ))
}

contrast_template <- function(df) {
  tidyr::expand_grid(
    date = as.Date(levels(df$date)),
    effect = TREATMENTS
  ) %>%
    mutate(
      contrast = unname(CONTRASTS[.data$effect]),
      coefficient = paste0(
        "date", format(.data$date, "%Y-%m-%d"), ":",
        .data$effect, unname(TREATED[.data$effect])
      ),
      key = paste(.data$date, .data$effect, sep = "||")
    )
}

extract_contrasts <- function(fit, template) {
  beta <- lme4::fixef(fit)
  if (!all(template$coefficient %in% names(beta))) {
    missing <- template$coefficient[!template$coefficient %in% names(beta)]
    stop("Missing fixed coefficients: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  stats::setNames(unname(beta[template$coefficient]), template$key)
}

fit_converged <- function(fit) {
  opt_code <- fit@optinfo$conv$opt %||% 0L
  messages <- fit@optinfo$conv$lme4$messages %||% character()
  identical(as.integer(opt_code), 0L) && length(messages) == 0L
}

bootstrap_design <- function(df) {
  row_index <- split(seq_len(nrow(df)), as.character(df$boxlabel))
  boxes_by_block <- split(
    unique(data.frame(boxlabel = as.character(df$boxlabel), block = as.character(df$block))),
    unique(data.frame(boxlabel = as.character(df$boxlabel), block = as.character(df$block)))$block
  )
  boxes_by_block <- lapply(boxes_by_block, function(x) as.character(x$boxlabel))
  list(row_index = row_index, boxes_by_block = boxes_by_block)
}

resample_containers <- function(df, design, replicate_id, seed) {
  set.seed(seed)
  selected <- unlist(lapply(design$boxes_by_block, function(boxes) {
    sample(boxes, size = length(boxes), replace = TRUE)
  }), use.names = FALSE)

  row_groups <- unname(design$row_index[selected])
  idx <- unlist(row_groups, use.names = FALSE)
  occurrence <- rep(seq_along(row_groups), lengths(row_groups))
  boot <- df[idx, , drop = FALSE]
  new_box <- sprintf("boot%04d_box%03d", replicate_id, occurrence)
  boot$tree_id <- factor(paste0(as.character(boot$tree_id), "__", new_box))
  boot$boxlabel <- factor(new_box)
  boot
}

run_one_boot <- function(attempt_id, seed, df, design, template) {
  tryCatch({
    boot_df <- resample_containers(df, design, attempt_id, seed)
    fit <- fit_temporal(boot_df)
    if (!fit_converged(fit)) {
      return(list(ok = FALSE, reason = "nonconverged", singular = lme4::isSingular(fit)))
    }
    estimates <- extract_contrasts(fit, template)
    if (length(estimates) != nrow(template) || any(!is.finite(estimates))) {
      return(list(ok = FALSE, reason = "incomplete_contrasts", singular = lme4::isSingular(fit)))
    }
    list(ok = TRUE, estimates = estimates, singular = lme4::isSingular(fit), reason = NA_character_)
  }, error = function(e) {
    list(ok = FALSE, reason = conditionMessage(e), singular = NA)
  })
}

run_bootstrap <- function(df, template, target, seed) {
  design <- bootstrap_design(df)
  set.seed(seed)
  attempt_seeds <- sample.int(.Machine$integer.max, target * 10L, replace = FALSE)
  successes <- list()
  failures <- list()
  attempt_start <- 1L

  while (length(successes) < target) {
    needed <- target - length(successes)
    # Attempt exactly the number still needed. This keeps the audit identity
    # attempts = accepted successes + recorded failures.
    batch_n <- min(needed, length(attempt_seeds) - attempt_start + 1L)
    if (batch_n <= 0L) stop("Exhausted bootstrap attempt seeds.", call. = FALSE)
    ids <- seq.int(attempt_start, length.out = batch_n)
    batch <- parallel::mclapply(
      ids,
      function(i) run_one_boot(i, attempt_seeds[[i]], df, design, template),
      mc.cores = N_CORES,
      mc.preschedule = TRUE,
      mc.set.seed = FALSE
    )
    for (j in seq_along(batch)) {
      item <- batch[[j]]
      item$attempt_id <- ids[[j]]
      if (isTRUE(item$ok) && length(successes) < target) {
        successes[[length(successes) + 1L]] <- item
      } else if (!isTRUE(item$ok)) {
        failures[[length(failures) + 1L]] <- item
      }
    }
    attempt_start <- max(ids) + 1L
  }

  draws <- do.call(rbind, lapply(successes, `[[`, "estimates"))
  rownames(draws) <- NULL
  singular_n <- sum(vapply(successes, function(x) isTRUE(x$singular), logical(1)))
  failure_reasons <- if (length(failures)) {
    sort(table(vapply(failures, `[[`, character(1), "reason")), decreasing = TRUE)
  } else {
    integer()
  }
  list(
    draws = draws,
    n_success = nrow(draws),
    n_attempts = attempt_start - 1L,
    n_failed = length(failures),
    n_singular = singular_n,
    failure_reasons = failure_reasons
  )
}

summarize_bootstrap <- function(draws, point, template) {
  stopifnot(nrow(draws) == B_TARGET, ncol(draws) == nrow(template))
  lower <- apply(draws, 2, stats::quantile, probs = 0.025, names = FALSE, type = 7)
  upper <- apply(draws, 2, stats::quantile, probs = 0.975, names = FALSE, type = 7)
  p_boot <- vapply(seq_len(ncol(draws)), function(j) {
    vals <- draws[, j]
    min(1, 2 * min(
      (sum(vals <= 0) + 1) / (length(vals) + 1),
      (sum(vals >= 0) + 1) / (length(vals) + 1)
    ))
  }, numeric(1))

  template %>%
    transmute(
      date, effect, contrast, key,
      estimate = unname(point[.data$key]),
      lower = unname(lower[.data$key]),
      upper = unname(upper[.data$key]),
      boot_mean = unname(colMeans(draws)[.data$key]),
      boot_median = unname(apply(draws, 2, stats::median)[.data$key]),
      boot_se = unname(apply(draws, 2, stats::sd)[.data$key]),
      p_boot = unname(p_boot[match(.data$key, names(point))]),
      display_sign = 1,
      estimate_raw = .data$estimate,
      lower_raw = .data$lower,
      upper_raw = .data$upper,
      inference = "container-cluster bootstrap percentile",
      bootstrap_replicates = B_TARGET
    )
}

as_draws_long <- function(draws, template, model_id) {
  tibble::as_tibble(draws) %>%
    mutate(bootstrap_id = dplyr::row_number(), .before = 1) %>%
    pivot_longer(-"bootstrap_id", names_to = "key", values_to = "estimate") %>%
    left_join(template %>% select("key", "date", "effect", "contrast"), by = "key") %>%
    mutate(model_id = model_id, .before = 1)
}

find_current_asymptotic <- function(data_name, resp_var, species) {
  basename_target <- paste0(
    "effect-tree-", data_name, "-", resp_var, "-", species,
    "-soil-both_without_soil_treatment-noCovars-swcMeas-effects.csv"
  )
  search_roots <- c(
    file.path(PROJECT_ROOT, "output"),
    file.path(PROJECT_ROOT, "_archive", "model-caches", "dated-output")
  )
  search_roots <- search_roots[dir.exists(search_roots)]
  files <- unlist(lapply(
    search_roots,
    list.files,
    recursive = TRUE,
    full.names = TRUE
  ), use.names = FALSE)
  files <- files[basename(files) == basename_target]
  if (!length(files)) stop("No current asymptotic effects file found: ", basename_target, call. = FALSE)
  date_match <- stringr::str_extract(files, "[0-9]{4}-[0-9]{2}-[0-9]{2}")
  ord <- order(as.Date(date_match), file.info(files)$mtime, decreasing = TRUE, na.last = TRUE)
  files[[ord[[1]]]]
}

run_model <- function(data_name, resp_var, species, seed) {
  model_id <- paste(data_name, resp_var, species, sep = "__")
  message("Fitting ", model_id, " with ", B_TARGET, " successful replicates...")
  df <- prepare_model_data(data_name, resp_var, species)
  template <- contrast_template(df)
  original_model <- fit_temporal(df)
  if (!fit_converged(original_model)) stop("Original model did not converge: ", model_id, call. = FALSE)
  point <- extract_contrasts(original_model, template)
  boot <- run_bootstrap(df, template, B_TARGET, seed)
  effects <- summarize_bootstrap(boot$draws, point, template) %>%
    mutate(
      model_id = model_id,
      data_name = data_name,
      response_var = resp_var,
      species = species,
      seed = seed,
      .before = 1
    )
  draws_long <- as_draws_long(boot$draws, template, model_id) %>%
    mutate(data_name = data_name, response_var = resp_var, species = species, seed = seed, .after = .data$model_id)

  current_file <- find_current_asymptotic(data_name, resp_var, species)
  current <- readr::read_csv(current_file, show_col_types = FALSE) %>%
    filter(.data$effect %in% TREATMENTS) %>%
    transmute(
      date = as.Date(.data$date), effect = as.character(.data$effect),
      current_contrast = .data$contrast,
      current_estimate = .data$estimate,
      current_se = .data$se,
      current_lower = .data$lower,
      current_upper = .data$upper,
      current_source_file = current_file
    )

  comparison <- effects %>%
    left_join(current, by = c("date", "effect")) %>%
    mutate(
      point_difference = .data$estimate - .data$current_estimate,
      lower_difference = .data$lower - .data$current_lower,
      upper_difference = .data$upper - .data$current_upper,
      current_width = .data$current_upper - .data$current_lower,
      bootstrap_width = .data$upper - .data$lower,
      width_ratio_boot_to_current = .data$bootstrap_width / .data$current_width,
      current_excludes_zero = .data$current_lower > 0 | .data$current_upper < 0,
      bootstrap_excludes_zero = .data$lower > 0 | .data$upper < 0,
      significance_agrees = .data$current_excludes_zero == .data$bootstrap_excludes_zero,
      sign_agrees = sign(.data$estimate) == sign(.data$current_estimate)
    )

  status <- tibble::tibble(
    model_id = model_id,
    data_name = data_name,
    response_var = resp_var,
    species = species,
    formula = paste(deparse(temporal_formula), collapse = " "),
    response_scaling = "z-scaled once on original species dataset; fixed in bootstrap",
    resampling = "containers with replacement within block; all nested rows retained",
    seed = seed,
    requested_successful_replicates = B_TARGET,
    successful_replicates = boot$n_success,
    attempts = boot$n_attempts,
    failed_attempts = boot$n_failed,
    singular_successful_fits = boot$n_singular,
    original_model_singular = lme4::isSingular(original_model),
    original_model_converged = fit_converged(original_model),
    n_obs = nrow(df),
    n_boxes = n_distinct(df$boxlabel),
    n_trees = n_distinct(df$tree_id),
    n_dates = n_distinct(df$date),
    n_blocks = n_distinct(df$block),
    current_asymptotic_source = current_file,
    failure_reasons = if (length(boot$failure_reasons)) {
      paste(names(boot$failure_reasons), as.integer(boot$failure_reasons), sep = "=", collapse = "; ")
    } else "none"
  )

  result <- list(
    model_id = model_id,
    formula = temporal_formula,
    seed = seed,
    settings = list(
      bootstrap_target = B_TARGET,
      resample_unit = "boxlabel",
      stratification = "block",
      interval = "percentile 95%",
      p_value = "empirical two-sided with +1 correction",
      response_scaling = "fixed original species-level z scale"
    ),
    original_model = original_model,
    original_effects = effects,
    bootstrap_draws = draws_long,
    status = status,
    comparison = comparison
  )
  saveRDS(result, file.path(MODEL_DIR, paste0(model_id, "-bootstrap.rds")), compress = "xz")
  list(effects = effects, draws = draws_long, status = status, comparison = comparison)
}

POSTPROCESS_ONLY <- tolower(Sys.getenv("ALINV_TEMPORAL_POSTPROCESS_ONLY", unset = "false")) %in% c("1", "true", "yes")
if (POSTPROCESS_ONLY) {
  message("Post-processing existing per-model bootstrap RDS files; no models will be refitted.")
  all_results <- lapply(seq_len(nrow(MODEL_SPECS)), function(i) {
    model_id <- paste(MODEL_SPECS$data_name[[i]], MODEL_SPECS$resp_var[[i]], MODEL_SPECS$species[[i]], sep = "__")
    result <- readRDS(file.path(MODEL_DIR, paste0(model_id, "-bootstrap.rds")))
    list(
      effects = result$original_effects,
      draws = result$bootstrap_draws,
      status = result$status,
      comparison = result$comparison
    )
  })
} else {
  all_results <- lapply(seq_len(nrow(MODEL_SPECS)), function(i) {
    run_model(
      MODEL_SPECS$data_name[[i]], MODEL_SPECS$resp_var[[i]],
      MODEL_SPECS$species[[i]], MODEL_SPECS$seed[[i]]
    )
  })
}

effects_all <- bind_rows(lapply(all_results, `[[`, "effects"))
draws_all <- bind_rows(lapply(all_results, `[[`, "draws"))
status_all <- bind_rows(lapply(all_results, `[[`, "status"))
comparison_all <- bind_rows(lapply(all_results, `[[`, "comparison"))

write_csv(effects_all, file.path(DATA_DIR, "temporal-lmm-bootstrap-effects.csv"))
write_csv(draws_all, file.path(DATA_DIR, "temporal-lmm-bootstrap-draws.csv"))
write_csv(comparison_all, file.path(DATA_DIR, "temporal-lmm-bootstrap-vs-current.csv"))
write_csv(status_all, file.path(STATUS_DIR, "temporal-lmm-bootstrap-status.csv"))
saveRDS(
  list(effects = effects_all, draws = draws_all, status = status_all, comparison = comparison_all),
  file.path(MODEL_DIR, "temporal-lmm-bootstrap-all.rds"), compress = "xz"
)

FIGURE_SELECTION <- tolower(Sys.getenv("ALINV_TEMPORAL_FIGURE", unset = "all"))
if (!FIGURE_SELECTION %in% c("all", "fig2", "fig3", "fig4")) {
  stop("ALINV_TEMPORAL_FIGURE must be one of: all, fig2, fig3, fig4.", call. = FALSE)
}

# Copy the descriptive Figure 2 exactly; it contains no model-based intervals.
if (FIGURE_SELECTION %in% c("all", "fig2")) {
  fig2_source <- file.path(PROJECT_ROOT, "output", "main_figures", "fig2_variation_timeseries.pdf")
  fig2_pdf <- file.path(FIGURE_DIR, "fig2_variation_timeseries.pdf")
  if (!file.copy(fig2_source, fig2_pdf, overwrite = TRUE)) stop("Could not copy Figure 2.", call. = FALSE)
  pdftoppm <- Sys.which("pdftoppm")
  if (!nzchar(pdftoppm)) stop("pdftoppm is required to render Figure 2 PNG.", call. = FALSE)
  status <- system2(pdftoppm, c("-png", "-singlefile", "-r", "300", shQuote(fig2_pdf), shQuote(file.path(FIGURE_DIR, "fig2_variation_timeseries"))))
  if (status != 0) stop("Failed to render Figure 2 PNG.", call. = FALSE)
}

EFFECT_COLORS <- c(culture = "#1B9E77", precipitation = "#D95F02", robinia = "#7570B3")
EFFECT_LABELS <- c(
  culture = "Culture (mono -> mixed)",
  precipitation = "Precipitation (control -> reduced)",
  robinia = "Robinia (without -> with)"
)
SUMMER_START <- as.Date("2025-06-20")
SUMMER_END <- as.Date("2025-09-01")
DROUGHT_WINDOWS <- tibble::tribble(
  ~start, ~end,
  as.Date("2025-06-20"), as.Date("2025-07-02"),
  as.Date("2025-08-12"), as.Date("2025-08-20")
)

theme_pub <- function(base_size = 7) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(family = "Helvetica", color = "black"),
      axis.text = element_text(color = "black"),
      axis.ticks = element_line(linewidth = 0.25, color = "black"),
      axis.line = element_line(linewidth = 0.25, color = "black"),
      panel.border = element_rect(fill = NA, color = "grey55", linewidth = 0.25),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.18),
      panel.grid.major.x = element_blank(), panel.grid.minor = element_blank()
    )
}

effect_panel <- function(df, title, limits, show_x) {
  drought_y <- limits[[1]] + diff(limits) * 0.04
  p <- ggplot(
    df %>% mutate(effect = factor(.data$effect, levels = TREATMENTS)),
    aes(x = .data$date, y = .data$estimate, ymin = .data$lower, ymax = .data$upper,
        color = .data$effect, fill = .data$effect, group = .data$effect)
  ) +
    geom_hline(yintercept = 0, linetype = "22", linewidth = 0.25, color = "grey35") +
    geom_vline(xintercept = SUMMER_START, linetype = "42", linewidth = 0.35) +
    geom_vline(xintercept = SUMMER_END, linetype = "42", linewidth = 0.35) +
    geom_ribbon(alpha = 0.16, color = NA) +
    geom_line(linewidth = 0.55) + geom_point(size = 0.8) +
    scale_color_manual(values = EFFECT_COLORS, breaks = TREATMENTS, labels = EFFECT_LABELS, name = NULL) +
    scale_fill_manual(values = EFFECT_COLORS, breaks = TREATMENTS, labels = EFFECT_LABELS, name = NULL) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b") +
    coord_cartesian(ylim = limits) +
    labs(x = NULL, y = "Effect size (SD units)", title = title) +
    theme_pub() +
    theme(
      axis.text.x = if (show_x) element_text() else element_blank(),
      axis.ticks.x = if (show_x) element_line(linewidth = 0.25) else element_blank(),
      # Generous explicit left padding protects facet titles and rotated
      # y-axis labels in PDF rasterizers that crop tightly to the media box.
      legend.position = "bottom", plot.margin = margin(2, 4, 2, 16)
    )
  for (i in seq_len(nrow(DROUGHT_WINDOWS))) {
    p <- p + annotate("segment", x = DROUGHT_WINDOWS$start[[i]], xend = DROUGHT_WINDOWS$end[[i]],
                      y = drought_y, yend = drought_y, color = "#E69F00", linewidth = 1.25, lineend = "round")
  }
  p
}

make_effect_figure <- function(resp_var, panel_suffix, filename) {
  dat <- effects_all %>% filter(.data$response_var == resp_var)
  vals <- c(dat$lower, dat$upper, dat$estimate, 0)
  rng <- range(vals[is.finite(vals)])
  span <- diff(rng)
  limits <- rng + c(-1, 1) * span * 0.08
  p1 <- effect_panel(filter(dat, .data$species == "fagus"), paste0("Fagus (", panel_suffix, ")"), limits, FALSE)
  p2 <- effect_panel(filter(dat, .data$species == "quercus"), paste0("Quercus (", panel_suffix, ")"), limits, TRUE)
  fig <- (p1 / p2) +
    plot_layout(guides = "collect") +
    plot_annotation(subtitle = "Current LMM estimates with 95% container-cluster bootstrap intervals (1,000 replicates)") &
    theme(
      plot.subtitle = element_text(size = 7, margin = margin(b = 3)),
      legend.position = "bottom",
      legend.text = element_text(size = 6.3),
      legend.key.width = grid::unit(7, "mm")
    )
  ggsave(file.path(FIGURE_DIR, paste0(filename, ".pdf")), fig, width = 160 / 25.4, height = 125 / 25.4, device = cairo_pdf)
  ggsave(file.path(FIGURE_DIR, paste0(filename, ".png")), fig, width = 160 / 25.4, height = 125 / 25.4, dpi = 300)
}

if (FIGURE_SELECTION %in% c("all", "fig3")) {
  make_effect_figure("diameter_inc_t0", "diameter incr.", "fig3_diameter_increment_effects")
}
if (FIGURE_SELECTION %in% c("all", "fig4")) {
  make_effect_figure("qy", "Fv/Fm", "fig4_quantum_yield_effects")
}

# Concise, reproducible numerical comparison report.
flip_rows <- comparison_all %>% filter(!.data$significance_agrees)
by_model <- comparison_all %>%
  group_by(.data$response_var, .data$species) %>%
  summarise(
    contrasts = n(),
    max_abs_point_difference = max(abs(.data$point_difference), na.rm = TRUE),
    median_width_ratio = median(.data$width_ratio_boot_to_current, na.rm = TRUE),
    min_width_ratio = min(.data$width_ratio_boot_to_current, na.rm = TRUE),
    max_width_ratio = max(.data$width_ratio_boot_to_current, na.rm = TRUE),
    significance_flips = sum(!.data$significance_agrees),
    sign_disagreements = sum(!.data$sign_agrees),
    .groups = "drop"
  )

fmt <- function(x, digits = 4) formatC(x, digits = digits, format = "fg", flag = "#")
report <- c(
  "# Temporal LMM bootstrap comparison",
  "",
  paste0("The comparison uses ", B_TARGET, " successful non-parametric container-cluster bootstrap replicates per species-response model. Containers were resampled within block, with all nested trees and dates retained. The response z-scaling and the current temporal LMM formula were held fixed."),
  "",
  "## Validation summary",
  "",
  paste0("- Successful models: ", sum(status_all$successful_replicates == B_TARGET), "/", nrow(status_all), "."),
  paste0("- Successful replicates: ", sum(status_all$successful_replicates), "; failed attempts: ", sum(status_all$failed_attempts), "; singular successful fits: ", sum(status_all$singular_successful_fits), "."),
  paste0("- Maximum absolute difference between refitted and saved current point estimates: ", fmt(max(abs(comparison_all$point_difference), na.rm = TRUE), 8), " SD."),
  paste0("- Point-estimate sign disagreements: ", sum(!comparison_all$sign_agrees), "."),
  paste0("- CI-based significance classifications changed for ", nrow(flip_rows), " of ", nrow(comparison_all), " date-specific contrasts."),
  paste0(
    "- Among those changes, ", sum(!flip_rows$current_excludes_zero & flip_rows$bootstrap_excludes_zero),
    " bootstrap intervals newly excluded zero and ", sum(flip_rows$current_excludes_zero & !flip_rows$bootstrap_excludes_zero),
    " no longer excluded zero."
  ),
  "",
  "## By response and species",
  "",
  "| Response | Species | Contrasts | Median width ratio | Width-ratio range | Significance flips |",
  "|---|---:|---:|---:|---:|---:|",
  vapply(seq_len(nrow(by_model)), function(i) paste0(
    "| ", by_model$response_var[[i]], " | ", by_model$species[[i]], " | ", by_model$contrasts[[i]],
    " | ", fmt(by_model$median_width_ratio[[i]], 3),
    " | ", fmt(by_model$min_width_ratio[[i]], 3), "–", fmt(by_model$max_width_ratio[[i]], 3),
    " | ", by_model$significance_flips[[i]], " |"
  ), character(1)),
  "",
  "## Changed CI classifications",
  ""
)

if (!nrow(flip_rows)) {
  report <- c(report, "None. Bootstrap and current asymptotic 95% intervals agree on whether zero is excluded for every date-specific contrast.")
} else {
  report <- c(
    report,
    "| Response | Species | Date | Effect | Current excludes zero | Bootstrap excludes zero | Current CI | Bootstrap CI |",
    "|---|---|---|---|---:|---:|---:|---:|",
    vapply(seq_len(nrow(flip_rows)), function(i) paste0(
      "| ", flip_rows$response_var[[i]], " | ", flip_rows$species[[i]], " | ", flip_rows$date[[i]],
      " | ", flip_rows$effect[[i]], " | ", flip_rows$current_excludes_zero[[i]],
      " | ", flip_rows$bootstrap_excludes_zero[[i]],
      " | [", fmt(flip_rows$current_lower[[i]], 3), ", ", fmt(flip_rows$current_upper[[i]], 3), "]",
      " | [", fmt(flip_rows$lower[[i]], 3), ", ", fmt(flip_rows$upper[[i]], 3), "] |"
    ), character(1))
  )
}

report <- c(
  report,
  "",
  "## Interpretation",
  "",
  "The two approaches estimate the same treatment contrasts: point estimates are unchanged apart from numerical fitting tolerance. Differences concern uncertainty only. The dominant temporal trajectories are therefore unchanged, and rows above are borderline inferential changes rather than reversals of the estimated biological direction. Date-specific intervals remain pointwise; no multiplicity correction across dates and treatments is applied.",
  "",
  paste0(
    sum(status_all$singular_successful_fits), " of ", sum(status_all$successful_replicates),
    " accepted bootstrap fits were singular, primarily the Quercus diameter model (",
    status_all$singular_successful_fits[status_all$response_var == "diameter_inc_t0" & status_all$species == "quercus"],
    "/1,000). These fits converged and returned complete contrasts, but at least one random-intercept variance was estimated at the boundary. This signals limited support for both random effects in some resamples and is a reason to interpret isolated date-specific CI flips cautiously."
  ),
  "",
  "## Files",
  "",
  "- `output/data/temporal-lmm-bootstrap-effects.csv`: plotted bootstrap summaries.",
  "- `output/data/temporal-lmm-bootstrap-draws.csv`: all successful replicate estimates.",
  "- `output/data/temporal-lmm-bootstrap-vs-current.csv`: row-level comparison.",
  "- `output/status/temporal-lmm-bootstrap-status.csv`: convergence and sample-size audit.",
  "- `output/models/`: complete per-model and combined RDS objects.",
  "- `output/figures/`: recreated Figures 2–4 in PDF and PNG."
)
writeLines(report, file.path(TEMPORAL_DIR, "comparison.md"))

message("Temporal bootstrap exploration completed: ", TEMPORAL_DIR)
