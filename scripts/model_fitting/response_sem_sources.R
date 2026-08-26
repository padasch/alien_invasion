#!/usr/bin/env Rscript

# Rebuild the response-model source samples used by the repeated-response and
# reduced-form bootstraps from the current tracked interim data. The additive
# manuscript formulas are specified here explicitly so the active analysis no
# longer depends on archived model-cache files.

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(lme4)
  library(purrr)
  library(readr)
  library(tibble)
  library(tidyr)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- gsub("~[+]~", " ", sub("^--file=", "", script_arg[[1]]))
script_file <- normalizePath(script_path, winslash = "/", mustWork = TRUE)
project_root <- normalizePath(file.path(dirname(script_file), "..", ".."), winslash = "/", mustWork = TRUE)
source_root <- file.path(
  project_root, "data", "final", "bootstrap", "repeated_response_sem",
  "source_models", format(Sys.Date(), "%Y-%m-%d")
)
dir.create(source_root, recursive = TRUE, showWarnings = FALSE)
setwd(project_root)

source(file.path(project_root, "scripts", "functions", "project_data.R"))
source(file.path(project_root, "scripts", "functions", "factorial_effects.R"))

response_specs <- tibble::tribble(
  ~data_name, ~resp_var, ~response_label,
  "growth", "volume", "Volume (total)",
  "growth", "volume_inc_phase_rel", "Volume (incr.)",
  "chlorophyll", "chl", "Chlorophyll",
  "condition", "condition", "Vitality",
  "quantum_yield", "qy", "Quantum yield",
  "senescence", "remaining_green", "Senescence (%)",
  "senescence", "chlavg", "Senescence (Chl)"
)

model_grid <- tidyr::crossing(
  species = c("fagus", "quercus"),
  response_specs
) |>
  dplyr::mutate(model_id = paste(.data$species, .data$resp_var, sep = "__"))

prepare_current_sem_data <- function(data_name, resp_var, species) {
  df <- prepare_df_generic(
    type = "tree",
    data_name = data_name,
    resp_var = resp_var,
    species_keep = species,
    standardize_response = FALSE,
    add_covars = FALSE,
    soil_type = "both",
    include_soil_treatment = FALSE
  ) |>
    dplyr::filter(!is.na(.data$y), !is.na(.data$swc), !is.na(.data$date)) |>
    dplyr::mutate(
      date = as.Date(as.character(.data$date)),
      tree_id = factor(.data$tree_id),
      boxlabel = factor(.data$boxlabel),
      robinia = factor(.data$robinia, levels = alinv_factor_levels("robinia")),
      extreme_event = factor(.data$extreme_event, levels = alinv_factor_levels("extreme_event")),
      precipitation = factor(.data$precipitation, levels = alinv_factor_levels("precipitation")),
      culture = factor(.data$culture, levels = alinv_factor_levels("culture")),
      soiltype = alinv_relevel_soiltype(.data$soiltype),
      doy = as.numeric(format(.data$date, "%j")),
      doy_c_raw = as.numeric(scale(.data$doy, center = TRUE, scale = FALSE)),
      swc_org = .data$swc,
      y_org = .data$y,
      swc = as.numeric(scale(.data$swc_org)),
      y = as.numeric(scale(.data$y_org)),
      doy_c = as.numeric(scale(.data$doy_c_raw)),
      doy_c2 = .data$doy_c^2
    ) |>
    dplyr::select(
      .data$tree_id, .data$boxlabel, .data$date,
      .data$swc, .data$swc_org, .data$y, .data$y_org,
      .data$robinia, .data$extreme_event, .data$precipitation,
      .data$culture, .data$soiltype, .data$doy, .data$doy_c, .data$doy_c2
    ) |>
    tidyr::drop_na()
  droplevels(df)
}

coefficient_info <- function(model, term) {
  beta <- lme4::fixef(model)
  hits <- if (identical(term, "swc")) {
    which(names(beta) == "swc")
  } else {
    which(grepl(term, names(beta), fixed = TRUE) & !grepl(":", names(beta), fixed = TRUE))
  }
  if (length(hits) != 1L) {
    return(list(est = NA_real_, se = NA_real_, p = NA_real_))
  }
  name <- names(beta)[hits[[1]]]
  est <- unname(beta[[name]])
  se <- sqrt(as.matrix(stats::vcov(model))[name, name])
  list(est = est, se = se, p = 2 * stats::pnorm(-abs(est / se)))
}

drop_collapsed_formula_terms <- function(formula, data) {
  collapsed <- names(data)[vapply(
    data,
    function(x) (is.factor(x) || is.character(x)) && dplyr::n_distinct(x, na.rm = TRUE) < 2L,
    logical(1)
  )]
  updated <- formula
  for (variable in collapsed) {
    updated <- stats::update(updated, stats::as.formula(paste(". ~ . -", variable)))
  }
  updated
}

extract_effects <- function(mod_swc, mod_resp, treatments) {
  b_info <- coefficient_info(mod_resp, "swc")
  purrr::map_dfr(treatments, function(treatment) {
    a_info <- coefficient_info(mod_swc, treatment)
    c_info <- coefficient_info(mod_resp, treatment)
    indirect_value <- a_info$est * b_info$est
    se_indirect <- sqrt((b_info$est^2) * (a_info$se^2) + (a_info$est^2) * (b_info$se^2))
    total_value <- c_info$est + indirect_value
    se_total <- sqrt(c_info$se^2 + se_indirect^2)
    tibble::tibble(
      factor = treatment,
      a = a_info$est, se_a = a_info$se, p_a = a_info$p,
      b = b_info$est, se_b = b_info$se, p_b = b_info$p,
      c_direct = c_info$est, se_c = c_info$se, p_c = c_info$p,
      indirect = indirect_value, se_ind = se_indirect,
      p_ind = 2 * stats::pnorm(-abs(indirect_value / se_indirect)),
      total = total_value, se_tot = se_total,
      p_tot = 2 * stats::pnorm(-abs(total_value / se_total))
    )
  })
}

make_matrix_data <- function(effects) {
  dplyr::bind_rows(
    effects |> dplyr::transmute(factor = .data$factor, path_type = "treatment_to_swc", estimate_raw = .data$a),
    effects |> dplyr::transmute(factor = .data$factor, path_type = "direct", estimate_raw = .data$c_direct),
    effects |> dplyr::transmute(factor = .data$factor, path_type = "indirect", estimate_raw = .data$indirect),
    effects |> dplyr::transmute(factor = .data$factor, path_type = "total", estimate_raw = .data$total)
  ) |>
    dplyr::mutate(estimate = .data$estimate_raw)
}

refresh_one <- function(row) {
  species <- row$species[[1]]
  data_name <- row$data_name[[1]]
  resp_var <- row$resp_var[[1]]
  response_label <- row$response_label[[1]]
  model_id <- row$model_id[[1]]
  source_stem <- paste0(
    "sem-tree-", data_name, "-", resp_var, "-", species,
    "-soil-both_without_soil_treatment-noInt-scaled-rfeAIC2-all-swcMeas"
  )
  source_rds_file <- file.path(source_root, paste0(source_stem, ".rds"))
  matrix_file <- sub("[.]rds$", "-matrix_data.csv", source_rds_file)

  data <- tryCatch(
    prepare_current_sem_data(data_name, resp_var, species),
    error = function(e) e
  )
  if (inherits(data, "error") || !nrow(data)) {
    note <- if (inherits(data, "error")) conditionMessage(data) else {
      "No current response rows were available after model-data filtering."
    }
    saveRDS(list(data = tibble::tibble(), note = note), source_rds_file)
    readr::write_csv(tibble::tibble(path_type = character()), matrix_file)
    return(tibble::tibble(
      species, resp_var, response_label, model_id, status = "unavailable",
      n_rows = 0L, source_rds_file, note
    ))
  }

  formula_swc <- drop_collapsed_formula_terms(
    swc ~ doy_c + doy_c2 + precipitation + robinia + culture + extreme_event +
      (1 | boxlabel),
    data
  )
  formula_resp <- drop_collapsed_formula_terms(
    y ~ swc + doy_c + doy_c2 + precipitation + robinia + culture + extreme_event +
      (1 | boxlabel) + (1 | tree_id),
    data
  )
  mod_swc <- suppressMessages(suppressWarnings(
    lme4::lmer(formula_swc, data = data, REML = TRUE)
  ))
  mod_resp <- suppressMessages(suppressWarnings(
    lme4::lmer(formula_resp, data = data, REML = TRUE)
  ))
  treatments <- c("precipitation", "robinia", "culture", "extreme_event")
  treatments <- treatments[vapply(treatments, function(treatment) {
    is.finite(coefficient_info(mod_swc, treatment)$est) &&
      is.finite(coefficient_info(mod_resp, treatment)$est)
  }, logical(1))]
  effects <- extract_effects(mod_swc, mod_resp, treatments)
  matrix_data <- make_matrix_data(effects)

  bundle <- list(
    data = data,
    mod_swc = mod_swc,
    mod_resp = mod_resp,
    effects = effects,
    effects_int = tibble::tibble(),
    matrix_data = matrix_data,
    modeled_factors = treatments,
    note = NA_character_,
    source_metadata = list(
      refreshed = Sys.time(),
      data_source = "current tracked interim data",
      formula_policy = "prespecified additive manuscript formula"
    )
  )
  saveRDS(bundle, source_rds_file, compress = "xz")
  readr::write_csv(matrix_data, matrix_file)
  readr::write_csv(data, sub("[.]rds$", "-data.csv", source_rds_file))

  tibble::tibble(
    species, resp_var, response_label, model_id, status = "complete",
    n_rows = nrow(data), source_rds_file, note = NA_character_
  )
}

inventory <- purrr::map_dfr(split(model_grid, seq_len(nrow(model_grid))), refresh_one)
readr::write_csv(inventory, file.path(source_root, "refreshed-source-inventory.csv"))
print(inventory)
cat("Refreshed repeated-response source models in:", source_root, "\n")
