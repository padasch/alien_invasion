# Shared accessors for the design-aligned bootstrap analyses used by the
# publication figures. Model fitting remains in the family-specific scripts;
# final and supplementary plotting code should use these functions rather than
# redefine resampling, interval, or significance logic.

ALINV_BOOTSTRAP_TARGET <- 1000L

alinv_bootstrap_analysis_root <- function(project_root = .alinv_project_root()) {
  file.path(project_root, "exploration", "2026-08-18 Testing bootstrapping")
}

alinv_bootstrap_paths <- function(project_root = .alinv_project_root()) {
  root <- alinv_bootstrap_analysis_root(project_root)
  list(
    root = root,
    temporal_script = file.path(root, "temporal_lmm", "scripts", "01_temporal_lmm_bootstrap.R"),
    temporal_effects = file.path(root, "temporal_lmm", "output", "data", "temporal-lmm-bootstrap-effects.csv"),
    temporal_status = file.path(root, "temporal_lmm", "output", "status", "temporal-lmm-bootstrap-status.csv"),
    biomass_script = file.path(root, "biomass_lmm", "scripts", "01-biomass-container-bootstrap.R"),
    biomass_effects = file.path(root, "biomass_lmm", "output", "biomass-wald-vs-bootstrap-comparison.csv"),
    biomass_status = file.path(root, "biomass_lmm", "output", "biomass-bootstrap-status.csv"),
    repeated_sem_script = file.path(root, "repeated_response_sem", "scripts", "01_repeated_response_sem_bootstrap.R"),
    repeated_sem_totals = file.path(root, "repeated_response_sem", "output", "figure6-ready-repeated-response-sem-path-summed-total.csv"),
    repeated_sem_effects = file.path(root, "repeated_response_sem", "output", "repeated-response-sem-bootstrap-effects.csv"),
    repeated_sem_status = file.path(root, "repeated_response_sem", "output", "repeated-response-sem-bootstrap-status.csv"),
    reduced_response_script = file.path(root, "reduced_response_lmm", "scripts", "01_reduced_response_bootstrap.R"),
    reduced_response_effects = file.path(root, "reduced_response_lmm", "output", "reduced-response-bootstrap-effects.csv"),
    reduced_response_status = file.path(root, "reduced_response_lmm", "output", "reduced-response-bootstrap-status.csv"),
    phenology_script = file.path(root, "phenology", "run_phenology_bootstrap.R"),
    phenology_totals = file.path(root, "phenology", "output", "figure6-ready-phenology-path-summed-total.csv"),
    phenology_primary = file.path(root, "phenology", "output", "primary-common-shift-bootstrap-effects.csv"),
    phenology_timing_index = file.path(root, "phenology", "output", "stage-centred-timing-index.csv"),
    phenology_sem_effects = file.path(root, "phenology", "output", "sem-effect-decomposition-block-stratified.csv"),
    phenology_sem_paths = file.path(root, "phenology", "output", "sem-constituent-paths-block-stratified.csv"),
    phenology_status = file.path(root, "phenology", "output", "phenology-bootstrap-model-status.csv")
  )
}

alinv_bootstrap_p <- function(draws) {
  draws <- draws[is.finite(draws)]
  if (!length(draws)) return(NA_real_)
  min(
    1,
    2 * min(
      (sum(draws <= 0) + 1) / (length(draws) + 1),
      (sum(draws >= 0) + 1) / (length(draws) + 1)
    )
  )
}

alinv_bootstrap_percentile <- function(draws, probs = c(0.025, 0.975)) {
  draws <- draws[is.finite(draws)]
  if (!length(draws)) return(rep(NA_real_, length(probs)))
  unname(stats::quantile(draws, probs = probs, type = 7, names = FALSE))
}

alinv_bootstrap_read_csv <- function(path, required_cols = character()) {
  if (!file.exists(path)) return(NULL)
  x <- tryCatch(
    readr::read_csv(path, show_col_types = FALSE),
    error = function(e) NULL
  )
  if (is.null(x) || !all(required_cols %in% names(x))) return(NULL)
  x
}

alinv_bootstrap_valid <- function(family, project_root = .alinv_project_root(), target = ALINV_BOOTSTRAP_TARGET) {
  paths <- alinv_bootstrap_paths(project_root)
  family <- match.arg(family, c("temporal", "biomass", "repeated_sem", "reduced_response", "phenology"))

  if (family == "temporal") {
    effects <- alinv_bootstrap_read_csv(paths$temporal_effects, c("model_id", "bootstrap_replicates", "lower", "upper", "p_boot"))
    status <- alinv_bootstrap_read_csv(paths$temporal_status, c("model_id", "successful_replicates"))
    return(!is.null(effects) && !is.null(status) && nrow(status) == 4L &&
      all(effects$bootstrap_replicates >= target) && all(status$successful_replicates >= target))
  }
  if (family == "biomass") {
    effects <- alinv_bootstrap_read_csv(paths$biomass_effects, c("species", "metric", "term", "ci_lo_boot", "ci_hi_boot", "p_boot", "n_boot"))
    status <- alinv_bootstrap_read_csv(paths$biomass_status, c("species", "metric", "successful"))
    return(!is.null(effects) && !is.null(status) && nrow(effects) == 18L && nrow(status) == 6L &&
      all(effects$n_boot >= target) && all(status$successful >= target))
  }
  if (family == "repeated_sem") {
    effects <- alinv_bootstrap_read_csv(paths$repeated_sem_effects, c("species", "resp_var", "treatment", "component", "lower", "upper", "p_boot", "n_boot"))
    status <- alinv_bootstrap_read_csv(paths$repeated_sem_status, c("species", "resp_var", "status", "n_boot_success"))
    estimable <- if (is.null(status)) logical() else status$status == "complete"
    return(!is.null(effects) && !is.null(status) &&
      all(effects$n_boot >= target) && all(status$n_boot_success[estimable] >= target))
  }
  if (family == "reduced_response") {
    effects <- alinv_bootstrap_read_csv(
      paths$reduced_response_effects,
      c("species", "resp_var", "treatment", "estimate", "lower", "upper", "p_boot", "n_boot")
    )
    status <- alinv_bootstrap_read_csv(
      paths$reduced_response_status,
      c("species", "resp_var", "status", "n_boot_success")
    )
    estimable <- if (is.null(status)) logical() else status$status == "complete"
    return(!is.null(effects) && !is.null(status) &&
      all(effects$n_boot >= target) && all(status$n_boot_success[estimable] >= target))
  }

  totals <- alinv_bootstrap_read_csv(paths$phenology_totals, c("species", "effect", "lower_95", "upper_95", "p_boot", "n_boot_success"))
  primary <- alinv_bootstrap_read_csv(paths$phenology_primary, c("species", "effect", "lower_95_oriented", "upper_95_oriented", "p_boot", "n_boot_success"))
  timing <- alinv_bootstrap_read_csv(paths$phenology_timing_index, c("species", "timing_index_days"))
  !is.null(totals) && !is.null(primary) && !is.null(timing) &&
    nrow(totals) == 6L && nrow(primary) == 6L &&
    all(totals$n_boot_success >= target) && all(primary$n_boot_success >= target)
}

alinv_run_bootstrap_family <- function(family,
                                       project_root = .alinv_project_root(),
                                       target = ALINV_BOOTSTRAP_TARGET,
                                       cores = 4L) {
  paths <- alinv_bootstrap_paths(project_root)
  family <- match.arg(family, c("temporal", "biomass", "repeated_sem", "reduced_response", "phenology"))
  script <- switch(
    family,
    temporal = paths$temporal_script,
    biomass = paths$biomass_script,
    repeated_sem = paths$repeated_sem_script,
    reduced_response = paths$reduced_response_script,
    phenology = paths$phenology_script
  )
  if (!file.exists(script)) stop("Missing bootstrap analysis script: ", script, call. = FALSE)

  args <- c("--vanilla", shQuote(script))
  env <- character()
  if (family == "temporal") {
    env <- c(
      paste0("ALINV_TEMPORAL_BOOT_B=", as.integer(target)),
      paste0("ALINV_TEMPORAL_BOOT_CORES=", as.integer(cores))
    )
  } else if (family %in% c("repeated_sem", "reduced_response", "phenology")) {
    args <- c(args, paste0("--bootstrap=", as.integer(target)), paste0("--cores=", as.integer(cores)))
  } else if (target != 1000L) {
    stop("The promoted biomass script currently supports 1,000 replicates.", call. = FALSE)
  }

  status <- system2(file.path(R.home("bin"), "Rscript"), args = args, env = env)
  if (!identical(status, 0L) || !alinv_bootstrap_valid(family, project_root, target)) {
    stop("Bootstrap analysis failed validation for family: ", family, call. = FALSE)
  }
  invisible(TRUE)
}

alinv_ensure_bootstrap_family <- function(family,
                                          project_root = .alinv_project_root(),
                                          target = ALINV_BOOTSTRAP_TARGET,
                                          cores = 4L) {
  if (!alinv_bootstrap_valid(family, project_root, target)) {
    alinv_run_bootstrap_family(family, project_root, target, cores)
  }
  invisible(alinv_bootstrap_paths(project_root))
}

alinv_read_temporal_bootstrap_effects <- function(data_name, resp_var, species,
                                                   project_root = .alinv_project_root()) {
  paths <- alinv_ensure_bootstrap_family("temporal", project_root)
  readr::read_csv(paths$temporal_effects, show_col_types = FALSE) |>
    dplyr::filter(
      .data$data_name == .env$data_name,
      .data$response_var == .env$resp_var,
      .data$species == .env$species
    ) |>
    dplyr::mutate(
      date = as.Date(.data$date),
      source_file = paths$temporal_effects,
      uncertainty_method = "block-stratified container-cluster bootstrap percentile"
    )
}

alinv_read_biomass_bootstrap_effects <- function(project_root = .alinv_project_root()) {
  paths <- alinv_ensure_bootstrap_family("biomass", project_root)
  readr::read_csv(paths$biomass_effects, show_col_types = FALSE) |>
    dplyr::mutate(
      source_file = paths$biomass_effects,
      uncertainty_method = "block-stratified container-cluster bootstrap percentile"
    )
}

alinv_read_repeated_sem_bootstrap_effects <- function(project_root = .alinv_project_root()) {
  paths <- alinv_ensure_bootstrap_family("repeated_sem", project_root)
  readr::read_csv(paths$repeated_sem_effects, show_col_types = FALSE) |>
    dplyr::mutate(
      source_file = paths$repeated_sem_effects,
      uncertainty_method = "block-stratified container-cluster bootstrap percentile"
    )
}

alinv_read_repeated_sem_bootstrap_totals <- function(project_root = .alinv_project_root()) {
  alinv_read_repeated_sem_bootstrap_effects(project_root) |>
    dplyr::filter(.data$component == "total") |>
    dplyr::transmute(
      species = .data$species,
      treatment = .data$treatment,
      response_label = .data$response_label,
      resp_var = .data$resp_var,
      estimate = .data$estimate,
      lower = .data$lower,
      upper = .data$upper,
      p_value = .data$p_boot,
      n_boot_success = .data$n_boot,
      source_file = .data$source_file,
      uncertainty_method = .data$uncertainty_method
    )
}

alinv_read_reduced_response_bootstrap_effects <- function(project_root = .alinv_project_root()) {
  paths <- alinv_ensure_bootstrap_family("reduced_response", project_root)
  readr::read_csv(paths$reduced_response_effects, show_col_types = FALSE) |>
    dplyr::mutate(
      source_file = paths$reduced_response_effects,
      uncertainty_method = "block-stratified container-cluster bootstrap percentile"
    )
}

alinv_factor_main_effect_coef <- function(model, factor_name) {
  hits <- grep(factor_name, names(lme4::fixef(model)), fixed = TRUE, value = TRUE)
  hits <- hits[!grepl(":", hits, fixed = TRUE)]
  if (length(hits) != 1L) return(NA_character_)
  hits[[1]]
}

alinv_resample_containers_within_block <- function(data, replicate_id,
                                                    container = "boxlabel",
                                                    tree = "tree_id") {
  dat <- data |>
    dplyr::mutate(
      .container_original = as.character(.data[[container]]),
      .tree_original = as.character(.data[[tree]]),
      .block = sub("-.*$", "", .data$.container_original)
    )
  containers <- dat |>
    dplyr::distinct(.data$.block, .data$.container_original) |>
    dplyr::arrange(.data$.block, .data$.container_original)
  selected <- containers |>
    dplyr::group_by(.data$.block) |>
    dplyr::group_modify(~tibble::tibble(
      .container_original = sample(
        .x$.container_original,
        size = nrow(.x),
        replace = TRUE
      ),
      .draw = seq_len(nrow(.x))
    )) |>
    dplyr::ungroup()

  pieces <- purrr::pmap(selected, function(.block, .container_original, .draw) {
    new_container <- paste0(
      "boot", replicate_id, "_", .block, "_", sprintf("%02d", .draw)
    )
    dat |>
      dplyr::filter(.data$.container_original == .env$.container_original) |>
      dplyr::mutate(
        "{container}" := new_container,
        "{tree}" := paste(new_container, .data$.tree_original, sep = "__")
      )
  })

  dplyr::bind_rows(pieces) |>
    dplyr::select(-.data$.container_original, -.data$.tree_original, -.data$.block) |>
    dplyr::mutate(
      "{container}" := factor(.data[[container]]),
      "{tree}" := factor(.data[[tree]])
    )
}

alinv_read_phenology_bootstrap_totals <- function(project_root = .alinv_project_root()) {
  paths <- alinv_ensure_bootstrap_family("phenology", project_root)
  readr::read_csv(paths$phenology_totals, show_col_types = FALSE) |>
    dplyr::mutate(
      source_file = paths$phenology_totals,
      uncertainty_method = "block-stratified container-cluster bootstrap percentile"
    )
}

alinv_read_phenology_bootstrap_primary <- function(project_root = .alinv_project_root()) {
  paths <- alinv_ensure_bootstrap_family("phenology", project_root)
  readr::read_csv(paths$phenology_primary, show_col_types = FALSE) |>
    dplyr::mutate(
      source_file = paths$phenology_primary,
      uncertainty_method = "block-stratified container-cluster bootstrap percentile"
    )
}

alinv_read_phenology_bootstrap_primary_standardized <- function(project_root = .alinv_project_root()) {
  paths <- alinv_ensure_bootstrap_family("phenology", project_root)
  scales <- readr::read_csv(paths$phenology_timing_index, show_col_types = FALSE) |>
    dplyr::group_by(.data$species) |>
    dplyr::summarise(
      timing_index_sd_days = stats::sd(.data$timing_index_days, na.rm = TRUE),
      .groups = "drop"
    )
  alinv_read_phenology_bootstrap_primary(project_root) |>
    dplyr::left_join(scales, by = "species") |>
    dplyr::mutate(
      estimate = .data$estimate_oriented / .data$timing_index_sd_days,
      lower = .data$lower_95_oriented / .data$timing_index_sd_days,
      upper = .data$upper_95_oriented / .data$timing_index_sd_days,
      response_scale = "species-specific SD of the stage-centred tree timing index"
    )
}

alinv_read_phenology_sem_bootstrap_effects <- function(project_root = .alinv_project_root()) {
  paths <- alinv_ensure_bootstrap_family("phenology", project_root)
  effects <- readr::read_csv(paths$phenology_sem_effects, show_col_types = FALSE)
  constituent <- readr::read_csv(paths$phenology_sem_paths, show_col_types = FALSE)
  list(
    effects = dplyr::mutate(effects, source_file = paths$phenology_sem_effects),
    constituent = dplyr::mutate(constituent, source_file = paths$phenology_sem_paths)
  )
}
