# Source shared functions robustly regardless of current working directory.
if (!exists("get_data_var_grid", mode = "function")) {
  resolve_source_path <- function() {
    candidates <- c(
      "functions/_source.R",
      "../functions/_source.R",
      "./functions/_source.R"
    )
    hit <- candidates[file.exists(candidates)][1]
    if (!is.na(hit)) return(hit)

    # Walk up a few levels to find project root containing functions/_source.R
    wd <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
    for (i in 0:5) {
      root_i <- wd
      if (i > 0) {
        for (j in seq_len(i)) root_i <- dirname(root_i)
      }
      cand <- file.path(root_i, "functions", "_source.R")
      if (file.exists(cand)) return(cand)
    }
    stop("Could not locate functions/_source.R from working directory: ", getwd())
  }

  source(resolve_source_path())
}

.make_failure_plot <- function(title, msg) {
  ggplot2::ggplot() +
    ggplot2::annotate(
      "text",
      x = 1,
      y = 1,
      label = paste0("Model failed:\n", msg),
      hjust = 0,
      vjust = 1,
      size = 4
    ) +
    ggplot2::xlim(0.8, 1.2) +
    ggplot2::ylim(0.8, 1.2) +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 11))
}

make_temporal_sem_combo <- function(
    type                = "tree",
    data_name,
    resp_var,
    species,
    soil_type           = "both",
    # SEM options
    include_interaction = TRUE,
    scale_all_numeric   = TRUE,
    do_rfe              = TRUE,
    aic_improve         = 2,
    # temporal GLMM options
    add_covars          = FALSE,
    covars_fun          = NULL,
    # SWC source
    swc_source          = "measured",
    # output
    outdir_root         = "./output",
    force_refit         = FALSE
) {
  # ------------------------------------------------------------
  # 0) Paths and cache file
  # ------------------------------------------------------------
  today  <- format(Sys.Date(), "%Y-%m-%d")
  outdir <- file.path(outdir_root, today, "temporal_factor+sem")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  base_name <- paste0(
    "temporal_sem_",
    type, "_", data_name,
    if (!is.null(resp_var) && nzchar(resp_var)) paste0("_", resp_var) else "",
    "_", species,
    "_soil-", soil_type,
    "_semInt-", ifelse(isTRUE(include_interaction), "yes", "no"),
    "_swc-", swc_source
  )
  
  png_path <- file.path(outdir, paste0(base_name, ".png"))
  rds_path <- file.path(outdir, paste0(base_name, ".rds"))

  # Search existing caches across dated output folders.
  cache_candidates <- list.files(
    path = outdir_root,
    pattern = "\\.rds$",
    recursive = TRUE,
    full.names = TRUE
  )
  cache_candidates <- cache_candidates[basename(cache_candidates) == paste0(base_name, ".rds")]
  cache_candidates <- cache_candidates[grepl("temporal_factor\\+sem", cache_candidates)]

  latest_cache <- NA_character_
  if (length(cache_candidates)) {
    mt <- file.info(cache_candidates)$mtime
    latest_cache <- cache_candidates[which.max(mt)]
  }
  
  # ------------------------------------------------------------
  # 1) Load cached models if available, otherwise run them
  # ------------------------------------------------------------
  if (!isTRUE(force_refit) && file.exists(rds_path)) {
    cache <- tryCatch(readRDS(rds_path), error = function(e) NULL)
    if (is.null(cache)) {
      warning("Failed to read cache file: ", rds_path, ". Re-fitting models.")
      res_temporal <- NULL
      res_sem      <- NULL
      cache_used   <- NA_character_
    } else {
      res_temporal <- cache$temporal
      res_sem      <- cache$sem
      cache_used   <- rds_path
    }
  } else if (!isTRUE(force_refit) && !is.na(latest_cache) && file.exists(latest_cache)) {
    cache <- tryCatch(readRDS(latest_cache), error = function(e) NULL)
    if (is.null(cache)) {
      warning("Failed to read cache file: ", latest_cache, ". Re-fitting models.")
      res_temporal <- NULL
      res_sem      <- NULL
      cache_used   <- NA_character_
    } else {
      res_temporal <- cache$temporal
      res_sem      <- cache$sem
      cache_used   <- latest_cache

      cache_age_days <- as.numeric(difftime(Sys.time(), file.info(latest_cache)$mtime, units = "days"))
      warning(
        sprintf(
          paste0(
            "Using cached temporal+SEM models from '%s' (%.1f days old). ",
            "To update models, rerun with: make_temporal_sem_combo(..., force_refit = TRUE) ",
            "or run_temporal_sem_over_grid(..., force_refit = TRUE)."
          ),
          latest_cache,
          cache_age_days
        ),
        call. = FALSE
      )
    }
  } else {
    res_temporal <- NULL
    res_sem      <- NULL
    cache_used   <- NA_character_
  }

  # ------------------------------------------------------------
  # 1b) Fit models when no usable cache was loaded
  # ------------------------------------------------------------
  if (is.null(res_temporal) || is.null(res_sem)) {
    temporal_err <- NULL
    sem_err <- NULL

    # Time-varying treatment effects (left plot)
    res_temporal <- tryCatch(
      make_effect_figure_generic(
        type           = type,
        data_name      = data_name,
        resp_var       = resp_var,
        target_species = species,
        soil_type      = soil_type,
        add_covars     = add_covars,
        covars_fun     = covars_fun,
        swc_source     = swc_source
      ),
      error = function(e) {
        temporal_err <<- conditionMessage(e)
        list(
          plot = .make_failure_plot("Temporal GLMM", temporal_err),
          error = temporal_err,
          effects = tibble::tibble()
        )
      }
    )

    # SEM (right plot)
    res_sem <- tryCatch(
      run_sem_for_trait(
        type                = type,
        data_name           = data_name,
        resp_var            = resp_var,
        species             = species,
        soil_type           = soil_type,
        include_interaction = include_interaction,
        scale_all_numeric   = scale_all_numeric,
        do_rfe              = do_rfe,
        aic_improve         = aic_improve,
        swc_source          = swc_source
      ),
      error = function(e) {
        sem_err <<- conditionMessage(e)
        list(
          plot = .make_failure_plot("SEM", sem_err),
          error = sem_err
        )
      }
    )

    saveRDS(
      list(temporal = res_temporal, sem = res_sem),
      file = rds_path
    )
    cache_used <- rds_path

    if (!is.null(temporal_err) || !is.null(sem_err)) {
      warn_parts <- c()
      if (!is.null(temporal_err)) warn_parts <- c(warn_parts, paste0("Temporal GLMM failed: ", temporal_err))
      if (!is.null(sem_err)) warn_parts <- c(warn_parts, paste0("SEM failed: ", sem_err))
      warning(paste(warn_parts, collapse = " | "), call. = FALSE)
    }
  }
  
  # ------------------------------------------------------------
  # 2) Extract plots and sanity checks
  # ------------------------------------------------------------
  p_temporal <- res_temporal$plot
  if (is.null(p_temporal)) {
    msg <- if (!is.null(res_temporal$error)) res_temporal$error else "No contrasts available"
    warning("Temporal effects plot is NULL - ", msg)
    p_temporal <- .make_failure_plot("Temporal GLMM", msg)
  }
  
  p_sem <- res_sem$plot
  if (is.null(p_sem)) {
    msg <- if (!is.null(res_sem$error)) res_sem$error else "No SEM output available"
    warning("SEM plot is NULL - ", msg)
    p_sem <- .make_failure_plot("SEM", msg)
  }
  
  # ------------------------------------------------------------
  # 3) Tidy individual plots before combining
  # ------------------------------------------------------------
  p_temporal <- p_temporal +
    ggplot2::theme(
      plot.margin = ggplot2::margin(t = 5.5, r = 15, b = 5.5, l = 5.5),
      plot.title  = ggplot2::element_text(size = 11)
    )
  
  # drop SEM title (we'll add a joint title) and add left margin
  p_sem <- p_sem +
    ggplot2::theme(
      plot.margin = ggplot2::margin(t = 5.5, r = 5.5, b = 5.5, l = 15),
      plot.title  = ggplot2::element_blank()
    )
  
  # Combined title
  joint_title <- paste0(
    "Temporal GLMM and SEM for ", species,
    " (soil: ", soil_type, ", data: ", data_name,
    if (!is.null(resp_var) && nzchar(resp_var)) paste0(", ", resp_var) else "", ")"
  )
  
  # ------------------------------------------------------------
  # 4) Combine side by side (patchwork)
  # ------------------------------------------------------------
  combined_plot <- patchwork::wrap_plots(
    p_temporal,
    p_sem,
    ncol   = 2,
    widths = c(1.6, 1)
  ) +
    patchwork::plot_annotation(
      title = joint_title
    ) &
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 12)
    )
  
  # ------------------------------------------------------------
  # 5) Save PNG
  # ------------------------------------------------------------
  ggplot2::ggsave(
    filename = png_path,
    plot     = combined_plot,
    width    = 15,
    height   = 10,
    dpi      = 300
  )
  
  # ------------------------------------------------------------
  # 6) Return everything useful
  # ------------------------------------------------------------
  list(
    temporal  = res_temporal,
    sem       = res_sem,
    combined  = combined_plot,
    file      = png_path,
    cache_rds = cache_used
  )
}

# ============================================================
# EXAMPLE CODE: Uncomment and run interactively when needed
# ============================================================
# combo_res <- make_temporal_sem_combo(
#   type      = "tree",
#   data_name = "growth",
#   resp_var  = "height",
#   species   = "fagus",
#   soil_type = "both"
# )
#
# combo_res$file   # path to PNG
# combo_res$combined


#' Run temporal + SEM combos over grid (optionally in parallel)
#'
#' @param species_vec character vector of species
#' @param soil_type   soil type passed to make_temporal_sem_combo()
#' @param parallel    logical; use furrr::future_pmap if TRUE
#' @param workers     number of parallel workers (if parallel = TRUE)
run_temporal_sem_over_grid <- function(
    species_vec = c("fagus", "quercus"),
    soil_type   = "both",
    swc_source  = "measured",
    parallel    = FALSE,
    workers     = max(1L, parallel::detectCores() - 1L),
    ...
) {
  var_grid <- get_data_var_grid() %>%
    dplyr::filter(.data$type == "tree")  # box vars usually don't go into SEM
  
  combos <- tidyr::crossing(
    var_grid,
    tibble::tibble(species = species_vec)
  )
  # combos has columns: type, data, resp_var, species
  
  run_one <- function(type, data, resp_var, species) {
    msg <- paste0("Running combo: type=", type,
                  ", data=", data,
                  ", resp_var=", resp_var,
                  ", species=", species)
    message(msg)
    
    tryCatch(
      make_temporal_sem_combo(
        type      = type,
        data_name = data,
        resp_var  = resp_var,
        species   = species,
        soil_type = soil_type,
        swc_source = swc_source,
        ...
      ),
      error = function(e) {
        warning("Failed: ", msg, " -> ", conditionMessage(e))
        NULL
      }
    )
  }
  
  if (!parallel) {
    # sequential version
    results <- purrr::pmap(
      combos,
      run_one
    )
  } else {
    # parallel version using furrr
    if (!requireNamespace("furrr", quietly = TRUE) ||
        !requireNamespace("future", quietly = TRUE)) {
      stop("Packages 'furrr' and 'future' are required for parallel = TRUE.")
    }
    
    future::plan(future::multisession, workers = workers)
    on.exit(future::plan(future::sequential), add = TRUE)
    
    results <- furrr::future_pmap(
      combos,
      run_one
    )
  }
  
  # return a tibble tying metadata to result objects
  dplyr::bind_cols(combos, tibble::tibble(result = results))
}

# run_temporal_sem_over_grid(
#     species_vec = c("fagus", "quercus"),
#     soil_type   = "both",
#     parallel    = FALSE,
#     workers     = max(1L, parallel::detectCores() - 1L),
# )
