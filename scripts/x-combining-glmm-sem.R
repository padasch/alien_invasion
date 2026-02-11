# library(patchwork)  # at top of your script

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
    # output
    outdir_root         = "./output"
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
    "_soil-", soil_type
  )
  
  png_path <- file.path(outdir, paste0(base_name, ".png"))
  rds_path <- file.path(outdir, paste0(base_name, ".rds"))
  
  # ------------------------------------------------------------
  # 1) Load cached models if available, otherwise run them
  # ------------------------------------------------------------
  if (file.exists(rds_path)) {
    cache        <- readRDS(rds_path)
    res_temporal <- cache$temporal
    res_sem      <- cache$sem
  } else {
    # Time-varying treatment effects (left plot)
    res_temporal <- make_effect_figure_generic(
      type           = type,
      data_name      = data_name,
      resp_var       = resp_var,
      target_species = species,
      soil_type      = soil_type,
      add_covars     = add_covars,
      covars_fun     = covars_fun
    )
    
    # SEM (right plot)
    res_sem <- run_sem_for_trait(
      type                = type,
      data_name           = data_name,
      resp_var            = resp_var,
      species             = species,
      soil_type           = soil_type,
      include_interaction = include_interaction,
      scale_all_numeric   = scale_all_numeric,
      do_rfe              = do_rfe,
      aic_improve         = aic_improve
    )
    
    saveRDS(
      list(temporal = res_temporal, sem = res_sem),
      file = rds_path
    )
  }
  
  # ------------------------------------------------------------
  # 2) Extract plots and sanity checks
  # ------------------------------------------------------------
  p_temporal <- res_temporal$plot
  if (is.null(p_temporal)) {
    warning("Temporal effects plot is NULL - no contrasts available?")
  }
  
  p_sem <- res_sem$plot
  if (is.null(p_sem)) {
    warning("SEM plot is NULL - check run_sem_for_trait output.")
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
    cache_rds = rds_path
  )
}

combo_res <- make_temporal_sem_combo(
  type      = "tree",
  data_name = "growth",
  resp_var  = "height",
  species   = "fagus",
  soil_type = "both"
)

combo_res$file   # path to PNG
combo_res$combined


#' Run temporal + SEM combos over grid (optionally in parallel)
#'
#' @param species_vec character vector of species
#' @param soil_type   soil type passed to make_temporal_sem_combo()
#' @param parallel    logical; use furrr::future_pmap if TRUE
#' @param workers     number of parallel workers (if parallel = TRUE)
run_temporal_sem_over_grid <- function(
    species_vec = c("fagus", "quercus"),
    soil_type   = "both",
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

run_temporal_sem_over_grid(
    species_vec = c("fagus", "quercus"),
    soil_type   = "both",
    parallel    = FALSE,
    workers     = max(1L, parallel::detectCores() - 1L),
)
