# Growth data preparation used by Figure 6.

GROWTH_FACTOR_LEVELS <- list(
  precipitation = alinv_factor_levels("precipitation"),
  culture = alinv_factor_levels("culture"),
  robinia = alinv_factor_levels("robinia"),
  species = alinv_factor_levels("species")
)

growth_metric_catalog <- function() {
  tibble::tribble(
    ~resp_var, ~within_phase,
    "height", FALSE,
    "height_inc_t0", FALSE,
    "height_inc_t0_rel", FALSE,
    "height_inc_phase_abs", TRUE,
    "height_inc_phase_rel", TRUE,
    "diameter", FALSE,
    "diameter_inc_t0", FALSE,
    "diameter_inc_t0_rel", FALSE,
    "diameter_inc_phase_abs", TRUE,
    "diameter_inc_phase_rel", TRUE,
    "volume", FALSE,
    "volume_inc_t0", FALSE,
    "volume_inc_t0_rel", FALSE,
    "volume_inc_phase_abs", TRUE,
    "volume_inc_phase_rel", TRUE
  )
}

growth_metric_spec <- function(resp_var) {
  spec <- growth_metric_catalog() %>%
    dplyr::filter(.data$resp_var == .env$resp_var) %>%
    dplyr::slice_head(n = 1)
  if (!nrow(spec)) {
    stop("Unknown growth response variable: ", resp_var, call. = FALSE)
  }
  spec
}

normalize_growth_factors <- function(df) {
  df %>%
    dplyr::mutate(
      precipitation = factor(.data$precipitation, GROWTH_FACTOR_LEVELS$precipitation),
      culture = factor(.data$culture, GROWTH_FACTOR_LEVELS$culture),
      robinia = factor(.data$robinia, GROWTH_FACTOR_LEVELS$robinia),
      species = factor(.data$species, GROWTH_FACTOR_LEVELS$species)
    )
}

add_growth_phase_restart_rows <- function(df_metric, resp_var) {
  phase_levels <- levels(df_metric$phase)
  if (is.null(phase_levels) || !"phase" %in% names(df_metric)) {
    return(df_metric %>% dplyr::mutate(is_phase_restart = FALSE))
  }

  df_base <- df_metric %>% dplyr::mutate(is_phase_restart = FALSE)
  restart_rows <- df_base %>%
    dplyr::group_by(.data$tree_id) %>%
    dplyr::group_modify(function(.x, .y) {
      purrr::map_dfr(phase_levels[-1], function(phase_i) {
        current_rows <- .x %>%
          dplyr::filter(as.character(.data$phase) == phase_i) %>%
          dplyr::arrange(.data$date)
        previous_rows <- .x %>%
          dplyr::filter(as.integer(.data$phase) < match(phase_i, phase_levels)) %>%
          dplyr::arrange(.data$date)

        if (!nrow(current_rows) || !nrow(previous_rows)) {
          return(tibble::tibble())
        }

        restart_row <- current_rows[1, , drop = FALSE]
        restart_row$date <- max(previous_rows$date, na.rm = TRUE)
        restart_row[[resp_var]] <- 0
        restart_row$is_phase_restart <- TRUE
        restart_row
      })
    }) %>%
    dplyr::ungroup()

  dplyr::bind_rows(df_base, restart_rows) %>%
    dplyr::arrange(.data$tree_id, .data$date, .data$is_phase_restart)
}

prepare_growth_metric_plot_data <- function(resp_var,
                                            soil_type = NULL,
                                            include_soil_treatment = NULL) {
  context <- alinv_get_analysis_context()
  soil_type <- soil_type %||% context$soil_filter %||% "both"
  include_soil_treatment <- alinv_resolve_include_soil_treatment(
    include_soil_treatment = include_soil_treatment,
    soil_filter = soil_type
  )
  within_phase <- isTRUE(growth_metric_spec(resp_var)$within_phase[[1]])

  df_metric <- get_data("tree", "growth") %>%
    alinv_apply_soil_context(
      soil_filter = soil_type,
      include_soil_treatment = include_soil_treatment
    ) %>%
    dplyr::filter(!is.na(.data[[resp_var]])) %>%
    normalize_growth_factors() %>%
    dplyr::mutate(
      phase = factor(
        .data$phase,
        levels = c("until June", "July-August", "September+")
      )
    )

  if (!nrow(df_metric)) return(df_metric)
  if (within_phase) {
    add_growth_phase_restart_rows(df_metric, resp_var)
  } else {
    df_metric %>% dplyr::mutate(is_phase_restart = FALSE)
  }
}

summarize_growth_metric_plot_data <- function(df_metric,
                                              resp_var,
                                              within_phase = FALSE) {
  keys <- c("species", "robinia", "precipitation", "culture")
  if (isTRUE(within_phase)) keys <- c(keys, "phase")

  response <- rlang::sym(resp_var)
  df_metric %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(keys)), .data$date) %>%
    dplyr::summarise(
      mean = mean(!!response, na.rm = TRUE),
      sd = stats::sd(!!response, na.rm = TRUE),
      n = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      se = .data$sd / sqrt(pmax(.data$n, 1)),
      line_group = if (isTRUE(within_phase)) {
        interaction(.data$precipitation, .data$culture, .data$phase, drop = TRUE)
      } else {
        interaction(.data$precipitation, .data$culture, drop = TRUE)
      }
    )
}
