# source("./functions/_source.R")

# !!! ####

# ================================================================
# Generalized per-date GLMM workflow for multiple data types
# ================================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(lme4)
  library(lubridate)
})
if (requireNamespace("emmeans", quietly = TRUE)) {
  emmeans::emm_options(lmer.df = "asymptotic")  # quiet, fast CIs
}

# ---------------------- CONFIG / MAPPING ------------------------

# Map data_name -> default response columns
.default_resp_map <- list(
  chlorophyll   = "chl",
  condition     = "condition",
  growth        = c(
    "height", "diameter", "volume",
    "height_inc_t0", "diameter_inc_t0", "volume_inc_t0",
    "height_inc_t0_rel", "diameter_inc_t0_rel", "volume_inc_t0_rel",
    "height_inc_phase_abs", "diameter_inc_phase_abs", "volume_inc_phase_abs",
    "height_inc_phase_rel", "diameter_inc_phase_rel", "volume_inc_phase_rel"
  ),
  quantum_yield = "qy",
  senescence    = c("percent_senesced", "chlavg"),
  phenology     = "doy"                          # special handling below
)
# ---------------------- HELPERS --------------------------------

# Compute increment since first measurement - ONLY for growth
.make_increment_since_first <- function(df, resp_var) {
  stopifnot(resp_var %in% names(df))
  df %>%
    arrange(tree_id, date) %>%
    group_by(tree_id) %>%
    mutate(
      resp0 = first(.data[[resp_var]]),
      y = .data[[resp_var]] - resp0
    ) %>%
    ungroup()
}

# Standardize factor baselines and drop NAs on required cols
.standardize_and_clean <- function(df, cols_needed) {
  df %>%
    dplyr::select(all_of(cols_needed)) %>%
    drop_na(date, robinia, precipitation, culture, species, tree_id, boxlabel, y) %>%
    mutate(
      date = as.Date(date),
      date = factor(date),
      robinia = factor(robinia, levels = c("without-robinia", "with-robinia")),
      precipitation = factor(precipitation, levels = c("control", "drought")),
      culture = factor(culture, levels = c("mono", "mixed")),
      soiltype = alinv_relevel_soiltype(soiltype),
      extreme_event = factor(extreme_event, levels = c("no", "yes")),
      species = factor(species)
    ) %>%
    droplevels()
}

# Optional covariate joiner (e.g., SWC) - supply a function returning a table
.maybe_add_covars <- function(df, add_covars, covars_fun) {
  if (!isTRUE(add_covars) || is.null(covars_fun)) return(df)
  cov_df <- covars_fun()
  by_cols <- intersect(c("boxlabel", "date"), names(df))
  by_cols <- union(by_cols, intersect(c("boxlabel", "date"), names(cov_df)))
  if (!length(by_cols)) stop("No common keys to join covariates.")
  df %>%
    left_join(cov_df, by = by_cols) %>%
    mutate(swc_sc = if ("swc" %in% names(.)) as.numeric(scale(swc)) else NA_real_)
}

# NEW: Pivot phenology (stage1..4, doy_s1..4) to long, bin dates to weeks,
# and return a frame with a single response column 'y' (DOY) and a factor 'date'
# NEW — robust phenology wrangling

# ---------------------- DATA PREP (GENERIC) --------------------

# Pull data, pick/validate response, do growth increment if needed,
# standardize, and return a modeling-ready frame with column 'y'
prepare_df_generic <- function(
    type = "tree",
    data_name = c("chlorophyll", "condition", "growth", "quantum_yield", "senescence", "phenology"),
    resp_var = NULL,                 # override mapping if you want a specific column
    species_keep = NULL,             # e.g., c("fagus","quercus") or "fagus"
    standardize_response = TRUE,
    add_covars = FALSE,
    covars_fun = NULL,                # function returning covariates (boxlabel+date)
    soil_type = "both",
    include_soil_treatment = NULL,
    swc_source = "measured"          # "measured" or "imputed_gam"
) {
  data_name <- match.arg(data_name)
  include_soil_treatment <- alinv_resolve_include_soil_treatment(
    include_soil_treatment = include_soil_treatment,
    soil_filter = soil_type
  )
  
  # Pull 
  df_raw <- get_data(type = type, data_name = data_name, swc_source = swc_source)
  
  # Early species filter if present
  if (!is.null(species_keep) && "species" %in% names(df_raw)) {
    df_raw <- dplyr::filter(df_raw, species %in% species_keep)
  }
  
  df_raw <- alinv_apply_soil_context(
    df_raw,
    soil_filter = soil_type,
    include_soil_treatment = include_soil_treatment
  )
  
  default_resp <- .default_resp_map[[data_name]]
  if (is.null(resp_var)) {
    if (length(default_resp) > 1) {
      stop(
        "Please set resp_var explicitly. Options for data_name='", data_name,
        "' are: ", paste(default_resp, collapse = ", ")
      )
    }
    resp_var <- default_resp
  }
  if (!resp_var %in% names(df_raw)) {
    stop("Requested resp_var '", resp_var, "' not found in dataset for ", data_name, ".")
  }
  
  base_cols <- c("tree_id", "boxlabel", "date", "date_num", "species", "culture", "soiltype", "extreme_event",
                 "precipitation", "robinia", "swc", resp_var)
  df <- dplyr::select(df_raw, dplyr::all_of(base_cols))
  
  # Response column 'y' 
  df <- dplyr::rename(df, y = !!resp_var)
  
  # Optional covariates (e.g., SWC)
  df <- .maybe_add_covars(df, add_covars, covars_fun)
  
  # Standardize factors/baselines and drop NAs on required columns
  cols_needed <- c("y", "date", "date_num", "robinia", "precipitation", "culture", "soiltype", "extreme_event",
                   "species", "tree_id", "boxlabel", "swc")
  cols_needed <- intersect(cols_needed, names(df))
  
  df <- .standardize_and_clean(df, cols_needed)
  df <- df |> drop_na(y)

  if (isTRUE(standardize_response)) {
    df <- df %>%
      mutate(
        y_org = y,
        y = as.numeric(scale(y_org))
      )
  }

  # Keep a standardized SWC covariate available for temporal GLMMs.
  # This allows add_covars=TRUE to include SWC regardless of an external covariate join.
  if ("swc" %in% names(df)) {
    df <- df %>% mutate(swc_sc = as.numeric(scale(swc)))
  }
  
  df
}

# ---------------------- MODEL (UNCHANGED) ----------------------

# Per-date GLMM with date main effect and interactions (same as before)
fit_glmm_per_date <- function(df,
                              include_covars = FALSE,
                              include_soil_treatment = NULL) {
  include_soil_treatment <- alinv_resolve_include_soil_treatment(
    include_soil_treatment = include_soil_treatment,
    soil_filter = if ("soiltype" %in% names(df) && dplyr::n_distinct(df$soiltype, na.rm = TRUE) > 1) "both" else "single"
  )

  # start with all candidate terms
  rhs_terms <- c("date", "date:robinia", "date:precipitation", "date:culture")
  if (isTRUE(include_soil_treatment)) {
    rhs_terms <- c(rhs_terms, "date:soiltype")
  }
  
  # mapping from interaction terms to the underlying factor to check
  term_factor_map <- c(
    "date:robinia"      = "robinia",
    "date:precipitation" = "precipitation",
    "date:culture"      = "culture",
    "date:soiltype"     = "soiltype",
    "extreme_event" = "extreme_event"
  )
  
  
  # check each factor and remove interaction term if < 2 levels or missing
  for (term in names(term_factor_map)) {
    var <- term_factor_map[[term]]
    
    if (!var %in% names(df)) {
      message("Variable '", var, "' not found in data, removing term '", term, "'.")
      rhs_terms <- setdiff(rhs_terms, term)
      next
    }
    
    n_levels <- length(unique(stats::na.omit(df[[var]])))
    
    if (n_levels < 2) {
      message(
        "Variable '", var, "' has only ", n_levels,
        " level in this subset, removing term '", term, "'."
      )
      rhs_terms <- setdiff(rhs_terms, term)
    }
  }
  
  # collapse remaining fixed-effect terms
  rhs <- paste(rhs_terms, collapse = " + ")
  
  # optional covariate
  if (isTRUE(include_covars) && "swc_sc" %in% names(df)) {
    rhs <- paste(rhs, "+ swc_sc")
  }
  
  fml <- as.formula(paste("y ~", rhs, "+ (1|boxlabel) + (1|tree_id)"))
  lme4::lmer(fml, data = df, REML = TRUE)
}

extract_per_date_contrasts_from_fixef <- function(fit, factor_name) {
  mf <- stats::model.frame(fit)
  if (!"date" %in% names(mf) || !factor_name %in% names(mf)) {
    stop("Model frame does not contain 'date' and factor '", factor_name, "'.", call. = FALSE)
  }

  factor_levels <- levels(mf[[factor_name]])
  if (length(factor_levels) < 2L) {
    return(tibble())
  }

  treatment_level <- factor_levels[[2]]
  date_levels <- levels(mf$date)
  beta <- lme4::fixef(fit)
  vc <- as.matrix(stats::vcov(fit))
  coef_names <- paste0("date", date_levels, ":", factor_name, treatment_level)
  keep <- coef_names %in% names(beta)

  if (!any(keep)) {
    stop(
      "Could not find date-specific coefficients for factor '", factor_name,
      "' in the fitted model.",
      call. = FALSE
    )
  }

  coef_names <- coef_names[keep]
  date_levels <- date_levels[keep]
  raw_est <- unname(beta[coef_names])
  se <- vapply(
    coef_names,
    function(term_i) sqrt(unname(vc[term_i, term_i])),
    numeric(1)
  )

  estimate <- raw_est
  stat <- estimate / se

  tibble::tibble(
    date = as.Date(date_levels),
    contrast = paste0(treatment_level, " - ", factor_levels[[1]]),
    estimate = estimate,
    se = se,
    stat = stat,
    df = NA_real_,
    lower = estimate - 1.96 * se,
    upper = estimate + 1.96 * se
  ) %>%
    dplyr::arrange(.data$date, .data$contrast)
}

swap_contrast_label <- function(x) {
  parts <- strsplit(x, " - ", fixed = TRUE)[[1]]
  if (length(parts) != 2L) {
    return(x)
  }
  paste(parts[[2]], parts[[1]], sep = " - ")
}

# Robust per-date contrasts for ANY categorical factor; convert emmeans
# revpairwise output to treatment − baseline for consistent display.
extract_per_date_contrasts_any <- function(fit, factor_name) {
  if (!requireNamespace("emmeans", quietly = TRUE)) {
    message(
      "Package 'emmeans' not available; extracting temporal contrasts from model coefficients for factor '",
      factor_name,
      "'."
    )
    return(extract_per_date_contrasts_from_fixef(fit, factor_name))
  }

  emm <- emmeans::emmeans(fit, specs = as.formula(paste("~", factor_name, "| date")), type = "link")
  ct <- emmeans::contrast(emm, method = "revpairwise", by = "date", adjust = "none")
  sm <- suppressMessages(suppressWarnings(summary(ct, infer = TRUE)))
  out <- tibble::as_tibble(sm)
  
  has_col <- function(nm) nm %in% names(out)
  zcol  <- if (has_col("z.ratio")) "z.ratio" else if (has_col("t.ratio")) "t.ratio" else NA_character_
  secol <- if (has_col("SE")) "SE" else "SE"
  lowc  <- if (has_col("lower.CL")) "lower.CL" else if (has_col("asymp.LCL")) "asymp.LCL" else "lower.CL"
  upc   <- if (has_col("upper.CL")) "upper.CL" else if (has_col("asymp.UCL")) "asymp.UCL" else "upper.CL"
  dfcol <- if (has_col("df")) "df" else NA_character_
  
  out %>%
    mutate(
      date_chr = as.character(.data$date),
      date = suppressWarnings(as.Date(date_chr)),
      date = dplyr::coalesce(date, as.Date(date_chr, tryFormats = c("%Y-%m-%d", "%d.%m.%Y")))
    ) %>%
    transmute(
      date,
      contrast = vapply(.data$contrast, swap_contrast_label, character(1)),
      estimate = -.data$estimate,
      se = .data[[secol]],
      stat = if (!is.na(zcol)) -.data[[zcol]] else -.data$estimate / .data[[secol]],
      df = if (!is.na(dfcol)) .data[[dfcol]] else NA_real_,
      lower = -.data[[upc]],
      upper = -.data[[lowc]]
    ) %>%
    filter(is.finite(estimate), !is.na(date)) %>%
    arrange(date, contrast)
}

# Combine contrasts for the three factors that have >=2 levels in the subset
compute_combined_effects <- function(fit, df, factors = NULL, resp_var = NULL) {
  if (is.null(factors)) {
    factors <- alinv_get_treatment_factors(
      include_soil_treatment = "soiltype" %in% names(df) &&
        dplyr::n_distinct(df$soiltype, na.rm = TRUE) > 1
    )
  }
  out <- list()
  for (fac in factors) {
    if (!fac %in% names(df)) next
    if (length(levels(df[[fac]])) < 2) next
    eff <- tryCatch(extract_per_date_contrasts_any(fit, fac), error = function(e) tibble())
    if (nrow(eff)) {
      eff$effect <- fac
      out[[fac]] <- eff
    }
  }
  if (!length(out)) {
    return(tibble())
  }

  bind_rows(out) %>%
    mutate(
      response_var = resp_var %||% NA_character_,
      contrast_basis = "treatment - baseline"
    ) %>%
    alinv_apply_response_orientation(
      resp_var = resp_var,
      estimate_col = "estimate",
      lower_col = "lower",
      upper_col = "upper",
      stat_col = "stat"
    )
}

# Combined plot of effect-size time series (Robinia, Precipitation, Culture)
plot_combined_effects <- function(effects_df,
                                  title,
                                  ylab = alinv_temporal_effect_y_label(),
                                  y_limits = NULL) {
  if (!nrow(effects_df)) stop("No contrasts to plot.")
  lab_map <- c(
    robinia = "Robinia (with - without)",
    precipitation = "Precipitation (drought - control)",
    culture = "Culture (mixed - mono)",
    soiltype = "Soil (drier beech soil - wetter robinia soil)",
    extreme_event = "Extreme event (yes - no)"
  )
  
  ylowest <- effects_df$estimate - effects_df$se
  ylowest <- min(ylowest, na.rm = TRUE)
  ylowest <- ylowest - abs(ylowest * 0.1)
  if (!is.null(y_limits) && length(y_limits) == 2L) {
    ylowest <- y_limits[1]
  }

  drought_periods <- list(
    c("2025-06-20", "2025-07-02"),
    c("2025-08-12", "2025-08-20")
  )
  date_range <- range(effects_df$date, na.rm = TRUE)
  visible_periods <- purrr::keep(drought_periods, function(win) {
    win_start <- as.Date(win[[1]])
    win_end <- as.Date(win[[2]])
    win_start <= date_range[2] && win_end >= date_range[1]
  })
  
  p <- effects_df %>%
    mutate(effect_lbl = dplyr::recode(effect, !!!lab_map)) %>%
    ggplot(aes(x = date, y = estimate, ymin = lower, ymax = upper,
               group = effect_lbl, color = effect_lbl, fill = effect_lbl)) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4) +
    geom_ribbon(alpha = 0.15, color = NA) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    scale_x_date(date_breaks = "1 month", date_labels = "%b") +
    labs(x = NULL, y = ylab, title = title, color = NULL, fill = NULL) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "bottom",
      # panel.grid.major = element_line(color = "lightgrey", linewidth = 0.3),
      # panel.grid.minor = element_line(color = "lightgrey", linewidth = 0.2),
      # panel.grid.major.x = element_line(),
      # panel.grid.major.y = element_line(),
      # panel.grid.minor.x = element_line(),
      # panel.grid.minor.y = element_line(),
      # panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4)
    ) +
    guides(color = guide_legend(nrow = 2, byrow = TRUE),
           fill  = guide_legend(nrow = 2, byrow = TRUE))

  for (win in visible_periods) {
    win_start <- max(as.Date(win[[1]]), date_range[1])
    win_end <- min(as.Date(win[[2]]), date_range[2])
    if (win_start > win_end) next

    p <- p + annotate(
      "segment",
      x = win_start,
      xend = win_end,
      y = ylowest,
      yend = ylowest,
      size = 2.5,
      color = "indianred",
      lineend = "round"
    )
  }

  if (!is.null(y_limits) && length(y_limits) == 2L) {
    p <- p + coord_cartesian(ylim = y_limits)
  }

  p
}

# ---------------------- ORCHESTRATORS --------------------------

temporal_effect_cache_path <- function(type = "tree",
                                       data_name,
                                       resp_var = NULL,
                                       target_species,
                                       soil_type = "both",
                                       include_soil_treatment = NULL,
                                       add_covars = FALSE,
                                       swc_source = "measured") {
  include_soil_treatment <- alinv_resolve_include_soil_treatment(
    include_soil_treatment = include_soil_treatment,
    soil_filter = soil_type
  )

  out_dir <- alinv_data_path("model-factorial-effect", create_dir = TRUE)
  resp_tag <- if (!is.null(resp_var)) resp_var else "default"
  covar_tag <- if (isTRUE(add_covars)) "withCovars" else "noCovars"
  swc_tag <- if (swc_source == "measured") "swcMeas" else "swcImputed"
  soil_mode_tag <- alinv_soil_mode_tag(
    soil_filter = soil_type,
    include_soil_treatment = include_soil_treatment
  )

  file_name <- paste0(
    "effect-",
    type, "-",
    data_name, "-",
    resp_tag, "-",
    target_species, "-",
    soil_mode_tag, "-",
    covar_tag, "-",
    swc_tag,
    ".rds"
  )

  file.path(out_dir, file_name)
}

temporal_effect_plot_meta <- function(type = "tree",
                                      data_name,
                                      resp_var = NULL,
                                      target_species,
                                      soil_type = "both",
                                      include_soil_treatment = NULL,
                                      add_covars = FALSE,
                                      swc_source = "measured",
                                      title = NULL,
                                      y_limits = alinv_temporal_effect_y_limits()) {
  tibble::tibble(
    type = type,
    data_name = data_name,
    resp_var = resp_var %||% "default",
    species = target_species,
    soil_filter = soil_type,
    include_soil_treatment = isTRUE(include_soil_treatment),
    add_covars = isTRUE(add_covars),
    swc_source = swc_source,
    title = title %||% paste0(
      "Time-varying treatment effects on ", target_species,
      " (soil: ", soil_type, ", data: ",
      data_name, if (!is.null(resp_var)) paste0(", ", resp_var) else "", ")"
    ),
    y_axis_label = alinv_temporal_effect_y_label(),
    y_limit_lower = y_limits[[1]],
    y_limit_upper = y_limits[[2]]
  )
}

extract_temporal_glmm_performance <- function(result,
                                              data_name,
                                              resp_var,
                                              target_species,
                                              soil_type = "both",
                                              include_soil_treatment = NULL,
                                              swc_source = "measured") {
  include_soil_treatment <- alinv_resolve_include_soil_treatment(
    include_soil_treatment = include_soil_treatment,
    soil_filter = soil_type
  )

  mod <- result$model %||% NULL
  data_df <- result$data_model %||% tibble::tibble()

  safe_stat <- function(expr) {
    tryCatch(expr, error = function(e) NA_real_)
  }

  r2_vals <- alinv_lmm_r2(mod)

  tibble::tibble(
    model_scope = "Temporal GLMM",
    data_name = data_name,
    resp_var = resp_var,
    species = target_species,
    soil_filter = soil_type,
    include_soil_treatment = isTRUE(include_soil_treatment),
    swc_source = swc_source,
    n_obs = nrow(data_df),
    n_boxes = if ("boxlabel" %in% names(data_df)) dplyr::n_distinct(data_df$boxlabel) else NA_integer_,
    n_trees = if ("tree_id" %in% names(data_df)) dplyr::n_distinct(data_df$tree_id) else NA_integer_,
    n_dates = if ("date" %in% names(data_df)) dplyr::n_distinct(data_df$date) else NA_integer_,
    aic = safe_stat(stats::AIC(mod)),
    bic = safe_stat(stats::BIC(mod)),
    logLik = safe_stat(as.numeric(stats::logLik(mod))),
    sigma = safe_stat(stats::sigma(mod)),
    r2_marginal = r2_vals$r2_marginal[[1]],
    r2_conditional = r2_vals$r2_conditional[[1]]
  )
}

write_temporal_effect_csv_bundle <- function(result,
                                             export_stem,
                                             meta = NULL) {
  effects_df <- result$effects %||% tibble::tibble()
  data_model_df <- result$data_model %||% tibble::tibble()
  meta_df <- meta %||% tibble::tibble()

  readr::write_csv(effects_df, paste0(export_stem, "-effects.csv"))
  readr::write_csv(data_model_df, paste0(export_stem, "-data_model.csv"))
  readr::write_csv(meta_df, paste0(export_stem, "-meta.csv"))
}

plot_temporal_effects_from_csv <- function(effects_df,
                                           meta_df = NULL,
                                           y_limits = alinv_temporal_effect_y_limits(),
                                           title = NULL) {
  meta_value <- function(col, default = NULL) {
    if (is.null(meta_df) || !nrow(meta_df) || !col %in% names(meta_df)) {
      return(default)
    }
    meta_df[[col]][[1]]
  }

  if (is.null(effects_df) || !nrow(effects_df)) {
    return(
      ggplot() +
        theme_void() +
        ggtitle(title %||% meta_value("title", "No temporal GLMM effect data"))
    )
  }

  plot_combined_effects(
    effects_df = effects_df,
    title = title %||% meta_value("title", "Temporal GLMM effects"),
    ylab = meta_value("y_axis_label", alinv_temporal_effect_y_label()),
    y_limits = y_limits
  )
}

# Make one figure for a target species and a given data_name/resp_var
make_effect_figure_generic <- function(
    type = "tree",
    data_name = c("chlorophyll", "condition", "growth", "quantum_yield", "senescence", "phenology"),
    resp_var = NULL,                 # set if data_name has multiple options
    target_species = "fagus",
    soil_type = "both",
    include_soil_treatment = NULL,
    add_covars = FALSE,
    covars_fun = NULL,
    y_limits = NULL,
    swc_source = "measured",         # "measured" or "imputed_gam"
    force_run = FALSE                # NEW: overwrite existing results for today
) {
  data_name <- match.arg(data_name)
  include_soil_treatment <- alinv_resolve_include_soil_treatment(
    include_soil_treatment = include_soil_treatment,
    soil_filter = soil_type
  )
  
  # --------------------------------------------------
  # 0) Build cache path and filename
  # --------------------------------------------------
  cache_path <- temporal_effect_cache_path(
    type = type,
    data_name = data_name,
    resp_var = resp_var,
    target_species = target_species,
    soil_type = soil_type,
    include_soil_treatment = include_soil_treatment,
    add_covars = add_covars,
    swc_source = swc_source
  )
  export_stem <- tools::file_path_sans_ext(cache_path)
  
  # --------------------------------------------------
  # 1) Load from cache if available and not forcing rerun
  # --------------------------------------------------
  if (file.exists(cache_path) && !force_run) {
    message("Loading cached results from: ", cache_path)
    cached_result <- readRDS(cache_path)
    cache_has_empty_effects <- {
      data_model <- cached_result$data_model %||% tibble::tibble()
      effects <- cached_result$effects %||% tibble::tibble()
      nrow(data_model) > 0 && nrow(effects) == 0
    }

    if (cache_has_empty_effects) {
      message("Cached temporal results contained no extracted contrasts. Re-running: ", cache_path)
    } else {
      cached_meta <- temporal_effect_plot_meta(
        type = type,
        data_name = data_name,
        resp_var = resp_var,
        target_species = target_species,
        soil_type = soil_type,
        include_soil_treatment = include_soil_treatment,
        add_covars = add_covars,
        swc_source = swc_source,
        y_limits = alinv_temporal_effect_y_limits()
      )
      cached_result$plot <- plot_temporal_effects_from_csv(
        effects_df = cached_result$effects %||% tibble::tibble(),
        meta_df = cached_meta,
        y_limits = y_limits %||% alinv_temporal_effect_y_limits()
      )
      write_temporal_effect_csv_bundle(
        result = cached_result,
        export_stem = export_stem,
        meta = cached_meta
      )
      return(cached_result)
    }
  }
  
  # --------------------------------------------------
  # 2) Run full pipeline
  # --------------------------------------------------
  dfm <- prepare_df_generic(
    type         = type,
    data_name    = data_name,
    resp_var     = resp_var,
    species_keep = target_species,
    add_covars   = add_covars,
    covars_fun   = covars_fun,
    soil_type    = soil_type,
    include_soil_treatment = include_soil_treatment,
    swc_source   = swc_source
  )

  if (!nrow(dfm)) {
    result <- list(
      plot = alinv_empty_plot(
        title = "No temporal GLMM data",
        subtitle = paste0("No model rows available for ", target_species, " / ", data_name, if (!is.null(resp_var)) paste0(" / ", resp_var) else "")
      ),
      model = NULL,
      effects = tibble::tibble(),
      data_model = dfm
    )
    saveRDS(result, cache_path)
    write_temporal_effect_csv_bundle(
      result = result,
      export_stem = export_stem,
      meta = temporal_effect_plot_meta(
        type = type,
        data_name = data_name,
        resp_var = resp_var,
        target_species = target_species,
        soil_type = soil_type,
        include_soil_treatment = include_soil_treatment,
        add_covars = add_covars,
        swc_source = swc_source,
        y_limits = y_limits %||% alinv_temporal_effect_y_limits()
      )
    )
    return(result)
  }
  
  fit <- fit_glmm_per_date(
    dfm,
    include_covars = add_covars,
    include_soil_treatment = include_soil_treatment
  )
  eff <- compute_combined_effects(
    fit,
    dfm,
    factors = alinv_get_treatment_factors(
      include_soil_treatment = include_soil_treatment,
      soil_filter = soil_type
    ),
    resp_var = resp_var %||% if (length(.default_resp_map[[data_name]]) == 1) .default_resp_map[[data_name]] else resp_var
  )
  if (nrow(eff)) {
    eff <- eff %>%
      mutate(
        scale_mode = "z_response",
        swc_source = swc_source,
        species = target_species,
        include_interaction = FALSE,
        include_soil_treatment = include_soil_treatment
      )
  }
  
  ttl <- paste0(
    "Time-varying treatment effects on ", target_species,
    " (soil: ", soil_type, ", data: ",
    data_name, if (!is.null(resp_var)) paste0(", ", resp_var) else "", ")"
  )
  p <- if (nrow(eff)) {
    plot_combined_effects(eff, ttl, y_limits = y_limits %||% alinv_temporal_effect_y_limits())
  } else {
    alinv_empty_plot(
      title = "No temporal GLMM contrasts",
      subtitle = paste0("No treatment contrasts could be extracted for ", target_species, " / ", data_name, if (!is.null(resp_var)) paste0(" / ", resp_var) else "")
    )
  }
  
  result <- list(plot = p, model = fit, effects = eff, data_model = dfm)
  
  # --------------------------------------------------
  # 3) Save to cache
  # --------------------------------------------------
  saveRDS(result, cache_path)
  write_temporal_effect_csv_bundle(
    result = result,
    export_stem = export_stem,
    meta = temporal_effect_plot_meta(
      type = type,
      data_name = data_name,
      resp_var = resp_var,
      target_species = target_species,
      soil_type = soil_type,
      include_soil_treatment = include_soil_treatment,
      add_covars = add_covars,
      swc_source = swc_source,
      title = ttl,
      y_limits = y_limits %||% alinv_temporal_effect_y_limits()
    )
  )
  message("Saved results to: ", cache_path)
  
  result
}

# Convenience: produce a list of figures for multiple species
make_figs_for_species <- function(
    species_vec = c("fagus", "quercus", "robinia"), ...
) {
  lapply(species_vec, function(sp) {
    res <- make_effect_figure_generic(target_species = sp, ...)
    res$species <- sp
    res
  })
}

# ---------------------- EXAMPLES --------------------------------
# # Chlorophyll on Fagus
# out_chl_fagus <- make_effect_figure_generic(
#   type = "tree", data_name = "chlorophyll",
#   target_species = "fagus"
# )
# out_chl_fagus$plot
#
# Growth - HEIGHT increments on Fagus
# out_grow_h_fagus <- make_effect_figure_generic(
#   type = "tree", data_name = "growth",
#   resp_var = "diameter_inc_t0", target_species = "fagus"
# )
# out_grow_h_fagus$plot
#
# # Growth - DIAMETER increments on Quercus
# out_grow_d_quer <- make_effect_figure_generic(
#   type = "tree", data_name = "growth",
#   resp_var = "diameter", target_species = "quercus"
# )
# out_grow_d_quer$plot
#
# # Quantum yield on Fagus
# out_qy_fagus <- make_effect_figure_generic(
#   type = "tree", data_name = "quantum_yield",
#   target_species = "fagus"
# )
# out_qy_fagus$plot
#
# # Senescence (percent) on Fagus
# out_sen_pct_fagus <- make_effect_figure_generic(
#   type = "tree", data_name = "senescence",
#   resp_var = "percent_senesced", target_species = "fagus"
# )
# out_sen_pct_fagus$plot
#
# Senescence (chlavg) on Fagus
# out_sen_chl_fagus <- make_effect_figure_generic(
#   type = "tree", data_name = "senescence",
#   resp_var = "percent_senesced", target_species = "fagus"
# )
# out_sen_chl_fagus$plot
# 
# # Effects on Fagus phenology (DOY; lower DOY = earlier stage)
# out_pheno_fagus <- make_effect_figure_generic(
#   type = "tree",
#   data_name = "phenology",
#   target_species = "fagus"
# )
# out_pheno_fagus$plot
