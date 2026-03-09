make_effect_figure_numeric <- function(
    type = "tree",
    data_name = c("chlorophyll", "condition", "growth", "quantum_yield", "senescence", "phenology"),
    resp_var = NULL,
    target_species = "fagus",
    soil_type = "both",
    add_covars = FALSE,
    covars_fun = NULL
) {
  data_name <- match.arg(data_name)
  
  dfm <- prepare_df_generic(
    type         = type,
    data_name    = data_name,
    resp_var     = resp_var,
    species_keep = target_species,
    add_covars   = add_covars,
    covars_fun   = covars_fun,
    soil_type    = soil_type
  )
  
  # global numeric-time model (all factors together)
  fit_global <- fit_glmm_time_numeric(dfm, include_covars = add_covars)
  
  # per-factor numeric-time contrasts
  eff <- compute_combined_effects_numeric(
    df  = dfm,
    include_covars = add_covars
  )
  
  ttl <- paste0(
    "Time-varying treatment effects (numeric time) on ", target_species,
    " (soil: ", soil_type, ", data: ",
    data_name, if (!is.null(resp_var)) paste0(", ", resp_var) else "", ")"
  )
  p <- if (nrow(eff)) plot_combined_effects(eff, ttl) else NULL
  
  list(
    plot       = p,
    model      = fit_global,  # global model
    effects    = eff,
    data_model = dfm
  )
}

fit_glmm_time_numeric <- function(df, include_covars = FALSE) {
  # base time + interactions
  rhs_terms <- c(
    "date_num",
    "date_num:robinia",
    "date_num:precipitation",
    "date_num:culture",
    "date_num:soiltype",
    "extreme_event"
  )
  
  # mapping from interaction terms to underlying factors
  term_factor_map <- c(
    "date_num:robinia"       = "robinia",
    "date_num:precipitation" = "precipitation",
    "date_num:culture"       = "culture",
    "date_num:soiltype"      = "soiltype",
    "extreme_event"          = "extreme_event"
  )
  
  # drop interaction/main terms whose factor has < 2 levels or is missing
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
  
  rhs <- paste(rhs_terms, collapse = " + ")
  
  if (isTRUE(include_covars) && "swc_sc" %in% names(df)) {
    rhs <- paste(rhs, "+ swc_sc")
  }
  
  fml <- as.formula(paste("y ~", rhs, "+ (1|boxlabel) + (1|tree_id)"))
  lme4::lmer(fml, data = df, REML = TRUE)
}

fit_glmm_time_numeric_factor <- function(df, factor_name, include_covars = FALSE) {
  if (!factor_name %in% names(df)) {
    stop("Factor '", factor_name, "' not found in data.")
  }
  if (!is.factor(df[[factor_name]])) {
    stop("'", factor_name, "' must be a factor.")
  }
  if (length(levels(df[[factor_name]])) < 2) {
    stop("Factor '", factor_name, "' has fewer than 2 levels.")
  }
  
  rhs <- paste0("date_num * ", factor_name)
  if (isTRUE(include_covars) && "swc_sc" %in% names(df)) {
    rhs <- paste(rhs, "+ swc_sc")
  }
  
  fml <- as.formula(paste("y ~", rhs, "+ (1|boxlabel) + (1|tree_id)"))
  lme4::lmer(fml, data = df, REML = TRUE)
}

extract_per_date_contrasts_numeric_factor <- function(fit, factor_name, df) {
  # grid of numeric times where we want contrasts
  time_vals <- sort(unique(df$date_num))
  
  emm <- emmeans::emmeans(
    fit,
    specs = as.formula(paste("~", factor_name, "| date_num")),
    at    = list(date_num = time_vals),
    type  = "link"
  )
  
  ct <- emmeans::contrast(
    emm,
    method = "revpairwise",   # baseline − treatment (because 'rev')
    by     = "date_num",
    adjust = "none"
  )
  
  sm  <- suppressMessages(suppressWarnings(summary(ct, infer = TRUE)))
  out <- tibble::as_tibble(sm)
  
  has_col <- function(nm) nm %in% names(out)
  zcol  <- if (has_col("z.ratio")) "z.ratio" else if (has_col("t.ratio")) "t.ratio" else NA_character_
  secol <- "SE"
  lowc  <- if (has_col("lower.CL")) "lower.CL" else if (has_col("asymp.LCL")) "asymp.LCL" else "lower.CL"
  upc   <- if (has_col("upper.CL")) "upper.CL" else if (has_col("asymp.UCL")) "asymp.UCL" else "upper.CL"
  dfcol <- if (has_col("df")) "df" else NA_character_
  
  out %>%
    dplyr::mutate(
      date = as.Date(.data$date_num, origin = "1970-01-01")
    ) %>%
    dplyr::transmute(
      date,
      contrast = .data$contrast,
      estimate = .data$estimate,
      se       = .data[[secol]],
      stat     = if (!is.na(zcol)) .data[[zcol]] else NA_real_,
      df       = if (!is.na(dfcol)) .data[[dfcol]] else NA_real_,
      lower    = .data[[lowc]],
      upper    = .data[[upc]]
    ) %>%
    dplyr::filter(is.finite(estimate), !is.na(date)) %>%
    dplyr::arrange(date, contrast)
}

compute_combined_effects_numeric <- function(
    df,
    include_covars = FALSE,
    factors = c("robinia", "precipitation", "culture", "soiltype", "extreme_event")
) {
  out <- list()
  
  for (fac in factors) {
    if (!fac %in% names(df)) next
    if (!is.factor(df[[fac]]) || length(levels(df[[fac]])) < 2) next
    
    eff <- tryCatch({
      fit_fac <- fit_glmm_time_numeric_factor(df, fac, include_covars = include_covars)
      extract_per_date_contrasts_numeric_factor(fit_fac, fac, df)
    }, error = function(e) {
      message("Numeric-time contrasts failed for factor '", fac, "': ", e$message)
      tibble::tibble()
    })
    
    if (nrow(eff)) {
      eff$effect <- fac
      out[[fac]] <- eff
    }
  }
  
  if (length(out)) dplyr::bind_rows(out) else tibble::tibble()
}
