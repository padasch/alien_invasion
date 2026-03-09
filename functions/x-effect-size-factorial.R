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
  library(lmerTest)
  library(emmeans)
  library(lubridate)
})
emmeans::emm_options(lmer.df = "asymptotic")  # quiet, fast CIs

# ---------------------- CONFIG / MAPPING ------------------------

# Map data_name -> default response columns
.default_resp_map <- list(
  chlorophyll   = "chl",
  condition     = "condition",
  growth        = c("height", "diameter"),       # requires increment!
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
      soiltype = factor(soiltype, levels = c("inoc-beech", "inoc-robinia")),
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
    add_covars = FALSE,
    covars_fun = NULL,                # function returning covariates (boxlabel+date)
    soil_type = "both",
    swc_source = "measured"          # "measured" or "imputed_gam"
) {
  data_name <- match.arg(data_name)
  
  # Pull 
  df_raw <- get_data(type = type, data_name = data_name, swc_source = swc_source)
  
  # Early species filter if present
  if (!is.null(species_keep) && "species" %in% names(df_raw)) {
    df_raw <- dplyr::filter(df_raw, species %in% species_keep)
  }
  
  # Early soil type filter
  if (soil_type != "both") {
    df_raw <- dplyr::filter(df_raw, soiltype == soil_type)
  }
  
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

  # Keep a standardized SWC covariate available for temporal GLMMs.
  # This allows add_covars=TRUE to include SWC regardless of an external covariate join.
  if ("swc" %in% names(df)) {
    df <- df %>% mutate(swc_sc = as.numeric(scale(swc)))
  }
  
  df
}

# ---------------------- MODEL (UNCHANGED) ----------------------

# Per-date GLMM with date main effect and interactions (same as before)
fit_glmm_per_date <- function(df, include_covars = FALSE) {
  # start with all candidate terms
  # rhs_terms <- c("date", "date:robinia", "date:precipitation", "date:culture", "date:soiltype", "extreme_event")
  rhs_terms <- c("date", "date:robinia", "date:precipitation", "date:culture", "date:soiltype")
  
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

# Robust per-date contrasts for ANY categorical factor; use revpairwise so
# effect = baseline − treatment (matches your interpretation)
extract_per_date_contrasts_any <- function(fit, factor_name) {
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
      contrast = .data$contrast,
      estimate = .data$estimate,
      se = .data[[secol]],
      stat = if (!is.na(zcol)) .data[[zcol]] else NA_real_,
      df = if (!is.na(dfcol)) .data[[dfcol]] else NA_real_,
      lower = .data[[lowc]],
      upper = .data[[upc]]
    ) %>%
    filter(is.finite(estimate), !is.na(date)) %>%
    arrange(date, contrast)
}

# Combine contrasts for the three factors that have >=2 levels in the subset
compute_combined_effects <- function(fit, df, factors = c("robinia", "precipitation", "culture", "soiltype", "extreme_event")) {
  out <- list()
  for (fac in factors) {
    if (length(levels(df[[fac]])) < 2) next
    eff <- tryCatch(extract_per_date_contrasts_any(fit, fac), error = function(e) tibble())
    if (nrow(eff)) {
      eff$effect <- fac
      out[[fac]] <- eff
    }
  }
  if (length(out)) bind_rows(out) else tibble()
}

# Combined plot of effect-size time series (Robinia, Precipitation, Culture)
plot_combined_effects <- function(effects_df, title, ylab = "Effect size (baseline − treatment)") {
  if (!nrow(effects_df)) stop("No contrasts to plot.")
  lab_map <- c(
    robinia = "Robinia (without − with)",
    precipitation = "Drought (control − drought)",
    culture = "Culture (mono − mixed)",
    soiltype = "Soil Type (normal - robinia inoc)",
    extreme_event = "Extreme Event (no - yes)"
  )
  
  ylowest <- effects_df$estimate - effects_df$se
  ylowest <- min(ylowest)
  ylowest <- ylowest - abs(ylowest * 0.1)
  
  effects_df %>%
    mutate(effect_lbl = dplyr::recode(effect, !!!lab_map)) %>%
    ggplot(aes(x = date, y = estimate, ymin = lower, ymax = upper,
               group = effect_lbl, color = effect_lbl, fill = effect_lbl)) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4) +
    geom_ribbon(alpha = 0.15, color = NA) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    annotate("segment",
             x = as.Date("2025-06-20"), xend = as.Date("2025-07-02"),
             y = ylowest, yend = ylowest,
             size = 2.5, color = "indianred", lineend = "round"
    ) +
    annotate("segment",
             x = as.Date("2025-08-12"), xend = as.Date("2025-08-20"),
             y = ylowest, yend = ylowest,
             size = 2.5, color = "indianred", lineend = "round"
    ) +
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
}

# ---------------------- ORCHESTRATORS --------------------------

# Make one figure for a target species and a given data_name/resp_var
make_effect_figure_generic <- function(
    type = "tree",
    data_name = c("chlorophyll", "condition", "growth", "quantum_yield", "senescence", "phenology"),
    resp_var = NULL,                 # set if data_name has multiple options
    target_species = "fagus",
    soil_type = "both",
    add_covars = FALSE,
    covars_fun = NULL,
    swc_source = "measured",         # "measured" or "imputed_gam"
    force_run = FALSE                # NEW: overwrite existing results for today
) {
  data_name <- match.arg(data_name)
  
  # --------------------------------------------------
  # 0) Build cache path and filename
  # --------------------------------------------------
  today_str <- format(Sys.Date(), "%Y-%m-%d")
  
  out_dir <- file.path("output", today_str, "model-factorial-effect")
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  resp_tag <- if (!is.null(resp_var)) resp_var else "default"
  covar_tag <- if (isTRUE(add_covars)) "withCovars" else "noCovars"
  swc_tag <- if (swc_source == "measured") "swcMeas" else "swcImputed"
  
  file_name <- paste0(
    "effect-",
    type, "-",
    data_name, "-",
    resp_tag, "-",
    target_species, "-",
    soil_type, "-",
    covar_tag, "-",
    swc_tag,
    ".rds"
  )
  
  cache_path <- file.path(out_dir, file_name)
  
  # --------------------------------------------------
  # 1) Load from cache if available and not forcing rerun
  # --------------------------------------------------
  if (file.exists(cache_path) && !force_run) {
    message("Loading cached results from: ", cache_path)
    return(readRDS(cache_path))
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
    swc_source   = swc_source
  )
  
  fit <- fit_glmm_per_date(dfm, include_covars = add_covars)
  eff <- compute_combined_effects(fit, dfm)
  
  ttl <- paste0(
    "Time-varying treatment effects on ", target_species,
    " (soil: ", soil_type, ", data: ",
    data_name, if (!is.null(resp_var)) paste0(", ", resp_var) else "", ")"
  )
  p <- if (nrow(eff)) plot_combined_effects(eff, ttl) else NULL
  
  result <- list(plot = p, model = fit, effects = eff, data_model = dfm)
  
  # --------------------------------------------------
  # 3) Save to cache
  # --------------------------------------------------
  saveRDS(result, cache_path)
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
