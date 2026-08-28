# Temporal model-data preparation shared by supplementary and bootstrap code.

# ================================================================
# Generalized model-data preparation for multiple response types
# ================================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

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
  senescence    = c("remaining_green", "chlavg"),
  phenology     = "doy"                          # special handling below
)
# ---------------------- HELPERS --------------------------------

# Standardize factor baselines and drop NAs on required cols
.standardize_and_clean <- function(df, cols_needed) {
  df %>%
    dplyr::select(all_of(cols_needed)) %>%
    drop_na(date, robinia, precipitation, culture, species, tree_id, boxlabel, y) %>%
    mutate(
      date = as.Date(date),
      date = factor(date),
      robinia = factor(robinia, levels = alinv_factor_levels("robinia")),
      precipitation = factor(precipitation, levels = alinv_factor_levels("precipitation")),
      culture = factor(culture, levels = alinv_factor_levels("culture")),
      soiltype = alinv_relevel_soiltype(soiltype),
      extreme_event = factor(extreme_event, levels = alinv_factor_levels("extreme_event")),
      species = factor(species, levels = alinv_factor_levels("species"))
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
    exclude_initial_growth_baseline = FALSE
) {
  data_name <- match.arg(data_name)
  include_soil_treatment <- alinv_resolve_include_soil_treatment(
    include_soil_treatment = include_soil_treatment,
    soil_filter = soil_type
  )

  # Pull
  df_raw <- get_data(type = type, data_name = data_name)

  # Early species filter if present
  if (!is.null(species_keep) && "species" %in% names(df_raw)) {
    df_raw <- dplyr::filter(df_raw, species %in% species_keep)
  }

  df_raw <- alinv_apply_soil_context(
    df_raw,
    soil_filter = soil_type,
    include_soil_treatment = include_soil_treatment
  )

  # Increment baselines are deterministic zeros, not growth observations.
  # Opt in only for the new phase-wise temporal analysis; existing models
  # retain their established data selection.
  if (isTRUE(exclude_initial_growth_baseline)) {
    if (data_name != "growth") stop("Baseline exclusion requires growth data.")
    df_raw <- df_raw %>%
      group_by(tree_id) %>%
      filter(date > min(date, na.rm = TRUE)) %>%
      ungroup()
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
                 "precipitation", "robinia", "swc", resp_var,
                 intersect("phase", names(df_raw)))
  df <- dplyr::select(df_raw, dplyr::all_of(base_cols))

  # Response column 'y'
  df <- dplyr::rename(df, y = !!resp_var)

  # Optional covariates (e.g., SWC)
  df <- .maybe_add_covars(df, add_covars, covars_fun)

  # Standardize factors/baselines and drop NAs on required columns
  cols_needed <- c("y", "date", "date_num", "robinia", "precipitation", "culture", "soiltype", "extreme_event",
                   "species", "tree_id", "boxlabel", "swc", "phase")
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
