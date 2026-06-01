
# Packages ------------------------------------------------------------------------------------

library(tidyverse)
library(ggplot2)
if (requireNamespace("ggalluvial", quietly = TRUE)) {
  library(ggalluvial)
}
library(readxl)
library(lubridate)
library(rlang)
library(dplyr)
library(conflicted)
if (requireNamespace("ggpubr", quietly = TRUE)) {
  library(ggpubr)
}
if (requireNamespace("ggpattern", quietly = TRUE)) {
  library(ggpattern)
}



# Conflicts -----------------------------------------------------------------------------------

conflicts_prefer(
  dplyr::select,
  dplyr::filter,
  dplyr::mutate
)

if (!exists("clean_names", mode = "function")) {
  clean_names <- function(data) {
    out <- data
    nm <- names(out)
    nm <- gsub("([a-z0-9])([A-Z])", "\\1_\\2", nm)
    nm <- tolower(nm)
    nm <- gsub("[^a-z0-9]+", "_", nm)
    nm <- gsub("^_+|_+$", "", nm)
    nm <- make.unique(nm, sep = "_")
    names(out) <- nm
    out
  }
}

if (!exists("remove_empty", mode = "function")) {
  remove_empty <- function(data, which = c("rows", "cols")) {
    which <- match.arg(which)
    is_empty <- function(x) {
      if (is.list(x)) {
        lengths(x) == 0
      } else if (is.character(x)) {
        is.na(x) | trimws(x) == ""
      } else {
        is.na(x)
      }
    }

    if (identical(which, "cols")) {
      keep <- vapply(data, function(col) !all(is_empty(col)), logical(1))
      data[, keep, drop = FALSE]
    } else {
      keep <- apply(data, 1, function(row) !all(is_empty(row)))
      data[keep, , drop = FALSE]
    }
  }
}

# Path helpers -------------------------------------------------------------------------------

# Resolve project root by finding a directory that contains data/interim.
.alinv_project_root <- function() {
  wd <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  candidates <- c(wd)
  for (i in 1:6) {
    candidates <- c(candidates, dirname(candidates[length(candidates)]))
  }
  for (root in unique(candidates)) {
    if (dir.exists(file.path(root, "data", "interim"))) return(root)
  }
  wd
}

# Escape regex metacharacters in literal paths.
.escape_regex <- function(x) {
  gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x)
}

# Resolve possibly relative paths against project root when needed.
.resolve_path <- function(path) {
  if (is.null(path) || !nzchar(path)) return(path)
  if (grepl("^(/|[A-Za-z]:[\\\\/])", path)) {
    return(normalizePath(path, winslash = "/", mustWork = FALSE))
  }

  project_root <- .alinv_project_root()
  project_path <- file.path(project_root, path)

  if (dir.exists(project_path) || file.exists(project_path)) {
    return(normalizePath(project_path, winslash = "/", mustWork = TRUE))
  }

  if (dir.exists(path) || file.exists(path)) {
    return(normalizePath(path, winslash = "/", mustWork = TRUE))
  }

  normalizePath(project_path, winslash = "/", mustWork = FALSE)
}

alinv_project_relative_path <- function(path) {
  if (is.null(path) || !nzchar(path)) return(path)

  project_root <- normalizePath(.alinv_project_root(), winslash = "/", mustWork = TRUE)
  path_abs <- normalizePath(.resolve_path(path), winslash = "/", mustWork = FALSE)

  if (identical(path_abs, project_root)) {
    return(".")
  }

  prefix <- paste0(project_root, "/")
  if (!startsWith(path_abs, prefix)) {
    return(path_abs)
  }

  sub(paste0("^", .escape_regex(prefix)), "", path_abs)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

alinv_config_path <- file.path(.alinv_project_root(), "config", "analysis-config.R")
if (file.exists(alinv_config_path)) {
  source(alinv_config_path, local = FALSE)
} else {
  stop("Missing analysis config: ", alinv_config_path, call. = FALSE)
}

alinv_clean_names <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("([a-z0-9])([A-Z])", "\\1_\\2", x)
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- tolower(x)
  x <- gsub("^_+|_+$", "", x)
  x[x == ""] <- "x"
  make.unique(x, sep = "_")
}

alinv_clean_names_df <- function(df) {
  names(df) <- alinv_clean_names(names(df))
  df
}

alinv_drop_empty_cols <- function(df) {
  keep <- vapply(
    df,
    function(col) {
      if (is.character(col)) {
        return(any(!is.na(col) & nzchar(trimws(col))))
      }
      any(!is.na(col))
    },
    logical(1)
  )
  df[, keep, drop = FALSE]
}

alinv_temporal_effect_y_limits <- function() {
  c(-2, 2)
}

alinv_temporal_effect_y_label <- function() {
  "Effect size (positive = beneficial, negative = harmful)"
}

alinv_empty_plot <- function(title, subtitle = NULL) {
  ggplot2::ggplot() +
    ggplot2::theme_void() +
    ggplot2::labs(title = title, subtitle = subtitle) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 10)
    )
}

alinv_lmm_r2 <- function(mod) {
  if (is.null(mod)) {
    return(tibble::tibble(r2_marginal = NA_real_, r2_conditional = NA_real_))
  }

  if (requireNamespace("performance", quietly = TRUE)) {
    r2_obj <- tryCatch(
      performance::r2_nakagawa(mod),
      error = function(e) NULL
    )
    if (!is.null(r2_obj)) {
      return(tibble::tibble(
        r2_marginal = unname(r2_obj$R2_marginal %||% NA_real_),
        r2_conditional = unname(r2_obj$R2_conditional %||% NA_real_)
      ))
    }
  }

  if (requireNamespace("MuMIn", quietly = TRUE)) {
    r2_obj <- tryCatch(
      MuMIn::r.squaredGLMM(mod),
      error = function(e) NULL
    )
    if (!is.null(r2_obj)) {
      return(tibble::tibble(
        r2_marginal = unname(r2_obj[1, "R2m"]),
        r2_conditional = unname(r2_obj[1, "R2c"])
      ))
    }
  }

  pred_fixed <- tryCatch(
    stats::predict(mod, re.form = NA),
    error = function(e) NULL
  )

  var_fixed <- if (is.null(pred_fixed)) {
    NA_real_
  } else {
    stats::var(as.numeric(pred_fixed), na.rm = TRUE)
  }

  vc_df <- tryCatch(
    as.data.frame(lme4::VarCorr(mod)),
    error = function(e) NULL
  )

  if (is.null(vc_df) || !nrow(vc_df) || !is.finite(var_fixed)) {
    return(tibble::tibble(r2_marginal = NA_real_, r2_conditional = NA_real_))
  }

  var_random <- sum(vc_df$vcov[vc_df$grp != "Residual"], na.rm = TRUE)
  var_residual <- vc_df$vcov[vc_df$grp == "Residual"][1]
  total_var <- var_fixed + var_random + var_residual

  if (!is.finite(total_var) || total_var <= 0) {
    return(tibble::tibble(r2_marginal = NA_real_, r2_conditional = NA_real_))
  }

  tibble::tibble(
    r2_marginal = var_fixed / total_var,
    r2_conditional = (var_fixed + var_random) / total_var
  )
}

ALINV_RESPONSE_DISPLAY_SIGN <- stats::setNames(numeric(), character())

alinv_response_display_sign <- function(resp_var) {
  signs <- ALINV_RESPONSE_DISPLAY_SIGN[as.character(resp_var)]
  signs[is.na(signs)] <- 1
  as.numeric(signs)
}

alinv_apply_response_orientation <- function(df,
                                             resp_var = NULL,
                                             resp_col = NULL,
                                             estimate_col = "estimate",
                                             lower_col = NULL,
                                             upper_col = NULL,
                                             stat_col = NULL,
                                             estimate_sig_col = NULL) {
  if (is.null(df) || !nrow(df)) {
    return(df)
  }

  if (!is.null(resp_col)) {
    if (!resp_col %in% names(df)) {
      stop("resp_col not found in data frame: ", resp_col, call. = FALSE)
    }
    resp_vals <- df[[resp_col]]
  } else if (!is.null(resp_var)) {
    resp_vals <- rep(resp_var, nrow(df))
  } else {
    stop("Provide either resp_var or resp_col.", call. = FALSE)
  }

  sign_vec <- alinv_response_display_sign(resp_vals)
  df$display_sign <- sign_vec

  orient_bounds <- function(lower_vals, upper_vals, signs) {
    tibble::tibble(
      lower = ifelse(signs < 0, upper_vals * signs, lower_vals * signs),
      upper = ifelse(signs < 0, lower_vals * signs, upper_vals * signs)
    )
  }

  if (!is.null(estimate_col) && estimate_col %in% names(df)) {
    raw_col <- paste0(estimate_col, "_raw")
    df[[raw_col]] <- df[[estimate_col]]
    df[[estimate_col]] <- df[[estimate_col]] * sign_vec
  }

  if (!is.null(lower_col) && !is.null(upper_col) &&
      lower_col %in% names(df) && upper_col %in% names(df)) {
    lower_raw_col <- paste0(lower_col, "_raw")
    upper_raw_col <- paste0(upper_col, "_raw")
    df[[lower_raw_col]] <- df[[lower_col]]
    df[[upper_raw_col]] <- df[[upper_col]]

    bounds <- orient_bounds(df[[lower_col]], df[[upper_col]], sign_vec)
    df[[lower_col]] <- bounds$lower
    df[[upper_col]] <- bounds$upper
  }

  if (!is.null(stat_col) && stat_col %in% names(df)) {
    raw_col <- paste0(stat_col, "_raw")
    df[[raw_col]] <- df[[stat_col]]
    df[[stat_col]] <- df[[stat_col]] * sign_vec
  }

  if (!is.null(estimate_sig_col) && estimate_sig_col %in% names(df)) {
    raw_col <- paste0(estimate_sig_col, "_raw")
    df[[raw_col]] <- df[[estimate_sig_col]]
    df[[estimate_sig_col]] <- df[[estimate_sig_col]] * sign_vec
  }

  df
}

ALINV_SOIL_LABELS <- alinv_level_labels("soiltype")

alinv_scenario_grid <- function() {
  ALINV_SCENARIO_GRID
}

alinv_get_analysis_context <- function() {
  getOption("alinv.analysis_context", NULL)
}

alinv_resolve_soil_filter <- function(soil_filter = NULL) {
  ctx <- alinv_get_analysis_context()
  soil_filter %||% ctx$soil_filter %||% "both"
}

alinv_resolve_include_soil_treatment <- function(include_soil_treatment = NULL,
                                                 soil_filter = NULL) {
  soil_filter <- alinv_resolve_soil_filter(soil_filter)
  if (!identical(soil_filter, "both")) {
    return(FALSE)
  }

  ctx <- alinv_get_analysis_context()
  include_soil_treatment %||% ctx$include_soil_treatment %||% TRUE
}

alinv_should_show_soil_panels <- function(soil_filter = NULL,
                                          include_soil_treatment = NULL) {
  soil_filter <- alinv_resolve_soil_filter(soil_filter)
  include_soil_treatment <- alinv_resolve_include_soil_treatment(
    include_soil_treatment = include_soil_treatment,
    soil_filter = soil_filter
  )

  identical(soil_filter, "both") && isTRUE(include_soil_treatment)
}

alinv_soil_mode_tag <- function(soil_filter = NULL,
                                include_soil_treatment = NULL) {
  soil_filter <- alinv_resolve_soil_filter(soil_filter)
  include_soil_treatment <- alinv_resolve_include_soil_treatment(
    include_soil_treatment = include_soil_treatment,
    soil_filter = soil_filter
  )

  if (!identical(soil_filter, "both")) {
    return(paste0("soil-", gsub("[^a-zA-Z0-9]+", "_", soil_filter)))
  }

  if (isTRUE(include_soil_treatment)) {
    "soil-both_with_soil_treatment"
  } else {
    "soil-both_without_soil_treatment"
  }
}

alinv_get_treatment_factors <- function(include_soil_treatment = NULL,
                                        soil_filter = NULL) {
  include_soil_treatment <- alinv_resolve_include_soil_treatment(
    include_soil_treatment = include_soil_treatment,
    soil_filter = soil_filter
  )

  cfg <- alinv_treatment_config(
    include_soil_treatment = include_soil_treatment,
    soil_filter = soil_filter
  )
  cfg$effect
}

alinv_treatment_config <- function(include_soil_treatment = NULL,
                                   soil_filter = NULL) {
  include_soil_treatment <- alinv_resolve_include_soil_treatment(
    include_soil_treatment = include_soil_treatment,
    soil_filter = soil_filter
  )

  cfg <- ALINV_TREATMENT_CONFIG
  if (!isTRUE(include_soil_treatment)) {
    cfg <- dplyr::filter(cfg, .data$effect != "soiltype")
  }
  cfg
}

alinv_set_analysis_context <- function(
    scenario_label = NULL,
    soil_filter = "both",
    include_soil_treatment = NULL,
    analysis_date = Sys.Date(),
    output_root = "output",
    create_dirs = TRUE
) {
  soil_filter <- alinv_resolve_soil_filter(soil_filter)
  include_soil_treatment <- alinv_resolve_include_soil_treatment(
    include_soil_treatment = include_soil_treatment,
    soil_filter = soil_filter
  )

  scenario_label <- scenario_label %||% if (!identical(soil_filter, "both")) {
    ALINV_SOIL_LABELS[[soil_filter]]
  } else if (isTRUE(include_soil_treatment)) {
    "both soils (with soil as treatment)"
  } else {
    "both soils (without soil as treatment)"
  }

  analysis_date <- as.character(as.Date(analysis_date))
  output_root <- .resolve_path(output_root)
  analysis_root <- file.path(output_root, analysis_date, scenario_label)
  analysis_data_root <- file.path(analysis_root, "data")
  notebooks_root <- file.path(analysis_root, "notebooks")

  ctx <- list(
    scenario_label = scenario_label,
    soil_filter = soil_filter,
    include_soil_treatment = include_soil_treatment,
    analysis_date = analysis_date,
    output_root = output_root,
    analysis_root = analysis_root,
    analysis_data_root = analysis_data_root,
    notebooks_root = notebooks_root
  )

  options(alinv.analysis_context = ctx)

  if (isTRUE(create_dirs)) {
    dir.create(analysis_root, recursive = TRUE, showWarnings = FALSE)
    dir.create(analysis_data_root, recursive = TRUE, showWarnings = FALSE)
    dir.create(notebooks_root, recursive = TRUE, showWarnings = FALSE)
  }

  ctx
}

alinv_init_notebook_context <- function(params = NULL,
                                        output_root = "output") {
  params <- params %||% list()

  alinv_set_analysis_context(
    scenario_label = params$scenario_label %||% NULL,
    soil_filter = params$soil_filter %||% "both",
    include_soil_treatment = params$include_soil_treatment %||% NULL,
    analysis_date = params$analysis_date %||% Sys.Date(),
    output_root = params$output_root %||% output_root
  )
}

alinv_analysis_path <- function(..., create_dir = FALSE) {
  ctx <- alinv_get_analysis_context()
  if (is.null(ctx)) {
    ctx <- alinv_set_analysis_context(create_dirs = FALSE)
  }

  path <- file.path(ctx$analysis_root, ...)
  if (isTRUE(create_dir)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  path
}

alinv_data_path <- function(..., create_dir = FALSE) {
  ctx <- alinv_get_analysis_context()
  if (is.null(ctx)) {
    ctx <- alinv_set_analysis_context(create_dirs = FALSE)
  }

  path <- file.path(ctx$analysis_data_root, ...)
  if (isTRUE(create_dir)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  path
}

alinv_resolve_output_date_dir <- function(output_root = "output",
                                          requested_date = NULL) {
  output_root <- .resolve_path(output_root)

  if (!is.null(requested_date) && nzchar(requested_date)) {
    date_dir <- file.path(output_root, requested_date)
    if (!dir.exists(date_dir)) {
      stop("Requested output date folder does not exist: ", date_dir, call. = FALSE)
    }
    return(date_dir)
  }

  date_dirs <- list.dirs(output_root, full.names = TRUE, recursive = FALSE)
  date_tbl <- tibble(
    path = date_dirs,
    label = basename(date_dirs),
    date = suppressWarnings(as.Date(basename(date_dirs)))
  ) %>%
    filter(!is.na(.data$date)) %>%
    arrange(desc(.data$date))

  if (!nrow(date_tbl)) {
    stop("No dated output folders found under ", output_root, call. = FALSE)
  }

  date_tbl$path[[1]]
}

alinv_cache_path <- function(subdir, ..., create_dir = TRUE) {
  dir_path <- alinv_data_path(subdir, create_dir = create_dir)
  file.path(dir_path, ...)
}

alinv_notebook_path <- function(..., create_dir = FALSE) {
  ctx <- alinv_get_analysis_context()
  if (is.null(ctx)) {
    ctx <- alinv_set_analysis_context(create_dirs = FALSE)
  }

  path <- file.path(ctx$notebooks_root, ...)
  if (isTRUE(create_dir)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  path
}

alinv_filter_by_soil <- function(df,
                                 soil_filter = NULL,
                                 soil_col = "soiltype") {
  soil_filter <- alinv_resolve_soil_filter(soil_filter)
  if (!soil_col %in% names(df)) {
    return(df)
  }

  if (!identical(soil_filter, "both")) {
    df <- df %>% dplyr::filter(.data[[soil_col]] == soil_filter)
  }

  df
}

alinv_relevel_soiltype <- function(x) {
  factor(x, levels = alinv_factor_levels("soiltype"))
}

alinv_apply_soil_context <- function(df,
                                     soil_filter = NULL,
                                     include_soil_treatment = NULL,
                                     soil_col = "soiltype") {
  soil_filter <- alinv_resolve_soil_filter(soil_filter)
  include_soil_treatment <- alinv_resolve_include_soil_treatment(
    include_soil_treatment = include_soil_treatment,
    soil_filter = soil_filter
  )

  df <- alinv_filter_by_soil(df, soil_filter = soil_filter, soil_col = soil_col)

  if (soil_col %in% names(df)) {
    df[[soil_col]] <- alinv_relevel_soiltype(df[[soil_col]])
    if (!isTRUE(include_soil_treatment) && identical(soil_filter, "both")) {
      df[[soil_col]] <- droplevels(df[[soil_col]])
    }
  }

  df
}

# Build a corrected site-level daily precipitation series from MeteoSwiss 10-min
# data, replacing the known corrupted 2025-07-13 to 2025-07-15 window with LWF.
get_site_precipitation_daily <- function(
  ms_file = "data/raw/sensor_data/precipitation_10m_MS.dat",
  lwf_file = "data/raw/sensor_data/precipitation_10m_LWF.csv",
  replacement_dates = seq(as.Date("2025-07-13"), as.Date("2025-07-15"), by = "day"),
  export_path = NULL
) {
  ms_file <- .resolve_path(ms_file)
  lwf_file <- .resolve_path(lwf_file)

  ms_daily <- read.table(
    ms_file,
    skip = 8,
    header = FALSE,
    fill = TRUE,
    stringsAsFactors = FALSE,
    na.strings = c("32767", "32767.0"),
    col.names = c("sta", "year", "month", "day", "hour", "minute", "precip_mm")
  ) %>%
    mutate(
      year = suppressWarnings(as.integer(year)),
      month = suppressWarnings(as.integer(month)),
      day = suppressWarnings(as.integer(day)),
      hour = suppressWarnings(as.integer(hour)),
      minute = suppressWarnings(as.integer(minute)),
      precip_mm = suppressWarnings(as.numeric(precip_mm))
    ) %>%
    filter(!is.na(year), !is.na(month), !is.na(day), !is.na(hour), !is.na(minute)) %>%
    transmute(
      date = as.Date(sprintf("%04d-%02d-%02d", year, month, day)),
      precip_mm
    ) %>%
    group_by(date) %>%
    summarise(
      precip_mm_ms = if (all(is.na(precip_mm))) NA_real_ else sum(precip_mm, na.rm = TRUE),
      ms_records = n(),
      ms_missing = sum(is.na(precip_mm)),
      .groups = "drop"
    )

  if (file.exists(lwf_file)) {
    lwf_raw <- readr::read_csv(lwf_file, show_col_types = FALSE)
    if (ncol(lwf_raw) < 2) {
      stop("LWF precipitation file must contain at least two columns.", call. = FALSE)
    }

    names(lwf_raw)[1:2] <- c("timestamp", "precip_mm")

    lwf_daily <- lwf_raw %>%
      transmute(
        datetime = lubridate::ymd_hms(timestamp, quiet = TRUE),
        date = as.Date(datetime),
        precip_mm = suppressWarnings(readr::parse_number(as.character(precip_mm)))
      ) %>%
      filter(!is.na(datetime)) %>%
      group_by(date) %>%
      summarise(
        precip_mm_lwf = if (all(is.na(precip_mm))) NA_real_ else sum(precip_mm, na.rm = TRUE),
        lwf_records = n(),
        lwf_missing = sum(is.na(precip_mm)),
        .groups = "drop"
      )
  } else {
    lwf_daily <- tibble(
      date = as.Date(character()),
      precip_mm_lwf = numeric(),
      lwf_records = integer(),
      lwf_missing = integer()
    )
  }

  precip_daily <- ms_daily %>%
    full_join(lwf_daily, by = "date") %>%
    arrange(date) %>%
    mutate(
      use_lwf = date %in% replacement_dates & !is.na(precip_mm_lwf),
      precip_mm = case_when(
        use_lwf ~ precip_mm_lwf,
        TRUE ~ precip_mm_ms
      ),
      precip_source = case_when(
        use_lwf ~ "lwf_replacement",
        !is.na(precip_mm_ms) ~ "ms",
        !is.na(precip_mm_lwf) ~ "lwf_only",
        TRUE ~ NA_character_
      ),
      # Keep a short cumulative variant for future lag sensitivity checks.
      precip_mm_2d = if_else(
        is.na(precip_mm),
        NA_real_,
        precip_mm + coalesce(dplyr::lag(precip_mm, 1), 0)
      )
    ) %>%
    select(
      date,
      precip_mm,
      precip_mm_2d,
      precip_source,
      precip_mm_ms,
      precip_mm_lwf,
      ms_records,
      ms_missing,
      lwf_records,
      lwf_missing
    )

  replacement_rows <- precip_daily %>%
    filter(date %in% replacement_dates) %>%
    select(date, precip_source)

  wrong_source_idx <- is.na(replacement_rows$precip_source) |
    replacement_rows$precip_source != "lwf_replacement"

  if (nrow(replacement_rows) != length(replacement_dates) ||
      any(wrong_source_idx)) {
    missing_dates <- setdiff(replacement_dates, replacement_rows$date)
    wrong_source_dates <- replacement_rows$date[wrong_source_idx]
    stop(
      paste(
        "Could not apply the required LWF precipitation replacements.",
        if (length(missing_dates) > 0) {
          paste("Missing date(s):", paste(as.character(missing_dates), collapse = ", "))
        } else {
          NULL
        },
        if (length(wrong_source_dates) > 0) {
          paste("Non-LWF source on:", paste(as.character(wrong_source_dates), collapse = ", "))
        } else {
          NULL
        }
      ),
      call. = FALSE
    )
  }

  if (!is.null(export_path) && nzchar(export_path)) {
    readr::write_csv(precip_daily, .resolve_path(export_path))
  }

  precip_daily
}

# Manual Functions ----------------------------------------------------------------------------

# Main Functions
get_meta <- function(which = c("tree", "box")) {
  which <- match.arg(which)
  interim_dir <- .resolve_path("data/interim")

  if (which == "tree") {
    df <- left_join(
        read_csv(file.path(interim_dir, "meta_tree.csv"), show_col_types = FALSE),
        read_csv(file.path(interim_dir, "meta_box.csv"), show_col_types = FALSE),
        by="boxlabel"
      ) |> 
      mutate(tree_id = as.character(tree_id))
  } else if (which == "box") {
    df <- read_csv(file.path(interim_dir, "meta_box.csv"), show_col_types = FALSE)
  } else {
    stop("Invalid 'which' argument. Use 'tree' or 'box'.")
  }
  return(df)
}

get_data <- function(type = c("tree", "box"), data_name, with_meta = TRUE, path = "./data/interim", swc_source = "measured") {
  type <- match.arg(type)
  swc_source <- match.arg(swc_source, choices = c("measured", "imputed_gam"))
  path <- .resolve_path(path)
  if (missing(data_name) || !nzchar(data_name)) {
    stop("Provide data_name (e.g., 'growth', 'respiration', 'chlorophyll').", call. = FALSE)
  }
  
  # Files in the target folder
  files <- list.files(path, pattern = "\\.csv$", full.names = FALSE)
  
  # Build expected names
  data_file <- sprintf("%s_%s.csv", type, data_name)   # e.g. tree_growth.csv
  meta_file <- sprintf("meta_%s.csv", type)            # e.g. meta_tree.csv
  
  # Helper: build a suggestion list from available files
  suggestion_text <- {
    m <- str_match(files, "^(tree|box)_(.+)\\.csv$")
    if (nrow(m) == 0) {
      paste("No tree_*.csv or box_*.csv found. Available CSVs:\n",
            paste(paste0("• ", files), collapse = "\n"))
    } else {
      df <- tibble(type = m[,2], data_name = m[,3]) %>% filter(!is.na(type))
      paste0(
        "Available data files:\n",
        paste(sprintf("• %s_%s.csv  -> get_data(type = \"%s\", data_name = \"%s\")",
                      df$type, df$data_name, df$type, df$data_name), collapse = "\n")
      )
    }
  }
  
  # Check main data file existence
  if (!data_file %in% files) {
    stop(paste0(
      "❌ '", data_file, "' not found in ", normalizePath(path, mustWork = FALSE), ".\n\n",
      suggestion_text
    ), call. = FALSE)
  }
  
  # Load data
  df <- read_csv(file.path(path, data_file), show_col_types = FALSE)

  if (identical(type, "tree") && identical(data_name, "senescence") &&
      "percent_senesced" %in% names(df) && !"remaining_green" %in% names(df)) {
    df <- df %>%
      dplyr::mutate(remaining_green = 100 - .data$percent_senesced)
  }
  
  # Optionally attach meta
  if (with_meta) {
    if (type == "tree"){
      df <- df |>
        mutate(tree_id = as.character(tree_id)) |> 
        left_join(get_meta("tree"), by="tree_id")
    } else {
      df <- df |>
        left_join(get_meta("box"), by="boxlabel")
    }
  }
  
  # Attach extreme event information (only if date column exists)
  if ("date" %in% names(df)) {
    buffer_days <- 14  # 3 weeks (to cover measurement measured after the extreme event)
    
    df <- df |> 
      mutate(
        extreme_event =
          dplyr::between(
            .data$date,
            as.Date("2025-06-20"),
            as.Date("2025-07-03") + buffer_days
          ) |
          dplyr::between(
            .data$date,
            as.Date("2025-08-12"),
            as.Date("2025-08-20") + buffer_days
          ),
        extreme_event = if_else(extreme_event, "yes", "no")
      )
    
    # Attach date_num 
    df <- df |> 
      mutate(date_num = as.numeric(.data$date))
  }
  
  # Attach SWC data (measured or imputed)
  if (type == "tree" && "date" %in% names(df)) {
    if (swc_source == "measured") {
      # Use closest measured SWC within 7-day window (original behavior)
      df_swc <- get_data("box", "soilwater", swc_source = "measured")
      
      df <- df %>%
        group_by(boxlabel) %>%
        mutate(
          swc = sapply(date, function(d) {
            # restrict to current box
            box <- unique(boxlabel)
            idx_box <- which(df_swc$boxlabel == box)
            if (length(idx_box) == 0) return(NA_real_)
            
            swc_box <- df_swc[idx_box, ]
            swc_box <- swc_box[order(swc_box$date), ]
            
            d_minus7 <- d - 7
            d_plus7  <- d + 7
            
            # 1) closest SWC within -7 days (past week, including same day)
            idx1 <- which(swc_box$date <= d & swc_box$date >= d_minus7)
            if (length(idx1) > 0) {
              # latest date in that interval
              i <- idx1[which.max(swc_box$date[idx1])]
              return(swc_box$swc[i])
            }
            
            # 2) if none, closest within +7 days (next week)
            idx2 <- which(swc_box$date > d & swc_box$date <= d_plus7)
            if (length(idx2) > 0) {
              # earliest date after d
              i <- idx2[which.min(swc_box$date[idx2])]
              return(swc_box$swc[i])
            }
            
            # 3) if still none, closest back in time (unbounded, as before)
            idx3 <- which(swc_box$date <= d)
            if (length(idx3) > 0) {
              i <- idx3[which.max(swc_box$date[idx3])]
              return(swc_box$swc[i])
            }
            
            # if there is no SWC at all before/after
            NA_real_
          })
        ) %>%
        ungroup()
    } else if (swc_source == "imputed_gam") {
      swc_candidates <- unique(c(
        file.path(path, "box_soilwater_daily_gam_agnostic.csv"),
        alinv_data_path("swc_interpolation", "box_soilwater_daily_gam_agnostic.csv"),
        .resolve_path("data/interim/box_soilwater_daily_gam_agnostic.csv")
      ))
      swc_file <- swc_candidates[file.exists(swc_candidates)][1]

      # Use imputed daily SWC (exact date match on measurement date)
      df_swc_daily <- tryCatch(
        read_csv(swc_file, show_col_types = FALSE),
        error = function(e) {
          warning(
            "Could not load imputed SWC; falling back to measured SWC. Error: ",
            conditionMessage(e),
            call. = FALSE
          )
          NULL
        }
      )

      if (is.null(df_swc_daily)) {
        # Fallback to measured SWC if imputed is not available
        return(get_data(type = type, data_name = data_name, with_meta = with_meta, 
                       path = path, swc_source = "measured"))
      }
      
      # Join imputed SWC on exact date match
      df <- df %>%
        left_join(
          df_swc_daily %>% 
            select(boxlabel, date, swc_hat) %>%
            rename(swc = swc_hat),
          by = c("boxlabel", "date")
        )
    }
  }
  
  message("✅ Loaded ", data_file, if (with_meta) paste0(" with ", meta_file) else "", 
          " (SWC source: ", swc_source, ")")
  return(df)
}


get_data_var_grid <- function(){
  tibble::tribble(
    ~type,   ~data,          ~resp_var,
    "tree",  "chlorophyll",  "chl",
    "tree",  "condition",    "condition",
    "tree",  "growth",       "height",
    "tree",  "growth",       "diameter",
    "tree",  "growth",       "volume",
    "tree",  "growth",       "height_inc_t0",
    "tree",  "growth",       "diameter_inc_t0",
    "tree",  "growth",       "volume_inc_t0",
    "tree",  "growth",       "height_inc_t0_rel",
    "tree",  "growth",       "diameter_inc_t0_rel",
    "tree",  "growth",       "volume_inc_t0_rel",
    "tree",  "growth",       "height_inc_phase_abs",
    "tree",  "growth",       "diameter_inc_phase_abs",
    "tree",  "growth",       "volume_inc_phase_abs",
    "tree",  "growth",       "height_inc_phase_rel",
    "tree",  "growth",       "diameter_inc_phase_rel",
    "tree",  "growth",       "volume_inc_phase_rel",
    "tree",  "phenology",    "stage",
    "tree",  "senescence",   "remaining_green",
    "tree",  "senescence",   "chlavg",
    "tree",  "quantum_yield","qy",
    "box",   "respiration",  "co2",
    "box",   "soilwater",    "swc"
  )
}
