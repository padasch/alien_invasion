
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

`%||%` <- function(x, y) if (is.null(x)) y else x

alinv_config_path <- file.path(.alinv_project_root(), "scripts", "config.R")
if (file.exists(alinv_config_path)) {
  source(alinv_config_path, local = FALSE)
} else {
  stop("Missing analysis config: ", alinv_config_path, call. = FALSE)
}

alinv_compute_volume_proxy <- function(diameter_mm, height_cm, species, method = NULL) {
  method <- alinv_volume_proxy_method(method)
  switch(
    method,
    "annighofer" = alinv_compute_volume_proxy_annighofer(diameter_mm, height_cm, species),
    "zianis_2005" = alinv_compute_volume_proxy_zianis(diameter_mm, height_cm, species),
    stop(
      "Unsupported volume-proxy method: ", method,
      ". Supported methods: ", paste(alinv_volume_proxy_methods(), collapse = ", "),
      call. = FALSE
    )
  )
}

alinv_compute_volume_proxy_annighofer <- function(diameter_mm, height_cm, species) {
  species_chr <- as.character(species)
  out <- rep(NA_real_, length(species_chr))

  species_levels <- unique(species_chr[!is.na(species_chr)])
  for (sp in species_levels) {
    idx <- which(species_chr == sp)
    spec_i <- tryCatch(
      alinv_volume_allometry_spec(sp),
      error = function(e) NULL
    )
    if (is.null(spec_i) || !nrow(spec_i)) {
      out[idx] <- NA_real_
      next
    }
    out[idx] <- spec_i$beta1[[1]] * (diameter_mm[idx]^2 * height_cm[idx])^spec_i$beta2[[1]]
  }

  out
}

alinv_compute_volume_proxy_zianis <- function(diameter_mm, height_cm, species) {
  species_chr <- as.character(species)
  out <- rep(NA_real_, length(species_chr))

  d_cm <- as.numeric(diameter_mm) / 10
  h_m <- as.numeric(height_cm) / 100

  species_levels <- unique(species_chr[!is.na(species_chr)])
  for (sp in species_levels) {
    idx <- which(species_chr == sp)
    spec_i <- tryCatch(
      alinv_volume_zianis_spec(sp),
      error = function(e) NULL
    )
    if (is.null(spec_i) || !nrow(spec_i)) {
      out[idx] <- NA_real_
      next
    }

    a <- spec_i$a[[1]]
    b <- spec_i$b[[1]]
    c <- spec_i$c[[1]]
    d <- spec_i$d[[1]]
    e <- spec_i$e[[1]]
    f <- spec_i$f[[1]]
    formula_type <- as.character(spec_i$formula_type[[1]])

    out[idx] <- switch(
      formula_type,
      "polynomial" = a + b * d_cm[idx] * h_m[idx]^2 + c * d_cm[idx]^3,
      "poly6" = a + b * d_cm[idx] + c * d_cm[idx]^2 + d * d_cm[idx]^3 + e * h_m[idx] + f * d_cm[idx]^2 * h_m[idx],
      "d2h" = a + b * d_cm[idx]^2 * h_m[idx],
      "power_product" = (d_cm[idx]^a) * (h_m[idx]^b) * exp(c),
      "log10_poly" = a * 10^(b * log10(d_cm[idx]) + c * log10(d_cm[idx])^2 + d * log10(h_m[idx]) + e * log10(h_m[idx])^2),
      stop(
        "Unsupported Zianis formula type: ", formula_type,
        ". Supported types: polynomial, poly6, d2h, power_product, log10_poly",
        call. = FALSE
      )
    )
  }

  out
}

alinv_volume_proxy_methods <- function() {
  ALINV_VOLUME_PROXY_METHODS
}

alinv_volume_proxy_method <- function(method = NULL) {
  if (is.null(method) || !nzchar(as.character(method))) {
    method <- Sys.getenv("ALINV_VOLUME_PROXY_METHOD", unset = ALINV_VOLUME_PROXY_DEFAULT)
  } else {
    method <- as.character(method)
  }

  method <- trimws(tolower(method))
  if (identical(method, "zianis")) {
    method <- "zianis_2005"
  }
  if (!method %in% alinv_volume_proxy_methods()) {
    stop(
      "Unsupported volume-proxy method: ", method,
      ". Supported methods: ", paste(alinv_volume_proxy_methods(), collapse = ", "),
      call. = FALSE
    )
  }
  method
}

alinv_volume_proxy_version <- function(method = NULL) {
  ALINV_VOLUME_PROXY_VERSION[[alinv_volume_proxy_method(method)]]
}

alinv_volume_allometry_range_report <- function(df,
                                                species_col = "species",
                                                diameter_col = "diameter",
                                                height_col = "height") {
  if (!all(c(species_col, diameter_col, height_col) %in% names(df))) {
    stop("Range report requires species, diameter, and height columns.", call. = FALSE)
  }

  report <- df %>%
    dplyr::mutate(
      rcd2h = .data[[diameter_col]]^2 * .data[[height_col]]
    ) %>%
    dplyr::group_by(.data[[species_col]]) %>%
    dplyr::summarise(
      n = dplyr::n(),
      min_rcd2h = min(.data$rcd2h, na.rm = TRUE),
      max_rcd2h = max(.data$rcd2h, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::rename(repo_species = 1L) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      source_species_key = alinv_volume_allometry_spec(.data$repo_species)$source_species_key[[1]],
      source_species_label = alinv_volume_allometry_spec(.data$repo_species)$source_species_label[[1]],
      rcd2h_min_calibration = alinv_volume_allometry_spec(.data$repo_species)$rcd2h_min[[1]],
      rcd2h_max_calibration = alinv_volume_allometry_spec(.data$repo_species)$rcd2h_max[[1]],
      within_calibration_range = .data$min_rcd2h >= .data$rcd2h_min_calibration &&
        .data$max_rcd2h <= .data$rcd2h_max_calibration
    ) %>%
    dplyr::ungroup()

  tibble::as_tibble(report)
}

alinv_volume_zianis_catalog <- function() {
  ALINV_VOLUME_ZIANIS$specs
}

alinv_volume_zianis_spec <- function(repo_species) {
  repo_species <- as.character(repo_species)
  specs <- alinv_volume_zianis_catalog()

  out <- specs[specs$repo_species == repo_species, , drop = FALSE]
  if (nrow(out) > 1L) {
    out <- out[1, , drop = FALSE]
  }

  if (!nrow(out)) {
    stop("No configured Zianis volume proxy for repo species: ", repo_species, call. = FALSE)
  }

  tibble::as_tibble(out)
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

ALINV_SOIL_LABELS <- alinv_level_labels("soiltype")

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
  replacement_dates = seq(as.Date("2025-07-13"), as.Date("2025-07-15"), by = "day")
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
      )
    ) %>%
    select(
      date,
      precip_mm,
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

  precip_daily %>%
    select(date, precip_mm)
}

# Load the shared site-level daily climate product generated by
# scripts/data_preparation/clean_sensor_data.R.
get_climate <- function(
  variables = NULL,
  path = "./data/interim"
) {
  path <- .resolve_path(path)
  climate_file <- file.path(path, "site_climate_daily.csv")

  if (!file.exists(climate_file)) {
    stop(
      paste0(
        "'site_climate_daily.csv' not found in ",
        normalizePath(path, winslash = "/", mustWork = FALSE),
        ". Run source(\"scripts/data_preparation/clean_sensor_data.R\") to create it."
      ),
      call. = FALSE
    )
  }

  climate <- readr::read_csv(climate_file, show_col_types = FALSE) %>%
    dplyr::mutate(date = as.Date(.data$date)) %>%
    dplyr::arrange(.data$date)

  if (!is.null(variables)) {
    variables <- unique(as.character(variables))
    unknown <- setdiff(variables, names(climate))

    if (length(unknown)) {
      stop(
        paste0(
          "Unknown climate variable(s): ", paste(unknown, collapse = ", "),
          ". Available variables: ",
          paste(setdiff(names(climate), "date"), collapse = ", "),
          "."
        ),
        call. = FALSE
      )
    }

    climate <- climate %>%
      dplyr::select("date", dplyr::all_of(setdiff(variables, "date")))
  }

  message(
    "Loaded site_climate_daily.csv (",
    nrow(climate),
    " daily records; ",
    as.character(min(climate$date, na.rm = TRUE)),
    " to ",
    as.character(max(climate$date, na.rm = TRUE)),
    ")"
  )

  climate
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

get_data <- function(type = c("tree", "box"), data_name, with_meta = TRUE, path = "./data/interim") {
  type <- match.arg(type)
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
    buffer_days <- 14  # Extend event attribution to measurements made within two weeks.
    
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
            as.Date("2025-08-05"),
            as.Date("2025-08-19") + buffer_days
          ),
        extreme_event = if_else(extreme_event, "yes", "no")
      )
    
    # Attach date_num 
    df <- df |> 
      mutate(date_num = as.numeric(.data$date))
  }
  
  # Attach measured SWC data.
  if (type == "tree" && "date" %in% names(df)) {
    # Use closest measured SWC within the established matching window.
    df_swc <- get_data("box", "soilwater")
      
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
  }
  
  message("✅ Loaded ", data_file, if (with_meta) paste0(" with ", meta_file) else "")
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
