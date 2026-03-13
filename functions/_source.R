
# Packages ------------------------------------------------------------------------------------

library(tidyverse)
library(ggplot2)
if (requireNamespace("ggalluvial", quietly = TRUE)) {
  library(ggalluvial)
}
library(readxl)
library(lubridate)
library(janitor)
library(ggpubr)
library(ggpattern)
library(rlang)
library(dplyr)
library(conflicted)



# Conflicts -----------------------------------------------------------------------------------

conflicts_prefer(
  dplyr::select,
  dplyr::filter,
  dplyr::mutate
)

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
  if (grepl("^(/|[A-Za-z]:[\\\\/])", path)) return(path)
  if (dir.exists(path) || file.exists(path)) return(path)
  file.path(.alinv_project_root(), path)
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
  if (type == "tree") {
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
      # Use imputed daily SWC (exact date match on measurement date)
      df_swc_daily <- tryCatch(
        read_csv(file.path(path, "box_soilwater_daily_gam_agnostic.csv"), show_col_types = FALSE),
        error = function(e) {
          warning("Could not load imputed SWC; falling back to measured SWC. Error: ", 
                  conditionMessage(e), call. = FALSE)
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
    # "tree",  "growth",       "xxx_inc_t0",
    # "tree",  "growth",       "xxx_inc_t0_rel",
    "tree",  "phenology",    "stage",
    "tree",  "senescence",   "percent_senesced",
    "tree",  "senescence",   "chlavg",
    "tree",  "quantum_yield","qy",
    "box",   "respiration",  "co2",
    "box",   "soilwater",    "swc"
  )
}
