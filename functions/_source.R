
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
