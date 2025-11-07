library(tidyverse)
library(ggplot2)
library(ggalluvial)
library(readxl)
library(lubridate)
library(janitor)
library(ggpubr)
library(ggpattern)
library(rlang)

# Main Functions
get_meta <- function(which = c("tree", "box")) {
  if (which == "tree") {
    return(
      left_join(
        read_csv("./data/interim/meta_tree.csv",show_col_types = FALSE),
        read_csv("./data/interim/meta_box.csv",show_col_types = FALSE),
        by="boxlabel"
      )
    )
  } else if (which == "box") {
    return(read_csv("./data/interim/meta_box.csv", show_col_types = FALSE))
  } else {
    stop("Invalid 'which' argument. Use 'tree' or 'box'.")
  }
}

get_data <- function(type = c("tree", "box"), data_name, with_meta = TRUE, path = "./data/interim") {
  type <- match.arg(type)
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
      df <- df %>% left_join(get_meta("tree"), by="tree_id")
    } else {
      df <- df %>% left_join(get_meta("box"), by="boxlabel")
    }
    
  }
  
  message("✅ Loaded ", data_file, if (with_meta) paste0(" with ", meta_file) else "")
  df
}
