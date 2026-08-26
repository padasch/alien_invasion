## Measurement-schedule data and plotting helpers.

# -------------------------------------------------------------------
# Packages
# -------------------------------------------------------------------
library(dplyr)
library(purrr)
library(ggplot2)
library(tibble)

# -------------------------------------------------------------------
# Measurement schedule settings
# -------------------------------------------------------------------

MEASUREMENT_SCHEDULE_START_DATE <- as.Date("2025-04-01")
MEASUREMENT_SCHEDULE_DATA_PATH <- "./data/interim"

MEASUREMENT_SCHEDULE_LABELS <- c(
  soil_water = "Soil water content",
  spring_phenology = "Spring phenology",
  height = "Height",
  diameter = "Diameter",
  chlorophyll_merged = "Chlorophyll content",
  chlorophyll_growing_season = "Chlorophyll (growing season)",
  quantum_yield = "Quantum yield",
  tree_condition = "Vitality",
  senescence_visual = "Senescence",
  senescence_visual_and_chlorophyll = "Senescence (visual and chlorophyll)",
  biomass = "Biomass",
  heatwave = "Heatwave"
)

MEASUREMENT_SCHEDULE_ROW_KEYS_MERGED <- c(
  "spring_phenology",
  "height",
  "diameter",
  "chlorophyll_merged",
  "quantum_yield",
  "tree_condition",
  "senescence_visual",
  "biomass",
  "soil_water",
  "heatwave"
)

MEASUREMENT_SCHEDULE_ROW_KEYS_SPLIT <- c(
  "spring_phenology",
  "height",
  "diameter",
  "chlorophyll_growing_season",
  "quantum_yield",
  "tree_condition",
  "senescence_visual_and_chlorophyll",
  "biomass",
  "soil_water",
  "heatwave"
)

MEASUREMENT_SCHEDULE_DATA_MAP <- tibble::tribble(
  ~type  , ~data           , ~resp_var         , ~merged_label_key    , ~split_label_key                    ,
  "box"  , "soilwater"     , "swc"             , "soil_water"         , "soil_water"                        ,
  "tree" , "phenology"     , "stage"           , "spring_phenology"   , "spring_phenology"                  ,
  "tree" , "growth"        , "height"          , "height"             , "height"                            ,
  "tree" , "growth"        , "diameter"        , "diameter"           , "diameter"                          ,
  "tree" , "chlorophyll"   , "chl"             , "chlorophyll_merged" , "chlorophyll_growing_season"        ,
  "tree" , "quantum_yield" , "qy"              , "quantum_yield"      , "quantum_yield"                     ,
  "tree" , "condition"     , "condition"       , "tree_condition"     , "tree_condition"                    ,
  "tree" , "senescence"    , "remaining_green" , "senescence_visual"  , "senescence_visual_and_chlorophyll" ,
  "tree" , "senescence"    , "chlavg"          , "chlorophyll_merged" , "senescence_visual_and_chlorophyll"
)

MEASUREMENT_SCHEDULE_GRID <- MEASUREMENT_SCHEDULE_DATA_MAP[c(
  "type",
  "data",
  "resp_var"
)]

MEASUREMENT_SCHEDULE_SOIL_WATER_COLOR <- "#1f77b4"
MEASUREMENT_SCHEDULE_HEATWAVE_COLOR <- "#8b1a1a"
MEASUREMENT_SCHEDULE_GRIDLINE_COLOR <- "grey88"
MEASUREMENT_SCHEDULE_COLORS <- c(
  soil_water = MEASUREMENT_SCHEDULE_SOIL_WATER_COLOR,
  spring_phenology = "#aec7e8",
  height = "#ff7f0e",
  diameter = "#ffbb78",
  chlorophyll_merged = "#2ca02c",
  chlorophyll_growing_season = "#2ca02c",
  quantum_yield = "#98df8a",
  tree_condition = "#d62728",
  senescence_visual = "#ff9896",
  senescence_visual_and_chlorophyll = "#ff9896",
  biomass = "#9467bd",
  heatwave = MEASUREMENT_SCHEDULE_HEATWAVE_COLOR
)

# Biomass has no date column; these dates come from workbook Read me extirpation notes:
# B2 on 2025-12-10 and B3 on 2025-12-11.
MEASUREMENT_SCHEDULE_BIOMASS_DATES <- as.Date(c("2025-12-10", "2025-12-11"))

MEASUREMENT_SCHEDULE_HEATWAVE_WINDOWS <- tibble::tribble(
  ~start                , ~end                  ,
  as.Date("2025-06-20") , as.Date("2025-07-02") ,
  as.Date("2025-08-05") , as.Date("2025-08-19")
)
MEASUREMENT_SCHEDULE_HEATWAVE_DATES <- unique(unlist(Map(
  seq.Date,
  MEASUREMENT_SCHEDULE_HEATWAVE_WINDOWS$start,
  MEASUREMENT_SCHEDULE_HEATWAVE_WINDOWS$end,
  by = "day"
)))

MEASUREMENT_SCHEDULE_TICK_DAYS <- c(1L, 15L)
MEASUREMENT_SCHEDULE_SHOW_END_DATE_TICK <- TRUE
MEASUREMENT_SCHEDULE_DATE_LABELS <- "%d %b"
MEASUREMENT_SCHEDULE_SHOW_SWC_LINES <- FALSE
MEASUREMENT_SCHEDULE_SHOW_HEATWAVE <- TRUE
MEASUREMENT_SCHEDULE_HEATWAVE_CONTINUOUS <- TRUE
MEASUREMENT_SCHEDULE_MERGE_SENESCENCE_CHLOROPHYLL <- FALSE
MEASUREMENT_SCHEDULE_CLUSTER_WINDOW_DAYS <- 3
MEASUREMENT_SCHEDULE_SQUARE_SIZE <- 2.8
MEASUREMENT_SCHEDULE_SQUARE_ALPHA <- 0.9
MEASUREMENT_SCHEDULE_GRIDLINE_WIDTH <- 0.3
MEASUREMENT_SCHEDULE_SWC_LINE_WIDTH <- 0.25
MEASUREMENT_SCHEDULE_SWC_LINE_ALPHA <- 1
MEASUREMENT_SCHEDULE_HEATWAVE_LINE_WIDTH <- 2.5
MEASUREMENT_SCHEDULE_BASE_SIZE <- 12
MEASUREMENT_SCHEDULE_X_TEXT_ANGLE <- 45
MEASUREMENT_SCHEDULE_X_TEXT_HJUST <- 1
MEASUREMENT_SCHEDULE_Y_EXPAND <- c(0.05, 0.05)

# -------------------------------------------------------------------
# 1) Table of datasets and response variables
# -------------------------------------------------------------------

get_sem_dataset_grid <- get_data_var_grid

get_measurement_schedule_grid <- function() {
  MEASUREMENT_SCHEDULE_GRID
}

measurement_schedule_date_breaks <- function(
  start_date,
  end_date,
  tick_days = MEASUREMENT_SCHEDULE_TICK_DAYS,
  include_end_date = MEASUREMENT_SCHEDULE_SHOW_END_DATE_TICK
) {
  first_month <- as.Date(format(start_date, "%Y-%m-01"))
  last_month <- as.Date(format(end_date, "%Y-%m-01"))
  month_starts <- seq.Date(first_month, last_month, by = "month")
  breaks <- unlist(
    lapply(tick_days, function(day) {
      as.Date(sprintf("%s-%02d", format(month_starts, "%Y-%m"), day))
    })
  )
  breaks <- sort(unique(as.Date(breaks, origin = "1970-01-01")))
  breaks <- breaks[breaks >= start_date & breaks <= end_date]
  if (isTRUE(include_end_date)) {
    breaks <- sort(unique(c(breaks, as.Date(end_date))))
  }
  breaks
}

measurement_schedule_resolve_path <- function(path) {
  if (exists(".resolve_path", mode = "function")) {
    return(.resolve_path(path))
  }
  path
}

measurement_schedule_read_data <- function(
  type,
  data_name,
  path = MEASUREMENT_SCHEDULE_DATA_PATH
) {
  path <- measurement_schedule_resolve_path(path)
  data_file <- file.path(path, sprintf("%s_%s.csv", type, data_name))

  if (!file.exists(data_file)) {
    stop(
      "Measurement schedule data file not found: ",
      normalizePath(data_file, mustWork = FALSE),
      call. = FALSE
    )
  }

  df <- readr::read_csv(data_file, show_col_types = FALSE)

  if (
    identical(type, "tree") &&
      identical(data_name, "senescence") &&
      "percent_senesced" %in% names(df) &&
      !"remaining_green" %in% names(df)
  ) {
    df <- df %>%
      dplyr::mutate(remaining_green = 100 - .data$percent_senesced)
  }

  df
}

measurement_schedule_cluster_dates <- function(
  meas_tbl,
  max_gap_days = MEASUREMENT_SCHEDULE_CLUSTER_WINDOW_DAYS
) {
  if (is.null(max_gap_days) || is.na(max_gap_days) || max_gap_days <= 0) {
    return(meas_tbl)
  }

  meas_tbl %>%
    mutate(date = as.Date(date)) %>%
    arrange(schedule_label, date) %>%
    group_by(schedule_label) %>%
    mutate(
      previous_date = dplyr::lag(date),
      cluster_id = cumsum(
        is.na(previous_date) | as.integer(date - previous_date) > max_gap_days
      )
    ) %>%
    group_by(schedule_label, cluster_id) %>%
    summarise(
      type = dplyr::first(type),
      data = dplyr::first(data),
      resp_var = dplyr::first(resp_var),
      date = as.Date(
        floor(stats::median(as.numeric(date)) + 0.5),
        origin = "1970-01-01"
      ),
      .groups = "drop"
    ) %>%
    select(type, data, resp_var, date, schedule_label)
}

# -------------------------------------------------------------------
# 2) Build table of unique measurement dates per dataset/resp_var
# -------------------------------------------------------------------

build_measurement_table <- function(
  grid = get_measurement_schedule_grid(),
  path = MEASUREMENT_SCHEDULE_DATA_PATH
) {
  grid <- grid %>%
    dplyr::distinct(type, data, resp_var)

  dataset_grid <- grid %>%
    dplyr::distinct(type, data)

  purrr::pmap_dfr(
    dataset_grid,
    function(type, data) {
      dataset_type <- type
      dataset_name <- data
      df <- measurement_schedule_read_data(
        type = dataset_type,
        data_name = dataset_name,
        path = path
      )

      if (!"date" %in% names(df)) {
        stop(
          "Measurement schedule data has no date column: ",
          dataset_type,
          "_",
          dataset_name,
          ".csv",
          call. = FALSE
        )
      }

      resp_vars <- grid %>%
        dplyr::filter(
          .data$type == dataset_type,
          .data$data == dataset_name
        ) %>%
        dplyr::pull(resp_var)

      purrr::map_dfr(
        resp_vars,
        function(resp_var) {
          if (!resp_var %in% names(df)) {
            warning(
              "Skipping missing measurement schedule column: ",
              dataset_type,
              "_",
              dataset_name,
              "$",
              resp_var,
              call. = FALSE
            )
            return(tibble::tibble(
              type = character(),
              data = character(),
              resp_var = character(),
              date = as.Date(character())
            ))
          }

          df %>%
            dplyr::mutate(date = as.Date(date)) %>%
            dplyr::filter(!is.na(.data[[resp_var]])) %>%
            dplyr::distinct(date) %>%
            dplyr::mutate(
              type = dataset_type,
              data = dataset_name,
              resp_var = resp_var
            ) %>%
            dplyr::select(type, data, resp_var, date)
        }
      )
    }
  )
}

# -------------------------------------------------------------------
# 3) Gantt-like plot of measurement schedule
#    - rows: data: resp_var
#    - x: dates (from April 2025)
#    - coloured by series
#    - right-side labels for each series
#    - optional SWC dates shown as vertical "reference bars" across all rows
# -------------------------------------------------------------------

plot_measurement_schedule <- function(
  meas_tbl,
  show_swc_lines = MEASUREMENT_SCHEDULE_SHOW_SWC_LINES,
  show_heatwave = MEASUREMENT_SCHEDULE_SHOW_HEATWAVE,
  heatwave_continuous = MEASUREMENT_SCHEDULE_HEATWAVE_CONTINUOUS,
  merge_senescence_chlorophyll = MEASUREMENT_SCHEDULE_MERGE_SENESCENCE_CHLOROPHYLL,
  biomass_dates = MEASUREMENT_SCHEDULE_BIOMASS_DATES,
  cluster_within_days = MEASUREMENT_SCHEDULE_CLUSTER_WINDOW_DAYS,
  heatwave_color = MEASUREMENT_SCHEDULE_HEATWAVE_COLOR,
  square_size = MEASUREMENT_SCHEDULE_SQUARE_SIZE,
  square_alpha = MEASUREMENT_SCHEDULE_SQUARE_ALPHA
) {
  library(dplyr)
  library(ggplot2)

  # 1) restrict to dates from April 2025 onward
  start_date <- MEASUREMENT_SCHEDULE_START_DATE
  row_keys <- if (isTRUE(merge_senescence_chlorophyll)) {
    MEASUREMENT_SCHEDULE_ROW_KEYS_MERGED
  } else {
    MEASUREMENT_SCHEDULE_ROW_KEYS_SPLIT
  }
  if (!isTRUE(show_heatwave)) {
    row_keys <- setdiff(row_keys, "heatwave")
  }
  row_order <- unname(MEASUREMENT_SCHEDULE_LABELS[row_keys])

  label_key_col <- if (isTRUE(merge_senescence_chlorophyll)) {
    "merged_label_key"
  } else {
    "split_label_key"
  }
  schedule_map <- MEASUREMENT_SCHEDULE_DATA_MAP
  schedule_map$schedule_label <- unname(MEASUREMENT_SCHEDULE_LABELS[schedule_map[[
    label_key_col
  ]]])
  schedule_map <- schedule_map[c("type", "data", "resp_var", "schedule_label")]

  meas_tbl <- meas_tbl %>%
    filter(date >= start_date) %>%
    inner_join(schedule_map, by = c("type", "data", "resp_var")) %>%
    distinct(schedule_label, date, .keep_all = TRUE)

  biomass_dates <- as.Date(biomass_dates)
  biomass_dates <- biomass_dates[!is.na(biomass_dates)]
  if (length(biomass_dates) > 0) {
    meas_tbl <- dplyr::bind_rows(
      meas_tbl,
      tibble::tibble(
        type = "tree",
        data = "biomass",
        resp_var = "biomass",
        date = unique(biomass_dates),
        schedule_label = MEASUREMENT_SCHEDULE_LABELS[["biomass"]]
      )
    ) %>%
      distinct(schedule_label, date, .keep_all = TRUE)
  }

  if (
    isTRUE(show_heatwave) &&
      !isTRUE(heatwave_continuous) &&
      length(MEASUREMENT_SCHEDULE_HEATWAVE_DATES) > 0
  ) {
    meas_tbl <- dplyr::bind_rows(
      meas_tbl,
      tibble::tibble(
        type = "event",
        data = "heatwave",
        resp_var = "heatwave",
        date = as.Date(
          MEASUREMENT_SCHEDULE_HEATWAVE_DATES,
          origin = "1970-01-01"
        ),
        schedule_label = MEASUREMENT_SCHEDULE_LABELS[["heatwave"]]
      )
    ) %>%
      distinct(schedule_label, date, .keep_all = TRUE)
  }

  if (nrow(meas_tbl) == 0) {
    stop("No measurement schedule rows remain after filtering.", call. = FALSE)
  }

  meas_tbl <- measurement_schedule_cluster_dates(
    meas_tbl,
    max_gap_days = cluster_within_days
  )

  meas_tbl <- meas_tbl %>%
    mutate(
      schedule_label = factor(schedule_label, levels = rev(row_order)),
      y_id = as.numeric(schedule_label)
    )

  # y breaks / labels (used for axis)
  y_breaks <- seq_along(row_order)
  y_labels <- rev(row_order)
  heatwave_y <- match(MEASUREMENT_SCHEDULE_LABELS[["heatwave"]], y_labels)
  heatwave_lines <- MEASUREMENT_SCHEDULE_HEATWAVE_WINDOWS %>%
    mutate(
      y = heatwave_y,
      yend = heatwave_y
    )

  # 4) tab20 palette, with soil water content fixed to tab20 blue
  pal <- unname(MEASUREMENT_SCHEDULE_COLORS[row_keys])
  names(pal) <- row_order
  pal[[MEASUREMENT_SCHEDULE_LABELS[[
    "soil_water"
  ]]]] <- MEASUREMENT_SCHEDULE_SOIL_WATER_COLOR
  pal[[MEASUREMENT_SCHEDULE_LABELS[["heatwave"]]]] <- heatwave_color

  swc_col <- pal[[MEASUREMENT_SCHEDULE_LABELS[["soil_water"]]]]

  # 5) x range
  x_dates <- meas_tbl$date
  if (isTRUE(show_heatwave) && nrow(heatwave_lines) > 0) {
    x_dates <- c(x_dates, heatwave_lines$end)
  }
  x_max <- max(x_dates, na.rm = TRUE)
  x_breaks <- measurement_schedule_date_breaks(start_date, x_max)
  x_break_lines <- tibble::tibble(date = x_breaks)

  # 6) SWC vertical reference lines across all rows, when requested
  y_min <- min(y_breaks) - 0.5
  y_max <- max(y_breaks) + 0.5

  swc_lines <- meas_tbl %>%
    filter(
      type == "box",
      data == "soilwater",
      resp_var == "swc"
    ) %>%
    distinct(date) %>%
    mutate(
      ymin = y_min,
      ymax = y_max
    )

  p <- ggplot(meas_tbl, aes(x = date, y = y_id, colour = schedule_label)) +
    geom_vline(
      data = x_break_lines,
      aes(xintercept = date),
      colour = MEASUREMENT_SCHEDULE_GRIDLINE_COLOR,
      linewidth = MEASUREMENT_SCHEDULE_GRIDLINE_WIDTH
    )

  if (isTRUE(show_swc_lines) && nrow(swc_lines) > 0) {
    p <- p +
      # SWC reference bars in SWC colour (background)
      geom_segment(
        data = swc_lines,
        aes(
          x = date,
          xend = date,
          y = ymin,
          yend = ymax
        ),
        inherit.aes = FALSE,
        colour = swc_col,
        linewidth = MEASUREMENT_SCHEDULE_SWC_LINE_WIDTH,
        alpha = MEASUREMENT_SCHEDULE_SWC_LINE_ALPHA
        # linewidth = 1.5,
        # alpha = 0.15
      )
  }

  if (
    isTRUE(show_heatwave) &&
      isTRUE(heatwave_continuous) &&
      nrow(heatwave_lines) > 0
  ) {
    p <- p +
      geom_segment(
        data = heatwave_lines,
        aes(
          x = start,
          xend = end,
          y = y,
          yend = yend
        ),
        inherit.aes = FALSE,
        colour = heatwave_color,
        linewidth = MEASUREMENT_SCHEDULE_HEATWAVE_LINE_WIDTH,
        lineend = "round"
      )
  }

  p +
    # measurement dates per series
    geom_point(
      shape = 15,
      size = square_size,
      alpha = square_alpha
    ) +
    scale_y_continuous(
      name = NULL,
      breaks = y_breaks,
      labels = y_labels,
      expand = expansion(mult = MEASUREMENT_SCHEDULE_Y_EXPAND)
    ) +
    scale_x_date(
      name = NULL,
      limits = c(start_date, x_max),
      breaks = x_breaks,
      date_labels = MEASUREMENT_SCHEDULE_DATE_LABELS
    ) +
    scale_colour_manual(values = pal) +
    labs(
      colour = NULL,
      title = NULL
    ) +
    theme_minimal(base_size = MEASUREMENT_SCHEDULE_BASE_SIZE) +
    theme(
      axis.text.y = element_text(face = "bold"), # left labels, bold
      panel.grid.major.x = element_blank(), # no major x grid
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.text.x = element_text(
        angle = MEASUREMENT_SCHEDULE_X_TEXT_ANGLE,
        hjust = MEASUREMENT_SCHEDULE_X_TEXT_HJUST
      ),
      legend.position = "none"
    )
}

# -------------------------------------------------------------------
# Example usage (if you want to test interactively)
# -------------------------------------------------------------------
# grid <- get_measurement_schedule_grid()
# meas_tbl <- build_measurement_table(grid)
# plot_measurement_schedule(meas_tbl)
# plot_measurement_schedule(meas_tbl, show_swc_lines = TRUE)
# plot_measurement_schedule(meas_tbl, merge_senescence_chlorophyll = FALSE)
