# -------------------------------------------------------------------
# Packages
# -------------------------------------------------------------------
library(dplyr)
library(purrr)
library(ggplot2)
library(tibble)

# -------------------------------------------------------------------
# 1) Table of datasets and response variables
# -------------------------------------------------------------------

get_sem_dataset_grid <- get_data_var_grid

# -------------------------------------------------------------------
# 2) Build table of unique measurement dates per dataset/resp_var
# -------------------------------------------------------------------

build_measurement_table <- function(grid = get_sem_dataset_grid()) {
  purrr::pmap_dfr(
    grid,
    function(type, data, resp_var) {
      df <- get_data(type = type, data_name = data)
      
      df %>%
        dplyr::mutate(date = as.Date(date)) %>%
        dplyr::filter(!is.na(.data[[resp_var]])) %>%
        dplyr::distinct(date) %>%
        dplyr::mutate(
          type     = type,
          data     = data,
          resp_var = resp_var
        ) %>%
        dplyr::select(type, data, resp_var, date)
    }
  )
}

# -------------------------------------------------------------------
# 3) Gantt-like plot of measurement schedule
#    - rows: data: resp_var
#    - x: dates (from April 2025)
#    - coloured by series
#    - right-side labels for each series
#    - SWC dates shown as vertical "reference bars" across all rows
# -------------------------------------------------------------------

plot_measurement_schedule <- function(meas_tbl) {
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
  
  # 1) restrict to dates from April 2025 onward
  start_date <- as.Date("2025-04-01")
  meas_tbl <- meas_tbl %>%
    filter(date >= start_date)
  
  # 2) combined label: "data: resp_var"
  meas_tbl <- meas_tbl %>%
    mutate(combo = paste0(data, ": ", resp_var))
  
  # 3) stable order for rows
  combo_levels <- meas_tbl %>%
    distinct(type, data, resp_var, combo) %>%
    arrange(type, data, resp_var) %>%
    pull(combo)
  
  meas_tbl <- meas_tbl %>%
    mutate(
      combo = factor(combo, levels = combo_levels),
      y_id  = as.numeric(combo)   # numeric row index
    )
  
  # y breaks / labels (used for axis)
  y_breaks <- seq_along(combo_levels)
  y_labels <- combo_levels
  
  # 4) palette so we can reuse SWC colour
  n_cols <- max(3, min(length(combo_levels), 12))
  pal    <- brewer.pal(n_cols, "Paired")[seq_along(combo_levels)]
  names(pal) <- combo_levels
  
  swc_name <- "soilwater: swc"
  swc_col  <- if (swc_name %in% names(pal)) pal[[swc_name]] else "grey80"
  
  # 5) x range
  x_max <- max(meas_tbl$date, na.rm = TRUE)
  
  # 6) SWC vertical reference lines across all rows
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
  
  ggplot(meas_tbl, aes(x = date, y = y_id, colour = combo)) +
    # SWC reference bars in SWC colour (background)
    geom_segment(
      data = swc_lines,
      aes(
        x    = date,
        xend = date,
        y    = ymin,
        yend = ymax
      ),
      inherit.aes = FALSE,
      colour = swc_col,
      linewidth = 0.25,
      alpha = 1
      # linewidth = 1.5,
      # alpha = 0.15
    ) +
    # measurement "ticks" per series
    geom_segment(
      aes(
        xend = date,
        y    = y_id - 0.35,
        yend = y_id + 0.35
      ),
      linewidth = 1.5
    ) +
    scale_y_continuous(
      name   = "Data: Response variable",
      breaks = y_breaks,
      labels = y_labels,
      expand = expansion(mult = c(0.05, 0.05))
    ) +
    scale_x_date(
      name        = "Date",
      limits      = c(start_date, x_max),
      date_breaks = "2 weeks",
      date_labels = "%d-%b"
    ) +
    scale_colour_manual(values = pal) +
    labs(
      colour = NULL,
      title  = "Measurement schedule across datasets"
    ) +
    annotate("segment",
             x = as.Date("2025-06-20"), xend = as.Date("2025-07-02"),
             y = 0, yend = 0,
             size = 2.5, color = "darkorange", lineend = "round"
    ) +
    annotate("segment",
             x = as.Date("2025-08-12"), xend = as.Date("2025-08-20"),
             y = 0, yend = 0,
             size = 2.5, color = "darkorange", lineend = "round"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.y        = element_text(face = "bold"),  # left labels, bold
      panel.grid.major.x = element_blank(),              # no major x grid
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.position    = "none"
    )
}

# -------------------------------------------------------------------
# Example usage (if you want to test)
# -------------------------------------------------------------------
# grid      <- get_sem_dataset_grid()
# meas_tbl  <- build_measurement_table(grid)
# plot_measurement_schedule(meas_tbl)

