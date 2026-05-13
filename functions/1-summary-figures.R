# Global plotting config & helpers ----

PLOT_STYLE <- list(
  levels = list(
    precipitation = c("control","drought"),
    culture       = c("mono","mixed"),
    robinia       = c("without-robinia","with-robinia"),
    species       = c("fagus","quercus","robinia"),
    soiltype      = c("inoc_robinia","inoc_beech")
  ),
  colors = list(
    precipitation = c(control = "#4F6674", drought = "indianred")
    # precipitation = c(control = "#4F6674", drought = "indianred")
  ),
  linetypes = list(
    culture = c(mono = "solid", mixed = "dotted")
  ),
  alphas = list(
    culture = c(mono = 0.8, mixed = 0.25)   # <- boxplots
  ),
  sizes = list(
    line = 1,
    ribbon_alpha = 0.15,
    error_width = 6,
    point = 1.6
  )
)

get_ylim <- function(var) {
  # Define a named list of y-limits for all variables
  lims <- list(
    diameter       = c(0, 8),
    height         = c(0, 60),
    diameter_boxplot = c(0, 12),
    height_boxplot = c(0, 100),
    qy             = c(0, 1),
    chlorophyll    = c(0, 25),
    senescence     = c(0, 100),
    condition      = c(0, 5),
    swc            = c(-2, 30),
    respiration    = c(0, 10),
    phenology_doy  = c(85, 130)
  )
  
  # Return the limit if found, otherwise NULL
  lims[[var]] %||% NULL
}

`%||%` <- function(x, y) if (is.null(x)) y else x

# -------- factor normalization (tree- vs soil-level) --------
normalize_factors_tree <- function(df) {
  df %>%
    mutate(
      precipitation = factor(precipitation, PLOT_STYLE$levels$precipitation),
      culture       = factor(culture,       PLOT_STYLE$levels$culture),
      robinia       = factor(robinia,       PLOT_STYLE$levels$robinia),
      species       = factor(species,       PLOT_STYLE$levels$species)
    )
}
normalize_factors_soil <- function(df) {
  df %>%
    mutate(
      precipitation = factor(precipitation, PLOT_STYLE$levels$precipitation),
      culture       = factor(culture,       PLOT_STYLE$levels$culture),
      robinia       = factor(robinia,       PLOT_STYLE$levels$robinia),
      soiltype      = factor(gsub("-", "_", soiltype), PLOT_STYLE$levels$soiltype)
    )
}

# -------- universal theme & scales --------
theme_common <- function() {
  theme_linedraw(base_size = 13) +
  theme(
    legend.position = "bottom",
    panel.grid.major = element_line(color = "lightgrey", linewidth = 0.3),
    panel.grid.minor = element_line(color = "lightgrey", linewidth = 0.2),
    panel.grid.major.x = element_line(),
    panel.grid.major.y = element_line(),
    panel.grid.minor.x = element_line(),
    strip.text = element_text(face = "bold"),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11),
    # panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4)  # keeps full borders in facets
  )
}
scale_precip_color <- function(name = "Precipitation") {
  scale_color_manual(values = PLOT_STYLE$colors$precipitation, name = name)
}
scale_precip_fill <- function(name = "Precipitation") {
  scale_fill_manual(values = PLOT_STYLE$colors$precipitation, name = name)
}
scale_culture_linetype <- function(name = "Culture") {
  scale_linetype_manual(values = PLOT_STYLE$linetypes$culture, name = name)
}
scale_culture_alpha <- function(name = "Culture") {
  scale_alpha_manual(values = PLOT_STYLE$alphas$culture, name = name)
}

# -------- soiltype display labels --------
SOILTYPE_LABELS <- c(
  inoc_beech   = "drier soil (beech soil)",
  inoc_robinia = "wetter soil (robinia soil)",
  `inoc-beech`   = "drier soil (beech soil)",
  `inoc-robinia` = "wetter soil (robinia soil)"
)

# -------- common facets --------
facet_tree <- function() {
  facet_grid(
    rows = vars(species),
    cols = vars(robinia),
    labeller = labeller(
      robinia = c(`without-robinia`="without robinia",
                  `with-robinia`   ="with robinia")
    )
  )
}
facet_soil <- function(show_soiltype = TRUE) {
  if (isTRUE(show_soiltype)) {
    facet_grid(
      rows = vars(soiltype),
      cols = vars(robinia),
      labeller = labeller(
        robinia  = c(`without-robinia`="without robinia",
                     `with-robinia`   ="with robinia"),
        soiltype = SOILTYPE_LABELS
      )
    )
  } else {
    facet_grid(
      . ~ robinia,
      labeller = labeller(
        robinia = c(`without-robinia`="without robinia",
                    `with-robinia`   ="with robinia")
      )
    )
  }
}

# -------- mean ± SE summarizer for time-series --------
summarize_ts <- function(df, y, keys) {
  y <- ensym(y)
  df %>%
    group_by(across(all_of(keys)), date) %>%
    summarise(
      mean = mean(!!y, na.rm = TRUE),
      sd   = sd(!!y,   na.rm = TRUE),
      n    = dplyr::n(),
      .groups = "drop"
    ) %>%
    mutate(se = sd / sqrt(pmax(n, 1)))
}

# -------- CI layers (modifiable in one place) --------
ci_layers <- function(style = c("band","errorbar")) {
  style <- match.arg(style)
  if (style == "band") {
    list(
      geom_ribbon(aes(ymin = mean - se, ymax = mean + se),
                  alpha = PLOT_STYLE$sizes$ribbon_alpha, color = NA),
      geom_line(size = PLOT_STYLE$sizes$line)
    )
  } else {
    list(
      geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                    width = PLOT_STYLE$sizes$error_width, alpha = 0.6),
      geom_point(size = PLOT_STYLE$sizes$point),
      geom_line(size = 0.8)
    )
  }
}

# ---- 1) Time-series builder (tree- or soil-level) ----
# df_in  : data.frame already normalized via normalize_factors_*()
# y      : unquoted y column to summarise (e.g., qy, inc, chl, swc, co2, ...)
# facet  : "tree" (rows=species) or "soil" (rows=soiltype)
# style  : "band" or "errorbar"
build_timeseries <- function(df_in, y, facet = c("tree","soil"),
                             style = "band",
                             ylab = NULL, title = NULL,
                             x_date_fmt = "%m/%y",
                             show_soiltype = TRUE) {
  facet <- match.arg(facet)
  keys  <- if (facet == "tree")
    c("species","robinia","precipitation","culture")
  else if (isTRUE(show_soiltype))
    c("soiltype","robinia","precipitation","culture")
  else
    c("robinia","precipitation","culture")
  
  df_sum <- summarize_ts(df_in, {{y}}, keys)
  
  p <- ggplot(
    df_sum,
    aes(x = date, y = mean,
        color = precipitation, fill = precipitation,
        linetype = culture,
        group = interaction(precipitation, culture))
  ) +
    ci_layers(style) +
    { if (facet == "tree") facet_tree() else facet_soil(show_soiltype = show_soiltype) } +
    scale_precip_color() +
    scale_precip_fill() +
    scale_culture_linetype() +
    scale_x_date(date_labels = x_date_fmt) +
    labs(x = "Date", y = ylab %||% as_name(ensym(y)), title = title) +
    theme_common()
  
  return(p)
}

# ---- 2) Boxplot builder (two boxes per culture: control vs drought) ----
# df_in : data.frame already normalized via normalize_factors_*()
# y     : unquoted y column to plot
# facet : "tree" or "soil"
build_boxplot <- function(df_in, y, facet = c("tree","soil"),
                          ylab = NULL, title = NULL,
                          show_soiltype = TRUE) {
  facet <- match.arg(facet)
  dodge_w <- 0.575
  
  p <- ggplot(
    df_in,
    aes(x = culture, y = {{y}},
        fill = precipitation, color = precipitation,
        alpha = culture)
  ) +
    geom_boxplot(
      aes(group = interaction(culture, precipitation)),
      width = 0.6,
      position = position_dodge2(width = dodge_w, preserve = "single"),
      outlier.shape = NA
    ) +
    geom_point(
      aes(group = interaction(culture, precipitation)),
      position = position_jitterdodge(jitter.width = 0.12, dodge.width = dodge_w, seed = 1),
      size = 1.2, show.legend = FALSE
    ) +
    { if (facet == "tree") facet_tree() else facet_soil(show_soiltype = show_soiltype) } +
    scale_precip_fill() +
    scale_precip_color() +
    scale_culture_alpha() +
    labs(x = "Culture", y = ylab %||% as_name(ensym(y)), title = title) +
    theme_common()
  
  return(p)
}

# ===============================
# Wrappers with drought + saving
# ===============================
library(ggplot2)
library(dplyr)
library(lubridate)

# ---- drought periods & helpers ----
DROUGHT_PERIODS <- list(c("2025-06-20","2025-07-07"),
                        c("2025-08-12","2025-08-20"))

y_limits_from_sum <- function(df_sum) {
  ymin <- min(df_sum$mean - df_sum$se, na.rm = TRUE)
  ymax <- max(df_sum$mean + df_sum$se, na.rm = TRUE)
  pad  <- 0.05 * (ymax - ymin)
  c(ymin, ymax + pad)
}

add_drought_bars_background <- function(p, df_sum, show = FALSE) {
  if (!isTRUE(show)) return(p)
  if (!"date" %in% names(df_sum) || !nrow(df_sum)) return(p)

  date_range <- range(df_sum$date, na.rm = TRUE)
  if (any(!is.finite(as.numeric(date_range)))) return(p)

  visible_periods <- purrr::keep(DROUGHT_PERIODS, function(win) {
    win_start <- as.Date(win[[1]])
    win_end <- as.Date(win[[2]])
    win_start <= date_range[2] && win_end >= date_range[1]
  })

  if (!length(visible_periods)) {
    return(p)
  }

  yl <- y_limits_from_sum(df_sum)
  ybar <- yl[1]

  for (win in visible_periods) {
    win_start <- max(as.Date(win[[1]]), date_range[1])
    win_end <- min(as.Date(win[[2]]), date_range[2])
    if (win_start > win_end) next

    p <- p + annotate(
      "segment",
      x = win_start,
      xend = win_end,
      y = ybar,
      yend = ybar,
      size = 2.5,
      color = "orange",
      lineend = "round",
      alpha = 0.9
    )
  }

  p + coord_cartesian(ylim = c(yl[1] - 0, yl[2]))
}

save_plot <- function(p, stub, soiltype = "both",
                      dir_root = "output", subdir = "timeseries",
                      device = c("pdf","png"), width = 8.5, height = 6, dpi = 300) {
  device <- match.arg(device)
  dirout <- alinv_data_path("figure_exports", subdir, paste0("soil_", gsub("-","_", soiltype)), create_dir = TRUE)
  if (!dir.exists(dirout)) dir.create(dirout, recursive = TRUE, showWarnings = FALSE)
  fpath <- file.path(dirout, paste0(stub, ".", device))
  if (device == "pdf") ggsave(fpath, p, width = width, height = height, device = cairo_pdf)
  else                ggsave(fpath, p, width = width, height = height, dpi = dpi)
  message("Saved: ", fpath)
  invisible(p)
}

# ---------------------------
# TREE-LEVEL: Growth (Height)
# ---------------------------
plot_ts_height_inc <- function(filter_soiltype = c("both", "inoc-beech", "inoc-robinia"),
                               style = "band",
                               drought_bars = FALSE,
                               save_fig = FALSE, device = "pdf") {
  
  df <- get_data("tree","growth") %>%
    { if (filter_soiltype != "both") filter(., soiltype == filter_soiltype) else . } %>%
    filter(!is.na(height)) %>%
    group_by(tree_id) %>% arrange(date, .by_group = TRUE) %>%
    mutate(inc = height - first(height)) %>% ungroup() %>%
    normalize_factors_tree()
  
  df_sum <- summarize_ts(df, inc, keys = c("species","robinia","precipitation","culture"))
  
  p <- build_timeseries(df, inc, facet = "tree", style = style,
                        ylab = "Height increment (cm)",
                        title = paste0("Height increment over time",
                                       if (filter_soiltype!="both") paste0(" (soiltype = ", filter_soiltype, ")")))
  
  yl <- get_ylim("height")
  if (!is.null(yl)) p <- p + coord_cartesian(ylim = yl)
  p <- add_drought_bars_background(p, df_sum, drought_bars)
  if (save_fig) save_plot(p, paste0("height_increment_", gsub("-","_", filter_soiltype)),
                          soiltype = filter_soiltype, device = device)
  p
}

# -----------------------------
# TREE-LEVEL: Growth (Diameter)
# -----------------------------
plot_ts_diam_inc <- function(filter_soiltype = c("both", "inoc-beech", "inoc-robinia"),
                             style = "band",
                             drought_bars = FALSE,
                             save_fig = FALSE, device = "pdf") {
  
  df <- get_data("tree","growth") %>%
    { if (filter_soiltype != "both") filter(., soiltype == filter_soiltype) else . } %>%
    filter(!is.na(diameter)) %>%
    group_by(tree_id) %>% arrange(date, .by_group = TRUE) %>%
    mutate(inc = diameter - first(diameter)) %>% ungroup() %>%
    normalize_factors_tree()
  
  df_sum <- summarize_ts(df, inc, keys = c("species","robinia","precipitation","culture"))
  
  p <- build_timeseries(df, inc, facet = "tree", style = style,
                        ylab = "Diameter increment (mm)",
                        title = paste0("Diameter increment over time",
                                       if (filter_soiltype!="both") paste0(" (soiltype = ", filter_soiltype, ")")))
  
  yl <- get_ylim("diameter")
  if (!is.null(yl)) p <- p + coord_cartesian(ylim = yl)
  p <- add_drought_bars_background(p, df_sum, drought_bars)
  if (save_fig) save_plot(p, paste0("diameter_increment_", gsub("-","_", filter_soiltype)),
                          soiltype = filter_soiltype, device = device)
  p
}

# ------------------------
# TREE-LEVEL: Quantum yield
# ------------------------
plot_ts_qy <- function(filter_soiltype = c("both", "inoc-beech", "inoc-robinia"),
                       style = "band",
                       drought_bars = FALSE,
                       save_fig = FALSE, device = "pdf") {
  
  df <- get_data("tree","quantum_yield") %>%
    { if (filter_soiltype != "both") filter(., soiltype == filter_soiltype) else . } %>%
    filter(!is.na(qy)) %>%
    normalize_factors_tree()
  
  df_sum <- summarize_ts(df, qy, keys = c("species","robinia","precipitation","culture"))
  
  p <- build_timeseries(df, qy, facet = "tree", style = style,
                        ylab = "Quantum yield (QY)",
                        title = paste0("Quantum yield over time",
                                       if (filter_soiltype!="both") paste0(" (soiltype = ", filter_soiltype, ")")))
  
  yl <- get_ylim("qy")
  if (!is.null(yl)) p <- p + coord_cartesian(ylim = yl)
  
  p <- add_drought_bars_background(p, df_sum, drought_bars)
  if (save_fig) save_plot(p, paste0("quantum_yield_", gsub("-","_", filter_soiltype)),
                          soiltype = filter_soiltype, device = device)
  p
}

# ------------------------
# TREE-LEVEL: Chlorophyll
# ------------------------
plot_ts_chl <- function(filter_soiltype = c("both", "inoc-beech", "inoc-robinia"),
                        style = "band",
                        drought_bars = FALSE,
                        save_fig = FALSE, device = "pdf") {
  df <- get_data("tree","chlorophyll") %>%
    { if (filter_soiltype != "both") filter(., soiltype == filter_soiltype) else . } %>%
    filter(!is.na(chl)) %>%
    normalize_factors_tree()
  
  df_sum <- summarize_ts(df, chl, keys = c("species","robinia","precipitation","culture"))
  
  p <- build_timeseries(df, chl, facet = "tree", style = style,
                        ylab = "Chlorophyll (a.u.)",
                        title = paste0("Chlorophyll over time",
                                       if (filter_soiltype!="both") paste0(" (soiltype = ", filter_soiltype, ")")))
  
  yl <- get_ylim("chlorophyll")
  if (!is.null(yl)) p <- p + coord_cartesian(ylim = yl)
  
  p <- add_drought_bars_background(p, df_sum, drought_bars)
  if (save_fig) save_plot(p, paste0("chlorophyll_", gsub("-","_", filter_soiltype)),
                          soiltype = filter_soiltype, device = device)
  p
}

# ------------------------
# TREE-LEVEL: Senescence (%)
# ------------------------
plot_ts_senescence_percent <- function(filter_soiltype = c("both", "inoc-beech", "inoc-robinia"),
                                       style = "band",
                                       drought_bars = FALSE,
                                       save_fig = FALSE, device = "pdf") {
  
  df <- get_data("tree","senescence") %>%
    { if (filter_soiltype != "both") filter(., soiltype == filter_soiltype) else . } %>%
    filter(!is.na(percent_senesced)) %>%
    normalize_factors_tree()
  
  df_sum <- summarize_ts(df, percent_senesced, keys = c("species","robinia","precipitation","culture"))
  
  p <- build_timeseries(df, percent_senesced, facet = "tree", style = style,
                        ylab = "Senescence (%)",
                        title = paste0("Senescence over time",
                                       if (filter_soiltype!="both") paste0(" (soiltype = ", filter_soiltype, ")")))
  
  yl <- get_ylim("senescence")
  if (!is.null(yl)) p <- p + coord_cartesian(ylim = yl)
  
  p <- add_drought_bars_background(p, df_sum, drought_bars)
  if (save_fig) save_plot(p, paste0("senescence_percent_", gsub("-","_", filter_soiltype)),
                          soiltype = filter_soiltype, device = device)
  p
}

# ------------------------
# TREE-LEVEL: Tree condition (-)
# ------------------------
plot_ts_condition <- function(filter_soiltype = c("both", "inoc-beech", "inoc-robinia"),
                              style = "band",
                              drought_bars = FALSE,
                              save_fig = FALSE, device = "pdf") {
  
  df <- get_data("tree","condition") %>%
    { if (filter_soiltype != "both") filter(., soiltype == filter_soiltype) else . } %>%
    filter(!is.na(condition)) %>%
    normalize_factors_tree()
  
  df_sum <- summarize_ts(df, condition, keys = c("species","robinia","precipitation","culture"))
  
  p <- build_timeseries(df, condition, facet = "tree", style = style,
                        ylab = "Tree stress (−)",
                        title = paste0("Tree condition over time",
                                       if (filter_soiltype!="both") paste0(" (soiltype = ", filter_soiltype, ")")))
  
  yl <- get_ylim("condition")
  if (!is.null(yl)) p <- p + coord_cartesian(ylim = yl)
  
  p <- add_drought_bars_background(p, df_sum, drought_bars)
  if (save_fig) save_plot(p, paste0("tree_condition_", gsub("-","_", filter_soiltype)),
                          soiltype = filter_soiltype, device = device)
  p
}

# ------------------------
# SOIL-LEVEL: Soil water content
# ------------------------
plot_ts_swc <- function(filter_soiltype = c("both","inoc-beech","inoc-robinia"),
                        style = "band",
                        drought_bars = FALSE,
                        include_soil_treatment = NULL,
                        save_fig = FALSE, device = "pdf") {
  filter_soiltype <- match.arg(filter_soiltype)
  include_soil_treatment <- alinv_resolve_include_soil_treatment(
    include_soil_treatment = include_soil_treatment,
    soil_filter = filter_soiltype
  )
  show_soiltype <- alinv_should_show_soil_panels(
    soil_filter = filter_soiltype,
    include_soil_treatment = include_soil_treatment
  )
  
  # load + optional soil filter
  df <- get_data("box","soilwater")
  df <- alinv_apply_soil_context(
    df,
    soil_filter = filter_soiltype,
    include_soil_treatment = include_soil_treatment
  )
  df <- df %>%
    dplyr::filter(!is.na(swc)) %>%
    normalize_factors_soil()
  
  # summarize for SE/drought bars
  swc_keys <- if (isTRUE(show_soiltype)) {
    c("soiltype","robinia","precipitation","culture")
  } else {
    c("robinia","precipitation","culture")
  }
  df_sum <- summarize_ts(df, swc, keys = swc_keys)
  
  # build plot
  p <- build_timeseries(df, swc, facet = "soil", style = style,
                        show_soiltype = show_soiltype,
                        ylab = "Soil water content (%)",
                        title = paste0("Soil water content over time",
                                       if (filter_soiltype != "both") paste0(" (soiltype = ", filter_soiltype, ")")))
  yl <- get_ylim("swc")
  if (!is.null(yl)) p <- p + coord_cartesian(ylim = yl)
  
  # optional drought bars & save
  p <- add_drought_bars_background(p, df_sum, show = drought_bars)
  if (save_fig) save_plot(p, "soilwater", soiltype = filter_soiltype, device = device)
  p
}

# ------------------------
# SOIL-LEVEL: Soil respiration
# ------------------------
plot_ts_resp <- function(filter_soiltype = c("both","inoc-beech","inoc-robinia"),
                         style = "band",
                         drought_bars = FALSE,
                         include_soil_treatment = NULL,
                         save_fig = FALSE, device = "pdf") {
  
  filter_soiltype <- match.arg(filter_soiltype)
  include_soil_treatment <- alinv_resolve_include_soil_treatment(
    include_soil_treatment = include_soil_treatment,
    soil_filter = filter_soiltype
  )
  show_soiltype <- alinv_should_show_soil_panels(
    soil_filter = filter_soiltype,
    include_soil_treatment = include_soil_treatment
  )

  # load + optional soil filter
  df <- get_data("box","respiration")
  df <- alinv_apply_soil_context(
    df,
    soil_filter = filter_soiltype,
    include_soil_treatment = include_soil_treatment
  )
  df <- df %>%
    dplyr::filter(!is.na(co2)) %>%
    normalize_factors_soil()

  if (!nrow(df)) {
    message("No soil respiration rows available after filtering for soil type: ", filter_soiltype)
    return(
      alinv_empty_plot(
        title = "No soil respiration data",
        subtitle = paste0("No rows available after filtering for soil type: ", filter_soiltype)
      )
    )
  }
  
  # summarize for SE/drought bars
  resp_keys <- if (isTRUE(show_soiltype)) {
    c("soiltype","robinia","precipitation","culture")
  } else {
    c("robinia","precipitation","culture")
  }
  df_sum <- summarize_ts(df, co2, keys = resp_keys)
  
  # build plot
  p <- build_timeseries(df, co2, facet = "soil", style = style,
                        show_soiltype = show_soiltype,
                        ylab = expression(CO[2]~flux~(µmol~m^-2~s^-1)),
                        title = paste0("Soil respiration over time",
                                       if (filter_soiltype != "both") paste0(" (soiltype = ", filter_soiltype, ")")))
  yl <- get_ylim("respiration")
  if (!is.null(yl)) p <- p + coord_cartesian(ylim = yl)
  
  # optional drought bars & save
  p <- add_drought_bars_background(p, df_sum, show = drought_bars)
  if (save_fig) save_plot(p, "soil_respiration", soiltype = filter_soiltype, device = device)
  p
}

# ------------------------
# ========================================
# SOIL-LEVEL: Soil water potential
# ========================================

# Wrangle soil water potential data from raw sensor file
# Loads 10-minute soil data, aggregates to daily means, joins with box metadata
wrangle_soil_data <- function() {
  soil_file <- .resolve_path("data/raw/sensor_data/soil_10min.dat")
  if (!file.exists(soil_file)) {
    stop("Soil water potential raw file not found: ", soil_file, call. = FALSE)
  }

  # Load and aggregate to daily mean
  soil10 <- utils::read.table(
    soil_file,
    sep = ",",
    header = TRUE,
    skip = 2,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::slice(-1)

  soil10_daily <- soil10 %>%
    dplyr::mutate(
      datetime = lubridate::ymd_hms(.data$TS),
      date = as.Date(.data$datetime)
    ) %>%
    dplyr::select(.data$date, dplyr::starts_with("kPa")) %>%
    dplyr::mutate(
      dplyr::across(
        -dplyr::all_of("date"),
        ~ as.numeric(dplyr::na_if(.x, "NAN"))
      )
    ) %>%
    dplyr::group_by(.data$date) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), mean, na.rm = TRUE), .groups = "drop")

  # Map sensors to box labels
  sensor_lookup <- tibble::tibble(
    sensor_nr = 1:24,
    boxlabel = c(
      "b2-p6-c1", "b2-p6-c2", "b2-p6-c3", "b2-p6-c4",
      "b2-p5-c1", "b2-p5-c2", "b2-p5-c3", "b2-p5-c4",
      "b2-p4-c1", "b2-p4-c2", "b2-p4-c3", "b2-p4-c4",
      "b2-p3-c1", "b2-p3-c2", "b2-p3-c3", "b2-p3-c4",
      "b2-p2-c1", "b2-p2-c2", "b2-p2-c3", "b2-p2-c4",
      "b2-p1-c1", "b2-p1-c2", "b2-p1-c3", "b2-p1-c4"
    )
  )

  # Pivot to long format and attach sensor lookup
  soil_mp_long <- soil10_daily %>%
    tidyr::pivot_longer(
      cols = dplyr::starts_with("kPa"),
      names_to = "sensor_raw",
      values_to = "soil_mp"
    ) %>%
    dplyr::mutate(
      sensor_nr = dplyr::if_else(
        .data$sensor_raw == "kPa",
        1L,
        as.integer(stringr::str_extract(.data$sensor_raw, "\\d+")) + 1L
      )
    ) %>%
    dplyr::select(.data$date, .data$sensor_nr, .data$soil_mp) %>%
    dplyr::left_join(sensor_lookup, by = "sensor_nr")

  # Get box metadata and join
  box_meta <- get_data("box", "soilwater") %>%
    dplyr::select(
      .data$boxlabel,
      .data$soiltype,
      .data$robinia,
      .data$culture,
      .data$precipitation,
      .data$species_mix
    ) %>%
    dplyr::distinct()

  soil_daily_long <- soil_mp_long %>%
    dplyr::left_join(box_meta, by = "boxlabel")

  return(soil_daily_long)
}

# Plot soil water potential from wrangled data
# Takes daily aggregated soil data with metadata and creates plot
# show: "both" = raw + smoothed, "raw" = raw trajectories only, "smooth" = smoothed trends only
plot_soil_mp <- function(df_soil, 
                         filter_soiltype = c("both","inoc-beech","inoc-robinia"),
                         show = c("both", "raw", "smooth"),
                         include_soil_treatment = NULL,
                         value_scale = c("log", "raw")) {
  filter_soiltype <- match.arg(filter_soiltype)
  show <- match.arg(show)
  value_scale <- match.arg(value_scale)
  include_soil_treatment <- alinv_resolve_include_soil_treatment(
    include_soil_treatment = include_soil_treatment,
    soil_filter = filter_soiltype
  )
  show_soiltype <- alinv_should_show_soil_panels(
    soil_filter = filter_soiltype,
    include_soil_treatment = include_soil_treatment
  )

  df_plot <- df_soil %>%
    dplyr::filter(!is.na(.data$soil_mp)) %>%
    {
      if (identical(value_scale, "log")) {
        dplyr::filter(., .data$soil_mp < 0)
      } else {
        .
      }
    } %>%
    dplyr::mutate(
      soil_mp_value = if (identical(value_scale, "log")) -log10(-.data$soil_mp) else .data$soil_mp,
      soiltype = gsub("-", "_", .data$soiltype),
      robinia = factor(.data$robinia, levels = c("without-robinia", "with-robinia")),
      precipitation = factor(.data$precipitation, levels = c("control", "drought")),
      culture = factor(.data$culture, levels = c("mono", "mixed"))
    )

  df_plot <- df_plot %>%
    dplyr::mutate(soiltype = gsub("-", "_", .data$soiltype)) %>%
    {
      if (filter_soiltype != "both") {
        dplyr::filter(., .data$soiltype == gsub("-", "_", filter_soiltype))
      } else {
        .
      }
    }

  if (!nrow(df_plot)) {
    message("No soil water potential rows available for soil type: ", filter_soiltype)
    return(
      alinv_empty_plot(
        title = "No soil water potential data",
        subtitle = paste0("No rows available after filtering for soil type: ", filter_soiltype)
      )
    )
  }

  y_min <- min(df_plot$soil_mp_value, na.rm = TRUE)
  y_max <- max(df_plot$soil_mp_value, na.rm = TRUE)
  y_bar <- y_min - 0.04 * (y_max - y_min)

  p <- ggplot2::ggplot(
    df_plot,
    ggplot2::aes(
      x = .data$date,
      y = .data$soil_mp_value,
      color = .data$precipitation,
      linetype = .data$culture,
      group = interaction(.data$precipitation, .data$culture)
    )
  )
  
  # Aggregate to one line per treatment group (mean across boxes per date)
  group_cols <- c("date", "precipitation", "culture", "robinia")
  if (isTRUE(show_soiltype)) {
    group_cols <- c(group_cols, "soiltype")
  }

  df_agg <- df_plot %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::summarise(soil_mp_value = mean(.data$soil_mp_value, na.rm = TRUE), .groups = "drop")

  # Add layers based on show parameter
  if (show %in% c("both", "raw")) {
    if (show == "raw") {
      # Single aggregated line per treatment
      p <- p + ggplot2::geom_line(data = df_agg, linewidth = 0.5)
    } else {
      # Thin faint per-box lines as background
      p <- p + ggplot2::geom_line(
        ggplot2::aes(group = interaction(.data$precipitation, .data$culture, .data$boxlabel)),
        alpha = 0.16,
        linewidth = 0.35
      )
    }
  }
  
  if (show %in% c("both", "smooth")) {
    p <- p +
      ggplot2::geom_line(
        data = df_agg,
        stat = "smooth",
        method = "loess",
        se = FALSE,
        linewidth = 1.1,
        span = 0.35
      )
  }
  
  visible_periods <- purrr::keep(DROUGHT_PERIODS, function(win) {
    win_start <- as.Date(win[[1]])
    win_end <- as.Date(win[[2]])
    win_start <= max(df_plot$date, na.rm = TRUE) && win_end >= min(df_plot$date, na.rm = TRUE)
  })

  p <- p +
    {
      if (isTRUE(show_soiltype)) {
        ggplot2::facet_grid(
          rows = ggplot2::vars(.data$soiltype),
          cols = ggplot2::vars(.data$robinia),
          labeller = ggplot2::labeller(soiltype = SOILTYPE_LABELS)
        )
      } else {
        ggplot2::facet_grid(. ~ .data$robinia)
      }
    } +
    scale_precip_color() +
    scale_culture_linetype() +
    ggplot2::scale_y_continuous(
      name = if (identical(value_scale, "log")) "log10(phi_soil) [kPa]" else "Soil water potential [kPa]",
      breaks = pretty(c(y_min, y_max), n = 6)
    ) +
    ggplot2::scale_x_date(date_labels = "%m/%y") +
    ggplot2::labs(
      x = "Date",
      title = paste0(
        if (identical(value_scale, "log")) "Soil water potential: " else "Soil water potential (raw kPa): ",
        switch(show,
               "both" = "raw trajectories + smoothed treatment means",
               "raw"  = "mean sensor trajectories per treatment",
               "smooth" = "smoothed treatment trends"),
        if (filter_soiltype != "both") paste0(" (soiltype = ", filter_soiltype, ")")
      )
    ) +
    theme_common() +
    ggplot2::coord_cartesian(ylim = c(y_bar, y_max + 0.03 * (y_max - y_min)))

  for (win in visible_periods) {
    win_start <- max(as.Date(win[[1]]), min(df_plot$date, na.rm = TRUE))
    win_end <- min(as.Date(win[[2]]), max(df_plot$date, na.rm = TRUE))
    if (win_start > win_end) next

    p <- p + ggplot2::annotate(
      "segment",
      x = win_start,
      xend = win_end,
      y = y_bar,
      yend = y_bar,
      size = 2.2,
      color = "orange",
      lineend = "round",
      alpha = 0.9
    )
  }

  return(p)
}

# Wrapper: wrangle and plot soil water potential data
plot_ts_soil_mp <- function(filter_soiltype = c("both","inoc-beech","inoc-robinia"),
                            show = c("both", "raw", "smooth"),
                            include_soil_treatment = NULL,
                            value_scale = c("log", "raw"),
                            save_fig = FALSE,
                            device = "pdf") {
  filter_soiltype <- match.arg(filter_soiltype)
  show <- match.arg(show)
  value_scale <- match.arg(value_scale)

  df_soil <- wrangle_soil_data()
  p <- plot_soil_mp(
    df_soil,
    filter_soiltype = filter_soiltype,
    show = show,
    include_soil_treatment = include_soil_treatment,
    value_scale = value_scale
  )

  if (!is.null(p) && save_fig) {
    save_plot(
      p,
      stub = paste0("soil_water_potential_", value_scale, "_", show, "_", gsub("-", "_", filter_soiltype)),
      soiltype = filter_soiltype,
      subdir = "timeseries",
      device = device
    )
  }

  p
}

# ------------------------
# TREE-LEVEL: Phenology (DOY across stages)
# ------------------------
plot_ts_phenology <- function(style = c("errorbar", "band"),
                              filter_soiltype = c("both","inoc-beech","inoc-robinia"),
                              save_fig = FALSE, device = "pdf") {
  style <- match.arg(style)
  filter_soiltype <- match.arg(filter_soiltype)
  
  df <- get_data("tree", "phenology_transitions") %>%
    alinv_filter_by_soil(soil_filter = filter_soiltype) %>%
    dplyr::filter(!is.na(.data$doy), !is.na(.data$stage_date), !is.na(.data$stage)) %>%
    dplyr::mutate(
      stage_num = as.integer(.data$stage),
      stage = factor(.data$stage_num,
                     levels = 1:4,
                     labels = paste("Stage", 1:4)),
      date = as.Date(.data$stage_date)
    ) %>%
    normalize_factors_tree()
  
  # Mean ± SE per stage × treatment
  df_sum <- df %>%
    dplyr::group_by(species, robinia, precipitation, culture, stage) %>%
    dplyr::summarise(
      mean = mean(doy, na.rm = TRUE),
      sd   = sd(doy,   na.rm = TRUE),
      n    = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      se = sd / sqrt(pmax(n, 1)),
      stage_num = as.numeric(stage),
      stage_lab = levels(stage)[stage_num]
    )

  err_alpha <- if (style == "band") 0.45 else 0.75
  err_size <- if (style == "band") 0.85 else 1.0
  cap_half_height <- if (style == "band") 0.055 else 0.075

  p <- ggplot(
    df_sum,
    aes(
      x = mean,
      y = stage_num,
      color = precipitation,
      linetype = culture,
      group = interaction(precipitation, culture)
    )
  ) +
    geom_segment(
      aes(
        x = mean - se,
        xend = mean + se,
        y = stage_num,
        yend = stage_num
      ),
      alpha = err_alpha,
      linewidth = err_size,
      lineend = "round"
    ) +
    geom_segment(
      aes(
        x = mean - se,
        xend = mean - se,
        y = stage_num - cap_half_height,
        yend = stage_num + cap_half_height
      ),
      alpha = err_alpha,
      linewidth = err_size,
      lineend = "round"
    ) +
    geom_segment(
      aes(
        x = mean + se,
        xend = mean + se,
        y = stage_num - cap_half_height,
        yend = stage_num + cap_half_height
      ),
      alpha = err_alpha,
      linewidth = err_size,
      lineend = "round"
    ) +
    geom_line(size = PLOT_STYLE$sizes$line) +
    geom_point(size = PLOT_STYLE$sizes$point) +
    facet_tree() +
    scale_precip_color() +
    scale_culture_linetype() +
    scale_y_continuous(
      breaks = 1:4,
      labels = paste("Stage", 1:4),
      limits = c(0.8, 4.2)
    ) +
    labs(
      x = "Mean DOY (± SE)",
      y = "Phenological stage",
      title = paste0("Phenology progression (mean ± SE)",
                     if (filter_soiltype!="both") paste0(" — soiltype: ", filter_soiltype))
    ) +
    theme_common()

  x_lim <- get_ylim("phenology_doy")
  data_x <- range(df_sum$mean + c(-1, 1) * df_sum$se, na.rm = TRUE)
  data_x <- c(data_x[1] - 2, data_x[2] + 3)
  if (!is.null(x_lim)) {
    x_final <- c(min(x_lim[1], data_x[1]), max(x_lim[2], data_x[2]))
  } else {
    x_final <- data_x
  }
  p <- p + coord_cartesian(xlim = x_final, ylim = c(0.8, 4.2))
  
  if (save_fig) save_plot(
    p,
    paste0("phenology_progression_", style, "_", gsub("-","_", filter_soiltype)),
    soiltype = filter_soiltype,
    device = device
  )
  
  p
}

# -------- RAW + SMOOTH builder (raw trajectories + smoothed treatment curves) --------
# Generic function: plot raw individual trajectories and smoothed treatment curves
# Useful for showing individual variation and treatment trends simultaneously
build_raw_smooth <- function(df_in, y, facet = c("tree","soil"),
                             ylab = NULL, title = NULL,
                             x_date_fmt = "%m/%y",
                             smooth_method = "loess",
                             smooth_span = 0.35,
                             raw_alpha = 0.16,
                             raw_linewidth = 0.35,
                             smooth_linewidth = 1.1) {
  facet <- match.arg(facet)
  
  # Grouping for facets
  if (facet == "tree") {
    facet_vars <- ggplot2::vars(.data$species, .data$robinia)
    labeller_spec <- ggplot2::labeller(
      robinia = c(`without-robinia`="without robinia", `with-robinia`="with robinia")
    )
    facet_grid <- ggplot2::facet_grid(rows = ggplot2::vars(.data$species), cols = ggplot2::vars(.data$robinia), labeller = labeller_spec)
  } else {
    facet_grid <- ggplot2::facet_grid(rows = ggplot2::vars(.data$soiltype), cols = ggplot2::vars(.data$robinia))
  }
  
  # Plot: raw trajectories (thin, semi-transparent) + smoothed treatment curves
  p <- ggplot(
    df_in,
    ggplot2::aes(
      x = .data$date,
      y = {{y}},
      color = .data$precipitation,
      linetype = .data$culture,
      group = interaction(.data$precipitation, .data$culture, if (facet == "tree") .data$tree_id else .data$boxlabel)
    )
  ) +
    ggplot2::geom_line(alpha = raw_alpha, linewidth = raw_linewidth) +
    ggplot2::geom_smooth(
      ggplot2::aes(group = interaction(.data$precipitation, .data$culture)),
      method = smooth_method,
      se = FALSE,
      linewidth = smooth_linewidth,
      span = smooth_span,
      color = NA  # Remove color override from geom_smooth to inherit from aes
    ) +
    facet_grid +
    scale_precip_color() +
    scale_culture_linetype() +
    scale_precip_fill() +
    ggplot2::scale_x_date(date_labels = x_date_fmt) +
    ggplot2::labs(x = "Date", y = ylab %||% as_name(ensym(y)), title = title) +
    theme_common()
  
  return(p)
}

plot_ts_senescence_chl <- function(filter_soiltype = c("both", "inoc-beech", "inoc-robinia"),
                                   style = c("band","errorbar"),
                                   drought_bars = FALSE,
                                   save_fig = FALSE,
                                   device = "pdf") {
  style <- match.arg(style)
  filter_soiltype <- match.arg(filter_soiltype)
  
  # load + optional soil filter
  df <- get_data("tree", "senescence")
  if (filter_soiltype != "both") {
    df <- dplyr::filter(df, soiltype == filter_soiltype)
  }
  df <- dplyr::filter(df, !is.na(chlavg))
  df <- normalize_factors_tree(df)
  
  # mean ± SE by treatment × date (for drought bars & limits)
  df_sum <- summarize_ts(df, chlavg,
                         keys = c("species","robinia","precipitation","culture"))
  
  # build plot via the common builder
  p <- build_timeseries(
    df_in = df, y = chlavg, facet = "tree", style = style,
    ylab = "Chlorophyll (a.u.)",
    title = paste0(
      "Chlorophyll during senescence",
      if (filter_soiltype != "both") paste0(" (soiltype = ", filter_soiltype, ")")
    )
  )
  
  # apply global y-limits if defined
  yl <- get_ylim("chlorophyll")
  if (!is.null(yl)) p <- p + coord_cartesian(ylim = yl)
  
  # optional drought bars (background) using your helper
  drought_bars <-  FALSE
  p <- add_drought_bars_background(p, df_sum, show = drought_bars)
  
  # optional save
  if (save_fig) {
    stub <- paste0("senescence_chlorophyll_", gsub("-", "_", filter_soiltype))
    save_plot(p, stub = stub, soiltype = filter_soiltype,
              subdir = "timeseries", device = device)
  }
  
  p
}

# ============================================================
# BOX PLOT WRAPPERS (SLA, Growth increment@date, Soil CN/iso)
# ============================================================

# Uses:
# - normalize_factors_tree(), normalize_factors_soil()
# - build_boxplot()
# - save_plot()  (we pass subdir="boxplots")
# - PLOT_STYLE (for global orders/colors/linetypes/alpha)

# -------------------------------
# 1) SLA (tree-level; no date)
# -------------------------------
boxplot_sla <- function(save_fig = FALSE, device = "pdf") {
  df <- get_data("tree", "sla") |>
    normalize_factors_tree()
  
  p <- build_boxplot(
    df, sla_mm2_per_mg, facet = "tree",
    ylab  = expression(SLA~(mm^2~mg^-1)),
    title = "Specific Leaf Area by precipitation (color) and culture (alpha)"
  )
  
  if (save_fig) save_plot(
    p, stub = "sla_boxplot", soiltype = "both",
    subdir = "boxplots", device = device
  )
  p
}

# ----------------------------------------------------------------
# 2) Growth increment on a chosen date (tree-level; diameter/height)
#    - metric: "diameter" | "height"
#    - target_date: Date or "YYYY-MM-DD"
#    - filter_soiltype: "both" | "inoc-beech" | "inoc-robinia"
#    - validates that metric exists on requested date; lists available
# ----------------------------------------------------------------
boxplot_growth_increment_on <- function(target_date,
                                        metric = c("diameter","height"),
                                        filter_soiltype = c("both", "inoc-beech", "inoc-robinia"),
                                        save_fig = FALSE, device = "pdf") {
  metric <- match.arg(metric)
  target_date <- as.Date(target_date)
  
  # load + optional soil filter
  df <- get_data("tree", "growth")
  if (filter_soiltype != "both") {
    df <- dplyr::filter(df, soiltype == filter_soiltype)
  }
  df <- dplyr::mutate(df, date = as.Date(date))
  
  # dates available for this metric
  available_dates <- df |>
    dplyr::filter(!is.na(.data[[metric]])) |>
    dplyr::distinct(date) |>
    dplyr::arrange(date) |>
    dplyr::pull(date)
  
  if (!(target_date %in% available_dates)) {
    stop(sprintf(
      "No '%s' measurements on %s.\nAvailable dates for '%s': %s",
      metric, format(target_date, "%Y-%m-%d"),
      metric, paste(format(available_dates, "%Y-%m-%d"), collapse = ", ")
    ))
  }
  
  # increment up to the chosen date (last - first per tree)
  df_inc <- df |>
    dplyr::filter(!is.na(.data[[metric]]), date <= target_date) |>
    dplyr::group_by(tree_id) |>
    dplyr::arrange(date, .by_group = TRUE) |>
    dplyr::summarise(
      inc           = dplyr::last(.data[[metric]]) - dplyr::first(.data[[metric]]),
      species       = dplyr::first(species),
      robinia       = dplyr::first(robinia),
      precipitation = dplyr::first(precipitation),
      culture       = dplyr::first(culture),
      .groups = "drop"
    ) |>
    normalize_factors_tree()
  
  ylab <- if (metric == "diameter") "Diameter increment (mm)" else "Height increment (cm)"
  ttl  <- sprintf("%s increment up to %s (precipitation = color, culture = alpha)",
                  tools::toTitleCase(metric), format(target_date, "%d.%m.%Y"))
  
  p <- build_boxplot(
    df_inc, inc, facet = "tree",
    ylab = ylab,
    title = ttl
  )
  
  yl <- get_ylim(paste0(metric, "_boxplot"))
  if (!is.null(yl)) p <- p + coord_cartesian(ylim = yl)
  
  if (save_fig) {
    stub <- sprintf("growth_%s_increment_on_%s_%s",
                    metric, format(target_date, "%Y%m%d"), gsub("-","_", filter_soiltype))
    save_plot(p, stub = stub, soiltype = filter_soiltype,
              subdir = "boxplots", device = device)
  }
  p
}

# -----------------------------------------------------------------------
# 3) Soil C/N & isotopes (box-level)
#    - facets 4 variables; color/fill by precipitation (consistent legend)
# -----------------------------------------------------------------------
boxplot_soil_cn_isotopes <- function(filter_soiltype = NULL,
                                     include_soil_treatment = NULL,
                                     save_fig = FALSE, device = "pdf") {
  filter_soiltype <- alinv_resolve_soil_filter(filter_soiltype)
  include_soil_treatment <- alinv_resolve_include_soil_treatment(
    include_soil_treatment = include_soil_treatment,
    soil_filter = filter_soiltype
  )

  df <- get_data("box", "cn_isotopes") |>
    alinv_apply_soil_context(
      soil_filter = filter_soiltype,
      include_soil_treatment = include_soil_treatment
    ) |>
    normalize_factors_soil()
  
  # Long form for facetting variables
  df_long <- df |>
    dplyr::select(soiltype, robinia, culture, precipitation,
                  c_perc, n_perc, d13c_permille, d15n_permille) |>
    tidyr::pivot_longer(
      c(c_perc, n_perc, d13c_permille, d15n_permille),
      names_to = "variable", values_to = "value"
    )
  
  var_labels <- c(
    c_perc        = "C (%)",
    n_perc        = "N (%)",
    d13c_permille = "δ13C (‰)",
    d15n_permille = "δ15N (‰)"
  )
  
  # Build baseline (soil facets) then facet variables on top
  p <- build_boxplot(
    df_long, value, facet = "soil",
    ylab = NULL,
    title = "Soil C/N and isotopes by treatment",
    show_soiltype = alinv_should_show_soil_panels(
      soil_filter = filter_soiltype,
      include_soil_treatment = include_soil_treatment
    )
  ) +
    facet_wrap(~ variable, scales = "free_y",
               labeller = labeller(variable = var_labels), nrow = 2)
  
  if (save_fig) {
    save_plot(p, stub = "soil_CN_isotopes",
              soiltype = "both", subdir = "boxplots", device = device)
  }
  p
}

boxplot_soil_cn_isotopes_simple <- function(save_fig = FALSE,
                                            device = "pdf",
                                            palette = c(inoc_beech = "purple",
                                                        inoc_robinia = "orange")) {
  # Load
  df <- get_data("box", "cn_isotopes")
  
  # Long format
  df_long <- df %>%
    dplyr::mutate(
      soiltype = gsub("-", "_", soiltype),
      soiltype = factor(soiltype, levels = c("inoc_beech", "inoc_robinia"))
    ) %>%
    dplyr::select(soiltype, c_perc, n_perc, d13c_permille, d15n_permille) %>%
    tidyr::pivot_longer(
      cols = c(c_perc, n_perc, d13c_permille, d15n_permille),
      names_to = "variable", values_to = "value"
    )
  
  var_labels <- c(
    c_perc        = "C (%)",
    n_perc        = "N (%)",
    d13c_permille = "d13C (\u2030)",
    d15n_permille = "d15N (\u2030)"
  )
  
  # Base plot
  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = soiltype, y = value, fill = soiltype, color = soiltype)) +
    ggplot2::geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.15, alpha = 0.6, size = 1.3) +
    ggplot2::facet_wrap(~ variable,
                        scales = "free_y",
                        labeller = ggplot2::labeller(variable = var_labels),
                        nrow = 2) +
    ggplot2::scale_fill_manual(values = palette, name = "Soil type") +
    ggplot2::scale_color_manual(values = palette, guide = "none") +
    ggplot2::labs(x = "Soil type", y = NULL,
                  title = "Isotopes and elemental % by soil type") +
    theme_common()
  
  # Optional t-test p-values if ggpubr is available
  if (requireNamespace("ggpubr", quietly = TRUE)) {
    p <- p + ggpubr::stat_compare_means(
      method = "t.test", label = "p.format",
      comparisons = list(c("inoc_beech", "inoc_robinia"))
    )
  }
  
  # Save if requested
  if (isTRUE(save_fig)) {
    save_plot(p, stub = "soil_CN_isotopes_simple",
              soiltype = "both", subdir = "boxplots", device = device)
  }
  
  p
}
