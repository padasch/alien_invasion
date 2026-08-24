# Shared data, plotting, and bootstrap helpers for Supplementary_v1 figures.

if (!exists("SUPP_SCRIPT_FILE", inherits = FALSE)) {
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(arg)) stop("Run the figure script with Rscript.", call. = FALSE)
  SUPP_SCRIPT_FILE <- normalizePath(
    gsub("~\\+~", " ", sub("^--file=", "", arg[[1]])),
    winslash = "/", mustWork = TRUE
  )
}

SUPP_SCRIPT_DIR <- dirname(SUPP_SCRIPT_FILE)
PROJECT_ROOT <- normalizePath(file.path(SUPP_SCRIPT_DIR, "..", ".."), winslash = "/", mustWork = TRUE)
SUPP_V1_OUTPUT <- file.path(PROJECT_ROOT, "output", "supplementary")
SUPP_MODEL_CACHE <- Sys.getenv(
  "ALINV_SUPP_MODEL_CACHE",
  unset = file.path(tempdir(), "alinv-supplementary-v1-model-cache")
)
dir.create(SUPP_V1_OUTPUT, recursive = TRUE, showWarnings = FALSE)
dir.create(SUPP_MODEL_CACHE, recursive = TRUE, showWarnings = FALSE)
setwd(PROJECT_ROOT)

renv_lib <- Sys.glob(file.path(PROJECT_ROOT, "renv", "library", "*", "R-*", "*"))
if (length(renv_lib)) .libPaths(c(normalizePath(renv_lib[[1]], winslash = "/", mustWork = TRUE), .libPaths()))

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(readr)
  library(ggplot2); library(lme4); library(patchwork); library(lubridate)
})
suppressMessages(suppressPackageStartupMessages({
  auxiliary_functions <- file.path(PROJECT_ROOT, "scripts", "auxiliary", "functions")
  source(file.path(auxiliary_functions, "_source.R"))
  source(file.path(auxiliary_functions, "2-biomass.R"))
  source(file.path(auxiliary_functions, "3-effect-size-factorial.R"))
  source(file.path(auxiliary_functions, "10-temporal-alignment.R"))
  source(file.path(auxiliary_functions, "11-bootstrap-inference.R"))
}))

SUPP_BOOT_B <- as.integer(Sys.getenv("ALINV_SUPP_BOOT_B", unset = "1000"))
if (!is.finite(SUPP_BOOT_B) || SUPP_BOOT_B < 1L) stop("ALINV_SUPP_BOOT_B must be positive.", call. = FALSE)
SUPP_BOOT_CORES <- as.integer(Sys.getenv("ALINV_SUPP_BOOT_CORES", unset = "4"))
if (!is.finite(SUPP_BOOT_CORES) || SUPP_BOOT_CORES < 1L) SUPP_BOOT_CORES <- 1L
supp_detected_cores <- suppressWarnings(parallel::detectCores(logical = TRUE))
if (!length(supp_detected_cores) || !is.finite(supp_detected_cores) || supp_detected_cores < 1L) supp_detected_cores <- 1L
SUPP_BOOT_CORES <- min(SUPP_BOOT_CORES, supp_detected_cores)

SUPP_SPECIES_LABELS <- c(fagus = "Fagus", quercus = "Quercus", robinia = "Robinia")
SUPP_PRECIP_COLORS <- c(control = "#4F6674", drought = "#D65F5F")
SUPP_EFFECT_COLORS <- c(culture = "#1B9E77", precipitation = "#D95F02", robinia = "#7570B3")
SUPP_EFFECT_LABELS <- c(
  culture = "Culture (mono -> mixed)",
  precipitation = "Precipitation (control -> reduced)",
  robinia = "Robinia (without -> with)"
)
SUPP_TREATMENTS <- c("culture", "precipitation", "robinia")
SUPP_TREATED <- c(culture = "mixed", precipitation = "drought", robinia = "with-robinia")
SUPP_CONTRASTS <- c(culture = "mixed - mono", precipitation = "drought - control", robinia = "with-robinia - without-robinia")
SUPP_SUMMER_START <- as.Date("2025-06-20")
SUPP_SUMMER_END <- as.Date("2025-09-01")
SUPP_DROUGHT <- tibble::tribble(
  ~start, ~end,
  as.Date("2025-06-20"), as.Date("2025-07-02"),
  as.Date("2025-08-12"), as.Date("2025-08-20")
)
SUPP_SEASONS <- tibble::tribble(
  ~start, ~end, ~fill,
  as.Date("2025-03-01"), SUPP_SUMMER_START, "#C9F7C5",
  SUPP_SUMMER_START, SUPP_SUMMER_END, "#A9C9AA",
  SUPP_SUMMER_END, as.Date("2025-12-31"), "#F8F6D6"
)

supp_theme <- function(base_size = 7) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(family = "Helvetica", color = "black"),
      axis.text = element_text(color = "black"),
      axis.line = element_line(linewidth = 0.25),
      axis.ticks = element_line(linewidth = 0.25),
      panel.border = element_rect(fill = NA, color = "grey55", linewidth = 0.25),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.18),
      panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey12", color = "grey12"),
      strip.text = element_text(color = "white", face = "bold"),
      legend.position = "bottom",
      plot.title.position = "plot",
      plot.title = element_text(face = "bold", size = base_size + 1),
      plot.subtitle = element_text(size = base_size)
    )
}

supp_output_dir <- function(figure_id) SUPP_V1_OUTPUT

supp_save <- function(plot, figure_id, width_mm = 160, height_mm = 120, width_limit_mm = 160) {
  if (width_mm > width_limit_mm) stop("Figure width exceeds its configured limit.", call. = FALSE)
  out <- supp_output_dir(figure_id)
  pdf <- file.path(out, paste0(figure_id, ".pdf"))
  ggsave(pdf, plot, width = width_mm, height = height_mm, units = "mm", device = cairo_pdf, bg = "white", limitsize = FALSE)
  invisible(pdf)
}

supp_block <- function(boxlabel) {
  z <- sub("^b([0-9]+)-.*$", "\\1", as.character(boxlabel))
  if (any(z == as.character(boxlabel))) stop("Could not derive block from boxlabel.", call. = FALSE)
  factor(z)
}

supp_boot_p <- alinv_bootstrap_p

# Descriptive treatment-group means and design-aligned container bootstrap CIs.
supp_group_bootstrap <- function(df, value_col, seed, B = SUPP_BOOT_B) {
  keys <- c("species", "robinia", "precipitation", "culture", "date")
  box_means <- df %>%
    filter(is.finite(.data[[value_col]])) %>%
    mutate(block = supp_block(.data$boxlabel)) %>%
    group_by(across(all_of(c(keys, "block", "boxlabel")))) %>%
    summarise(value = mean(.data[[value_col]]), n_trees = n_distinct(.data$tree_id), .groups = "drop")
  groups <- group_split(group_by(box_means, across(all_of(keys))))
  purrr::imap_dfr(groups, function(g, i) {
    by_block <- split(g$value, g$block)
    set.seed(seed + i)
    draws <- replicate(B, mean(unlist(lapply(by_block, function(v) sample(v, length(v), replace = TRUE))), na.rm = TRUE))
    bind_cols(
      g[1, keys, drop = FALSE],
      tibble(
        mean = mean(g$value), lower = alinv_bootstrap_percentile(draws)[[1]],
        upper = alinv_bootstrap_percentile(draws)[[2]], n_containers = nrow(g),
        n_trees = sum(g$n_trees), bootstrap_replicates = B,
        interval = "95% block-stratified container-bootstrap percentile"
      )
    )
  })
}

supp_prepare_series <- function(data_name, value_col, species) {
  df <- get_data("tree", data_name) %>% filter(.data$species %in% species)
  if (!value_col %in% names(df)) stop("Missing response: ", value_col, call. = FALSE)
  df %>%
    filter(!is.na(.data[[value_col]])) %>%
    mutate(
      species = factor(.data$species, levels = species, labels = unname(SUPP_SPECIES_LABELS[species])),
      robinia = factor(.data$robinia, levels = c("without-robinia", "with-robinia"), labels = c("without robinia", "with robinia")),
      precipitation = factor(.data$precipitation, levels = c("control", "drought")),
      culture = factor(.data$culture, levels = c("mono", "mixed")), date = as.Date(.data$date)
    )
}

supp_timeseries_figure <- function(figure_id, data_name, value_col, species, y_label, seed, limits = NULL) {
  df <- supp_prepare_series(data_name, value_col, species)
  src <- supp_group_bootstrap(df, value_col, seed)
  src$species <- factor(src$species, levels = unname(SUPP_SPECIES_LABELS[species]))
  src$robinia <- factor(src$robinia, levels = c("without robinia", "with robinia"))
  season_layers <- pmap(SUPP_SEASONS, function(start, end, fill) {
    annotate("rect", xmin = start, xmax = end, ymin = -Inf, ymax = Inf,
             fill = fill, alpha = .42, color = NA)
  })
  p <- ggplot(src, aes(date, mean, ymin = lower, ymax = upper, color = precipitation,
                        fill = precipitation, linetype = culture,
                        group = interaction(precipitation, culture))) +
    season_layers +
    geom_ribbon(alpha = .14, color = NA) + geom_line(linewidth = .5) +
    facet_grid(species ~ robinia, drop = FALSE) +
    scale_color_manual(values = SUPP_PRECIP_COLORS, labels = c(control = "control", drought = "reduced"), name = "Precipitation") +
    scale_fill_manual(values = SUPP_PRECIP_COLORS, guide = "none") +
    scale_linetype_manual(values = c(mono = "solid", mixed = "42"), labels = c(mono = "mono", mixed = "mixed"), name = "Culture") +
    scale_x_date(date_breaks = "3 months", date_labels = "%m/%y") +
    labs(x = "Date", y = y_label, subtitle = "Means and 95% block-stratified container-bootstrap intervals (1,000 replicates)") +
    supp_theme(7)
  yr <- range(c(src$lower, src$upper), na.rm = TRUE)
  if (!is.null(limits)) { p <- p + coord_cartesian(ylim = limits); ybar <- limits[1] + diff(limits) * .02 } else ybar <- yr[1] - diff(yr) * .03
  for (i in seq_len(nrow(SUPP_DROUGHT))) p <- p + annotate("segment", x = SUPP_DROUGHT$start[i], xend = SUPP_DROUGHT$end[i], y = ybar, yend = ybar, color = "#E69F00", linewidth = 1.2, lineend = "round")
  supp_save(p, figure_id, height_mm = if (length(species) == 3) 145 else if (length(species) == 1) 78 else 112)
}

supp_soil_elemental <- function(figure_id = "fig-s1-1-soil-elemental", seed = 2026082101L) {
  df <- get_data("box", "cn_isotopes") %>%
    mutate(block = supp_block(.data$boxlabel), soiltype = factor(gsub("-", "_", .data$soiltype), levels = c("inoc_beech", "inoc_robinia")), cn_ratio = .data$c_perc / .data$n_perc) %>%
    select(boxlabel, block, soiltype, c_perc, n_perc, cn_ratio)
  long <- pivot_longer(df, c(c_perc, n_perc, cn_ratio), names_to = "metric", values_to = "value")
  tests <- map_dfr(unique(long$metric), function(metric_i) {
    z <- filter(long, .data$metric == metric_i, is.finite(.data$value))
    set.seed(seed + match(metric_i, unique(long$metric)))
    draws <- replicate(SUPP_BOOT_B, {
      zz <- group_split(group_by(z, block)) %>% map_dfr(~.x[sample(seq_len(nrow(.x)), nrow(.x), TRUE), ])
      mean(zz$value[zz$soiltype == "inoc_robinia"], na.rm = TRUE) - mean(zz$value[zz$soiltype == "inoc_beech"], na.rm = TRUE)
    })
    tibble(metric = metric_i, estimate = mean(z$value[z$soiltype == "inoc_robinia"], na.rm = TRUE) - mean(z$value[z$soiltype == "inoc_beech"], na.rm = TRUE),
           lower = alinv_bootstrap_percentile(draws)[[1]], upper = alinv_bootstrap_percentile(draws)[[2]], p_boot = supp_boot_p(draws), bootstrap_replicates = SUPP_BOOT_B)
  })
  labels <- c(c_perc = "C (%)", n_perc = "N (%)", cn_ratio = "C:N")
  ypos <- long %>% group_by(metric) %>% summarise(y = max(value, na.rm = TRUE) + .08 * diff(range(value, na.rm = TRUE)), .groups = "drop") %>% left_join(tests, by = "metric")
  p <- ggplot(filter(long, is.finite(.data$value)), aes(soiltype, value, fill = soiltype, color = soiltype)) +
    geom_boxplot(alpha = .65, outlier.shape = NA, width = .62) + geom_jitter(width = .12, alpha = .65, size = .9) +
    geom_text(data = ypos, aes(x = 1.5, y = y, label = sprintf("bootstrap P = %.3f", p_boot)), inherit.aes = FALSE, size = 2.3) +
    facet_wrap(~metric, scales = "free_y", nrow = 1, labeller = labeller(metric = labels)) +
    scale_fill_manual(values = c(inoc_beech = "#B35CFF", inoc_robinia = "#8DD3F5"), labels = c(inoc_beech = "beech-inoculated", inoc_robinia = "Robinia-inoculated"), name = "Soil type") +
    scale_color_manual(values = c(inoc_beech = "#B35CFF", inoc_robinia = "#8DD3F5"), guide = "none") +
    scale_x_discrete(labels = c(inoc_beech = "Beech", inoc_robinia = "Robinia")) +
    labs(x = "Soil inoculum", y = NULL, subtitle = "Block-stratified container bootstrap (1,000 replicates)") + supp_theme(7)
  supp_save(p, figure_id, height_mm = 82)
}

supp_measurement_schedule <- function(figure_id = "fig-s1-2-measurement-schedule") {
  src <- build_measurement_table(get_measurement_schedule_grid(), file.path(PROJECT_ROOT, "data", "interim"))
  p <- plot_measurement_schedule(
    src,
    show_heatwave = TRUE,
    heatwave_continuous = TRUE,
    merge_senescence_chlorophyll = TRUE
  )
  supp_save(p, figure_id, width_mm = 180, height_mm = 80, width_limit_mm = 180)
}

supp_phenology_progression <- function(figure_id = "fig-s2-6-spring-phenology") {
  species <- c("fagus", "quercus", "robinia")
  d <- get_data("tree", "phenology_transitions") %>% filter(.data$species %in% species, .data$stage %in% 1:4, !is.na(.data$doy)) %>%
    mutate(species = factor(.data$species, levels = species, labels = unname(SUPP_SPECIES_LABELS[species])),
           robinia = factor(.data$robinia, levels = c("without-robinia", "with-robinia"), labels = c("without robinia", "with robinia")),
           precipitation = factor(.data$precipitation, levels = c("control", "drought")), culture = factor(.data$culture, levels = c("mono", "mixed")))
  src <- d %>% group_by(species, robinia, precipitation, culture, stage) %>% summarise(mean_doy = mean(doy), se = sd(doy) / sqrt(n()), n = n(), .groups = "drop")
  p <- ggplot(src, aes(mean_doy, stage, color = precipitation, linetype = culture, group = interaction(precipitation, culture))) +
    geom_segment(aes(x = mean_doy - se, xend = mean_doy + se, yend = stage), alpha = .5) + geom_line(linewidth = .55) + geom_point(size = 1.2) +
    facet_grid(species ~ robinia, drop = FALSE) + scale_color_manual(values = SUPP_PRECIP_COLORS, labels = c(control = "control", drought = "reduced"), name = "Precipitation") +
    scale_linetype_manual(values = c(mono = "solid", mixed = "42"), name = "Culture") + scale_y_continuous(breaks = 1:4, labels = paste("Stage", 1:4)) +
    labs(x = "Mean day of year (+/- SE)", y = "Spring phenology") + supp_theme(7)
  supp_save(p, figure_id, height_mm = 105)
}

supp_biomass <- function(figure_id, species) {
  fp <- file.path(PROJECT_ROOT, "data", "raw", "data_20260309.xlsx")
  d <- wrangle_tree_biomass(fp)
  metric_labels <- c(
    root_biomass = "Root biomass (g)",
    shoot_biomass = "Shoot biomass (g)",
    root_shoot_biomass = "Root:shoot biomass (-)"
  )
  plot_data <- d %>%
    filter(.data$species == species) %>%
    pivot_longer(all_of(names(metric_labels)), names_to = "metric", values_to = "value") %>%
    filter(is.finite(.data$value), !is.na(.data$precipitation), !is.na(.data$culture), !is.na(.data$robinia)) %>%
    mutate(
      metric = factor(.data$metric, levels = names(metric_labels), labels = unname(metric_labels)),
      culture = factor(.data$culture, levels = c("mono", "mixed")),
      precipitation = factor(.data$precipitation, levels = c("control", "drought")),
      robinia = factor(.data$robinia, levels = c("without-robinia", "with-robinia"), labels = c("without robinia", "with robinia"))
    )
  p <- ggplot(plot_data, aes(culture, value, fill = precipitation, color = precipitation, alpha = culture)) +
    geom_boxplot(outlier.shape = NA, width = .65, position = position_dodge2(width = .75, preserve = "single")) +
    geom_point(position = position_jitterdodge(jitter.width = .12, dodge.width = .75), size = 1.25, stroke = 0) +
    facet_grid(metric ~ robinia, scales = "free_y", switch = "y") +
    scale_fill_manual(values = c(control = "grey50", drought = "indianred"), name = "Precipitation") +
    scale_color_manual(values = c(control = "grey50", drought = "indianred"), guide = "none") +
    scale_alpha_manual(values = c(mono = .55, mixed = .9), guide = "none") +
    labs(title = paste("Tree biomass by treatments -", species), x = "Culture", y = NULL) +
    theme_bw(base_size = 9) +
    theme(
      strip.background.x = element_rect(fill = "black", color = "black"),
      strip.text.x = element_text(color = "white", face = "bold"),
      strip.background.y = element_blank(),
      strip.text.y.left = element_text(color = "black", face = "plain", angle = 90),
      strip.placement = "outside", legend.position = "bottom",
      plot.title.position = "plot", plot.title = element_text(face = "bold"),
      plot.margin = margin(5, 5, 5, 10)
    )
  supp_save(p, figure_id, height_mm = 145)
}

# Species-specific temporal LMM with block-stratified container bootstrap.
supp_temporal_formula <- y ~ date + date:robinia + date:precipitation + date:culture + (1 | boxlabel) + (1 | tree_id)

supp_prepare_temporal <- function(data_name, resp_var, species) {
  prepare_df_generic(type = "tree", data_name = data_name, resp_var = resp_var, species_keep = species,
                     standardize_response = TRUE, add_covars = FALSE, soil_type = "both",
                     include_soil_treatment = FALSE, swc_source = "measured") %>%
    mutate(block = supp_block(.data$boxlabel), date = factor(as.character(.data$date), levels = sort(unique(as.character(.data$date)))),
           culture = factor(.data$culture, levels = c("mono", "mixed")), precipitation = factor(.data$precipitation, levels = c("control", "drought")),
           robinia = factor(.data$robinia, levels = c("without-robinia", "with-robinia")), boxlabel = factor(.data$boxlabel), tree_id = factor(.data$tree_id))
}

supp_fit_temporal <- function(d) suppressMessages(suppressWarnings(lmer(supp_temporal_formula, data = d, REML = TRUE,
  control = lmerControl(check.rankX = "stop.deficient", check.conv.singular = "ignore"))))

supp_temporal_template <- function(d) expand_grid(date = as.Date(levels(d$date)), effect = SUPP_TREATMENTS) %>%
  mutate(contrast = unname(SUPP_CONTRASTS[effect]), coefficient = paste0("date", format(date, "%Y-%m-%d"), ":", effect, unname(SUPP_TREATED[effect])), key = paste(date, effect, sep = "||"))

supp_extract_temporal <- function(fit, template) {
  b <- fixef(fit); if (!all(template$coefficient %in% names(b))) stop("Missing temporal coefficients.", call. = FALSE)
  setNames(unname(b[template$coefficient]), template$key)
}

supp_resample_containers <- function(d, seed, rep_id) {
  set.seed(seed); boxes <- unique(select(d, boxlabel, block)); selected <- unlist(lapply(split(as.character(boxes$boxlabel), boxes$block), function(x) sample(x, length(x), TRUE)))
  rows <- split(seq_len(nrow(d)), as.character(d$boxlabel)); groups <- unname(rows[selected]); idx <- unlist(groups); occ <- rep(seq_along(groups), lengths(groups)); z <- d[idx, , drop = FALSE]
  new_box <- sprintf("boot%04d_%03d", rep_id, occ); z$tree_id <- factor(paste0(as.character(z$tree_id), "__", new_box)); z$boxlabel <- factor(new_box); z
}

supp_run_temporal_model <- function(data_name, resp_var, species, seed) {
  id <- paste(data_name, resp_var, species, sep = "__")
  cache <- file.path(SUPP_MODEL_CACHE, paste0(id, "__B", SUPP_BOOT_B, ".rds"))
  if (file.exists(cache)) return(readRDS(cache))
  d <- supp_prepare_temporal(data_name, resp_var, species); template <- supp_temporal_template(d); fit <- supp_fit_temporal(d); point <- supp_extract_temporal(fit, template)
  set.seed(seed); seeds <- sample.int(.Machine$integer.max, SUPP_BOOT_B * 10L); success <- list(); failures <- character(); singular <- 0L; attempt <- 1L
  while (length(success) < SUPP_BOOT_B) {
    need <- SUPP_BOOT_B - length(success); ids <- seq.int(attempt, length.out = need)
    batch <- parallel::mclapply(ids, function(i) tryCatch({
      f <- supp_fit_temporal(supp_resample_containers(d, seeds[i], i)); conv <- (f@optinfo$conv$opt %||% 0L) == 0L && !length(f@optinfo$conv$lme4$messages %||% character())
      if (!conv) stop("nonconverged"); list(est = supp_extract_temporal(f, template), singular = isSingular(f))
    }, error = function(e) e), mc.cores = SUPP_BOOT_CORES, mc.set.seed = FALSE)
    for (x in batch) if (inherits(x, "error")) failures <- c(failures, conditionMessage(x)) else { success[[length(success) + 1L]] <- x$est; singular <- singular + as.integer(x$singular) }
    attempt <- max(ids) + 1L
  }
  draws <- do.call(rbind, success[seq_len(SUPP_BOOT_B)])
  lower <- apply(draws, 2, function(x) alinv_bootstrap_percentile(x)[[1]])
  upper <- apply(draws, 2, function(x) alinv_bootstrap_percentile(x)[[2]])
  effects <- template %>% transmute(date, effect, contrast, key, estimate = point[key], lower = lower[key], upper = upper[key],
    p_boot = vapply(seq_len(ncol(draws)), function(j) supp_boot_p(draws[, j]), numeric(1))[match(key, names(point))], bootstrap_replicates = SUPP_BOOT_B,
    species = species, data_name = data_name, response_var = resp_var)
  status <- tibble(model_id = id, formula = paste(deparse(supp_temporal_formula), collapse = " "), successful = SUPP_BOOT_B,
    attempts = attempt - 1L, failures = length(failures), singular_successful = singular, original_singular = isSingular(fit),
    n_obs = nrow(d), n_containers = n_distinct(d$boxlabel), n_trees = n_distinct(d$tree_id), n_dates = n_distinct(d$date), n_blocks = n_distinct(d$block), seed = seed)
  out <- list(effects = effects, status = status, draws = draws, model = fit); saveRDS(out, cache, compress = "xz"); out
}

supp_temporal_effect_figure <- function(figure_id, data_name, resp_var, species, response_label, seed) {
  results <- map2(species, seed + seq_along(species), ~supp_run_temporal_model(data_name, resp_var, .x, .y))
  effects <- bind_rows(map(results, "effects")) %>%
    mutate(panel = factor(
      .data$species,
      levels = species,
      labels = paste0(unname(SUPP_SPECIES_LABELS[species]), " (", response_label, ")")
    ))
  status <- bind_rows(map(results, "status"))
  vals <- c(effects$lower, effects$upper, 0)
  lim <- range(vals[is.finite(vals)]); lim <- lim + c(-1, 1) * diff(lim) * .08
  drought_marks <- crossing(panel = levels(effects$panel), SUPP_DROUGHT) %>% mutate(y = lim[1] + diff(lim) * .035)
  fig <- ggplot(effects, aes(date, estimate, ymin = lower, ymax = upper, color = effect, fill = effect, group = effect)) +
    geom_hline(yintercept = 0, linetype = "22", linewidth = .25) +
    geom_vline(xintercept = c(SUPP_SUMMER_START, SUPP_SUMMER_END), linetype = "42", linewidth = .3) +
    geom_ribbon(alpha = .16, color = NA) + geom_line(linewidth = .55) + geom_point(size = .75) +
    geom_segment(data = drought_marks, aes(x = start, xend = end, y = y, yend = y), inherit.aes = FALSE,
                 color = "#E69F00", linewidth = 1.2, lineend = "round") +
    facet_wrap(~panel, ncol = 1) +
    scale_color_manual(values = SUPP_EFFECT_COLORS, labels = SUPP_EFFECT_LABELS, breaks = SUPP_TREATMENTS, name = NULL) +
    scale_fill_manual(values = SUPP_EFFECT_COLORS, labels = SUPP_EFFECT_LABELS, breaks = SUPP_TREATMENTS, name = NULL) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b") + coord_cartesian(ylim = lim) +
    labs(x = NULL, y = "Effect size (SD units)", subtitle = "95% container-bootstrap intervals; 1,000 successful block-stratified replicates") +
    supp_theme(7) + theme(plot.margin = margin(4, 4, 4, 10))
  supp_save(fig, figure_id, height_mm = if (length(species) == 1) 82 else 125)
}
