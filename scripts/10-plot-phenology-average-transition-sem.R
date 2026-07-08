#!/usr/bin/env Rscript

renv_lib <- Sys.glob(file.path("renv", "library", "*", "R-*", "*"))
if (length(renv_lib)) {
  .libPaths(c(normalizePath(renv_lib[[1]], winslash = "/", mustWork = TRUE), .libPaths()))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(lme4)
  library(patchwork)
  library(readr)
  library(tidyr)
})

suppressMessages(suppressPackageStartupMessages({
  source("functions/_source.R")
  source("functions/1-summary-figures.R")
  source("functions/7-phenology-transition-models.R")
  source("functions/4-structural-equation-model.R")
}))

# Standalone phenology SEM verification.
# Response: per-tree average DOY of reaching stages 2-4.
# Mediator: per-tree average SWC matched to those transition dates.
analysis_date <- as.character(Sys.Date())
species_vec <- c("fagus", "quercus")
soil_type <- "both"
include_soil_treatment <- FALSE
stages_keep <- 2:4
resp_var <- "phenology_transition_doy"
ALINV_RESPONSE_LABELS[resp_var] <- "Avg. transition DOY"

match_swc_to_transition_dates <- function(df_transition) {
  swc_df <- get_data("box", "soilwater", swc_source = "measured") %>%
    dplyr::mutate(date = as.Date(.data$date)) %>%
    dplyr::arrange(.data$boxlabel, .data$date)

  closest_swc <- function(box_i, date_i) {
    swc_box <- swc_df %>% dplyr::filter(.data$boxlabel == box_i)
    if (!nrow(swc_box) || is.na(date_i)) {
      return(NA_real_)
    }

    d_minus7 <- date_i - 7
    d_plus7 <- date_i + 7

    idx_past_week <- which(swc_box$date <= date_i & swc_box$date >= d_minus7)
    if (length(idx_past_week)) {
      return(swc_box$swc[idx_past_week[which.max(swc_box$date[idx_past_week])]])
    }

    idx_next_week <- which(swc_box$date > date_i & swc_box$date <= d_plus7)
    if (length(idx_next_week)) {
      return(swc_box$swc[idx_next_week[which.min(swc_box$date[idx_next_week])]])
    }

    idx_prior <- which(swc_box$date <= date_i)
    if (length(idx_prior)) {
      return(swc_box$swc[idx_prior[which.max(swc_box$date[idx_prior])]])
    }

    NA_real_
  }

  df_transition %>%
    dplyr::mutate(
      stage_date = as.Date(.data$stage_date),
      swc_transition = mapply(
        closest_swc,
        as.character(.data$boxlabel),
        .data$stage_date,
        USE.NAMES = FALSE
      )
    )
}

prepare_average_transition_sem_data <- function(species_keep) {
  df_transition <- prepare_phenology_transition_data(
    species_keep = species_keep,
    soil_type = soil_type,
    include_soil_treatment = include_soil_treatment,
    stages_keep = stages_keep
  ) %>%
    match_swc_to_transition_dates()

  df_transition %>%
    dplyr::group_by(
      .data$tree_id,
      .data$boxlabel,
      .data$species,
      .data$precipitation,
      .data$robinia,
      .data$culture,
      .data$soiltype
    ) %>%
    dplyr::summarise(
      y = mean(.data$doy, na.rm = TRUE),
      swc = mean(.data$swc_transition, na.rm = TRUE),
      n_stages = dplyr::n(),
      n_swc = sum(!is.na(.data$swc_transition)),
      first_transition_date = min(.data$stage_date, na.rm = TRUE),
      last_transition_date = max(.data$stage_date, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::filter(!is.na(.data$y), !is.na(.data$swc), .data$n_swc > 0) %>%
    dplyr::mutate(
      y_org = .data$y,
      swc_org = .data$swc,
      y = as.numeric(scale(.data$y_org)),
      swc = as.numeric(scale(.data$swc_org)),
      tree_id = factor(.data$tree_id),
      boxlabel = factor(.data$boxlabel),
      species = factor(.data$species, levels = alinv_factor_levels("species")),
      precipitation = factor(.data$precipitation, levels = alinv_factor_levels("precipitation")),
      robinia = factor(.data$robinia, levels = alinv_factor_levels("robinia")),
      culture = factor(.data$culture, levels = alinv_factor_levels("culture")),
      soiltype = alinv_relevel_soiltype(.data$soiltype)
    ) %>%
    droplevels()
}

fit_average_transition_sem <- function(df_sem) {
  rhs_terms <- phenology_transition_terms(
    df_sem,
    include_soil_treatment = include_soil_treatment,
    soil_filter = soil_type
  )
  rhs <- if (length(rhs_terms)) paste(rhs_terms, collapse = " + ") else "1"

  fml_swc <- stats::as.formula(paste("swc ~", rhs, "+ (1 | boxlabel)"))
  fml_response <- stats::as.formula(paste("y ~ swc +", rhs, "+ (1 | boxlabel)"))

  list(
    mod_swc = lme4::lmer(fml_swc, data = df_sem, REML = TRUE),
    mod_response = lme4::lmer(fml_response, data = df_sem, REML = TRUE),
    fml_swc = fml_swc,
    fml_response = fml_response,
    modeled_factors = rhs_terms
  )
}

coef_info <- function(model, coef_name) {
  beta <- lme4::fixef(model)
  vc <- as.matrix(stats::vcov(model))
  if (!coef_name %in% names(beta)) {
    return(list(est = NA_real_, se = NA_real_, p = NA_real_))
  }

  est <- unname(beta[[coef_name]])
  se <- sqrt(unname(vc[coef_name, coef_name]))
  p <- 2 * stats::pnorm(-abs(est / se))
  list(est = est, se = se, p = p)
}

extract_average_transition_sem_effects <- function(species_keep, df_sem, sem_fit) {
  trt_info <- phenology_transition_treatment_info(
    include_soil_treatment = include_soil_treatment,
    soil_filter = soil_type
  ) %>%
    dplyr::filter(.data$effect %in% sem_fit$modeled_factors)

  b_info <- coef_info(sem_fit$mod_response, "swc")
  b <- b_info$est
  se_b <- b_info$se

  dplyr::bind_rows(lapply(seq_len(nrow(trt_info)), function(i) {
    info_i <- trt_info[i, ]
    coef_name <- paste0(info_i$effect[[1]], info_i$treatment_level[[1]])

    a_info <- coef_info(sem_fit$mod_swc, coef_name)
    c_info <- coef_info(sem_fit$mod_response, coef_name)

    indirect <- a_info$est * b
    var_indirect <- (b^2) * (a_info$se^2) + (a_info$est^2) * (se_b^2)
    se_indirect <- sqrt(var_indirect)
    p_indirect <- 2 * stats::pnorm(-abs(indirect / se_indirect))

    total <- c_info$est + indirect
    se_total <- sqrt(c_info$se^2 + var_indirect)
    p_total <- 2 * stats::pnorm(-abs(total / se_total))

    tibble::tibble(
      species = species_keep,
      effect = info_i$effect[[1]],
      treatment_label = info_i$treatment_label[[1]],
      contrast_label = info_i$contrast_label[[1]],
      a_treatment_to_swc = a_info$est,
      se_a = a_info$se,
      p_a = a_info$p,
      b_swc_to_transition = b,
      se_b = se_b,
      p_b = b_info$p,
      direct = c_info$est,
      se_direct = c_info$se,
      p_direct = c_info$p,
      indirect = indirect,
      se_indirect = se_indirect,
      p_indirect = p_indirect,
      total = total,
      se_total = se_total,
      p_total = p_total,
      n_trees = dplyr::n_distinct(df_sem$tree_id),
      n_boxes = dplyr::n_distinct(df_sem$boxlabel),
      y_sd = stats::sd(df_sem$y_org, na.rm = TRUE),
      swc_sd = stats::sd(df_sem$swc_org, na.rm = TRUE),
      effect_scale = "standardized_y_swc",
      singular_swc = lme4::isSingular(sem_fit$mod_swc, tol = 1e-4),
      singular_response = lme4::isSingular(sem_fit$mod_response, tol = 1e-4),
      plot_order = info_i$plot_order[[1]]
    )
  }))
}

as_existing_sem_effect_table <- function(effects_df) {
  effects_df %>%
    dplyr::transmute(
      factor = .data$effect,
      a = .data$a_treatment_to_swc,
      se_a = .data$se_a,
      p_a = .data$p_a,
      b = .data$b_swc_to_transition,
      se_b = .data$se_b,
      p_b = .data$p_b,
      c_direct = .data$direct,
      se_c = .data$se_direct,
      p_c = .data$p_direct,
      indirect = .data$indirect,
      se_ind = .data$se_indirect,
      p_ind = .data$p_indirect,
      total = .data$total,
      se_tot = .data$se_total,
      p_tot = .data$p_total
    )
}

sem_results <- lapply(species_vec, function(species_i) {
  df_sem <- prepare_average_transition_sem_data(species_i)
  fit <- fit_average_transition_sem(df_sem)
  effects <- extract_average_transition_sem_effects(species_i, df_sem, fit)
  list(data = df_sem, fit = fit, effects = effects)
})
names(sem_results) <- species_vec

sem_effects <- dplyr::bind_rows(lapply(sem_results, `[[`, "effects"))

if (!nrow(sem_effects)) {
  stop("No phenology SEM effects could be estimated.", call. = FALSE)
}

output_dir <- file.path(
  .alinv_project_root(),
  "output",
  analysis_date,
  "phenology-average-transition-sem"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

effects_file <- file.path(output_dir, "phenology-average-transition-sem-effects.csv")
matrix_file <- file.path(output_dir, "phenology-average-transition-sem-matrix-data.csv")

readr::write_csv(sem_effects, effects_file)
unlink(file.path(output_dir, "phenology-average-transition-sem.png"))
unlink(file.path(output_dir, "phenology-average-transition-sem-heatmap-*.png"))

matrix_data <- dplyr::bind_rows(lapply(species_vec, function(species_i) {
  effects_i <- sem_effects %>%
    dplyr::filter(.data$species == species_i) %>%
    as_existing_sem_effect_table()

  build_sem_matrix_data(
    effects_main = effects_i,
    effects_int = NULL,
    resp_var = resp_var,
    species = species_i,
    include_interaction = FALSE,
    swc_source = "measured",
    phase = "average_transition"
  )
})) %>%
  dplyr::mutate(effect_scale = "standardized_y_swc")
readr::write_csv(matrix_data, matrix_file)

graph_plots <- list()
for (species_i in species_vec) {
  effects_i <- sem_effects %>%
    dplyr::filter(.data$species == species_i) %>%
    as_existing_sem_effect_table()

  fit_i <- sem_results[[species_i]]$fit
  graph_bits <- build_sem_graph_components(
    effects_main = effects_i,
    effects_int = NULL,
    sem_mod = NULL,
    mod_swc = fit_i$mod_swc,
    mod_resp = fit_i$mod_response,
    resp_var = resp_var,
    species = species_i,
    soil_type = soil_type,
    include_interaction = FALSE,
    modeled_factors = fit_i$modeled_factors,
    p_sig = 0.1
  )
  graph_bits$metrics$title <- paste0(
    "SEM: average phenology transition timing (species: ",
    species_i,
    ", soil: ",
    soil_type,
    ")"
  )

  readr::write_csv(graph_bits$nodes, file.path(output_dir, paste0("phenology-average-transition-sem-", species_i, "-graph-nodes.csv")))
  readr::write_csv(graph_bits$edges, file.path(output_dir, paste0("phenology-average-transition-sem-", species_i, "-graph-edges.csv")))
  readr::write_csv(graph_bits$metrics, file.path(output_dir, paste0("phenology-average-transition-sem-", species_i, "-graph-metrics.csv")))

  graph_plot <- plot_sem_graph_components(
    nodes_df = graph_bits$nodes,
    edges_df = graph_bits$edges,
    metrics_df = graph_bits$metrics
  )
  graph_plots[[species_i]] <- graph_plot

  ggsave(
    file.path(output_dir, paste0("phenology-average-transition-sem-", species_i, "-graph.png")),
    graph_plot,
    width = 9,
    height = 5.6,
    dpi = 300
  )
}

combined_graph <- patchwork::wrap_plots(graph_plots, ncol = 1) +
  patchwork::plot_annotation(tag_levels = "A")
summary_graph <- patchwork::wrap_plots(graph_plots, ncol = 1)
ggsave(
  file.path(output_dir, "phenology-average-transition-sem-both-species-graph.png"),
  combined_graph,
  width = 9,
  height = 11.2,
  dpi = 300
)

heatmap_specs <- sem_heatmap_specs()
shared_heatmap_limit <- compute_shared_sem_heatmap_limit(matrix_data)
treat_labels <- sem_heatmap_treatment_labels()
resp_labels <- ALINV_RESPONSE_LABELS
heatmap_plots <- list()

for (spec_i in seq_len(nrow(heatmap_specs))) {
  spec <- heatmap_specs[spec_i, ]
  species_plots <- lapply(species_vec, function(species_i) {
    panel_df <- prepare_sem_heatmap_panel(
      matrix_df = matrix_data %>% dplyr::filter(.data$species == species_i),
      path_type = spec$path_type,
      resp_labels = resp_labels,
      treat_labels = treat_labels,
      resp_var_order = resp_var
    )

    plot_sem_heatmap_panel(
      panel_df = panel_df,
      limit = shared_heatmap_limit,
      title = paste0(
        stringr::str_to_title(species_i),
        ": ",
        spec$panel_title
      )
    )
  })
  names(species_plots) <- species_vec
  heatmap_plots[[spec$file_stub]] <- species_plots

  for (species_i in species_vec) {
    ggsave(
      file.path(output_dir, paste0("phenology-average-transition-sem-", species_i, "-heatmap-", spec$file_stub, ".png")),
      species_plots[[species_i]],
      width = 9,
      height = 5.6,
      dpi = 300
    )
  }

  combined_heatmap <- (
    patchwork::wrap_plots(species_plots, ncol = 1) +
      patchwork::plot_layout(guides = "collect") +
      patchwork::plot_annotation(
        title = paste0("Phenology average-transition SEM: ", spec$panel_title),
        tag_levels = "A"
      )
  ) &
    ggplot2::theme(legend.position = "right")

  ggsave(
    file.path(output_dir, paste0("phenology-average-transition-sem-both-species-heatmap-", spec$file_stub, ".png")),
    combined_heatmap,
    width = 9,
    height = 8,
    dpi = 300
  )
}

for (path_type_i in c("direct", "indirect", "total")) {
  spec_stub <- heatmap_specs$file_stub[match(path_type_i, heatmap_specs$path_type)]

  for (species_i in species_vec) {
    ggsave(
      file.path(output_dir, paste0("phenology-average-transition-sem-", species_i, "-heatmap-", path_type_i, ".png")),
      heatmap_plots[[spec_stub]][[species_i]],
      width = 9,
      height = 5.6,
      dpi = 300
    )
  }

  combined_alias <- (
    patchwork::wrap_plots(heatmap_plots[[spec_stub]], ncol = 1) +
      patchwork::plot_layout(guides = "collect") +
      patchwork::plot_annotation(
        title = paste0(
          "Phenology average-transition SEM: ",
          heatmap_specs$panel_title[match(path_type_i, heatmap_specs$path_type)]
        ),
        tag_levels = "A"
      )
  ) &
    ggplot2::theme(legend.position = "right")

  ggsave(
    file.path(output_dir, paste0("phenology-average-transition-sem-both-species-heatmap-", path_type_i, ".png")),
    combined_alias,
    width = 9,
    height = 8,
    dpi = 300
  )
}

three_panel_specs <- tibble::tribble(
  ~path_type, ~panel_title,
  "direct", "Direct effects",
  "indirect", "Indirect effects",
  "treatment_to_swc", "SWC effects"
)

four_panel_specs <- tibble::tribble(
  ~path_type, ~panel_title,
  "direct", "Direct effects",
  "indirect", "Indirect effects",
  "total", "Total effects",
  "treatment_to_swc", "SWC effects"
)

make_path_heatmap <- function(species_i, path_type_i, panel_title_i) {
  panel_df <- prepare_sem_heatmap_panel(
    matrix_df = matrix_data %>% dplyr::filter(.data$species == species_i),
    path_type = path_type_i,
    resp_labels = resp_labels,
    treat_labels = treat_labels,
    resp_var_order = resp_var
  )

  plot_sem_heatmap_panel(
    panel_df = panel_df,
    limit = shared_heatmap_limit,
    title = paste0(stringr::str_to_title(species_i), ": ", panel_title_i)
  ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      plot.margin = ggplot2::margin(t = 8, r = 22, b = 10, l = 22),
      axis.text.x = ggplot2::element_text(margin = ggplot2::margin(t = 5))
    )
}

three_panel_plots <- list()
for (species_i in species_vec) {
  species_three_panel_plots <- lapply(seq_len(nrow(three_panel_specs)), function(i) {
    make_path_heatmap(
      species_i = species_i,
      path_type_i = three_panel_specs$path_type[[i]],
      panel_title_i = three_panel_specs$panel_title[[i]]
    )
  })

  three_panel_plots[[species_i]] <- species_three_panel_plots

  species_three_panel <- (
    patchwork::wrap_plots(species_three_panel_plots, ncol = 3) +
      patchwork::plot_layout(guides = "collect") +
      patchwork::plot_annotation(
        title = paste0(
          "Phenology average-transition SEM heatmap: ",
          stringr::str_to_title(species_i)
        )
      )
  ) &
    ggplot2::theme(legend.position = "right")

  ggsave(
    file.path(output_dir, paste0("phenology-average-transition-sem-", species_i, "-heatmap-three-panel.png")),
    species_three_panel,
    width = 13,
    height = 5.2,
    dpi = 300
  )
}

combined_three_panel <- (
  patchwork::wrap_plots(unlist(three_panel_plots, recursive = FALSE), ncol = 3) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(
      title = "Phenology average-transition SEM heatmap"
    )
) &
  ggplot2::theme(legend.position = "right")

ggsave(
  file.path(output_dir, "phenology-average-transition-sem-both-species-heatmap-three-panel.png"),
  combined_three_panel,
  width = 13,
  height = 8.2,
  dpi = 300
)

four_panel_plots <- list()
for (species_i in species_vec) {
  species_four_panel_plots <- lapply(seq_len(nrow(four_panel_specs)), function(i) {
    make_path_heatmap(
      species_i = species_i,
      path_type_i = four_panel_specs$path_type[[i]],
      panel_title_i = four_panel_specs$panel_title[[i]]
    )
  })

  four_panel_plots[[species_i]] <- species_four_panel_plots

  species_four_panel <- (
    patchwork::wrap_plots(species_four_panel_plots, ncol = 4) +
      patchwork::plot_layout(guides = "collect") +
      patchwork::plot_annotation(
        title = paste0(
          "Phenology average-transition SEM heatmap: ",
          stringr::str_to_title(species_i)
        )
      )
  ) &
    ggplot2::theme(legend.position = "right")

  ggsave(
    file.path(output_dir, paste0("phenology-average-transition-sem-", species_i, "-heatmap-four-panel.png")),
    species_four_panel,
    width = 16.5,
    height = 5.2,
    dpi = 300
  )
}

combined_four_panel <- (
  patchwork::wrap_plots(unlist(four_panel_plots, recursive = FALSE), ncol = 4) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(
      title = "Phenology average-transition SEM heatmap"
    )
) &
  ggplot2::theme(legend.position = "right")

ggsave(
  file.path(output_dir, "phenology-average-transition-sem-both-species-heatmap-four-panel.png"),
  combined_four_panel,
  width = 17,
  height = 8.2,
  dpi = 300
)

original_phenology_plot <- plot_ts_phenology(
  style = "band",
  filter_soiltype = soil_type,
  save_fig = FALSE
) +
  ggplot2::labs(
    title = "Original phenology progression assessment",
    subtitle = "Mean transition DOY by treatment group; intervals are mean ± SE."
  ) +
  ggplot2::theme(plot.margin = ggplot2::margin(t = 8, r = 16, b = 8, l = 16))

transition_results <- purrr::map(
  species_vec,
  ~ run_phenology_transition_analysis(
    species_keep = .x,
    soil_type = soil_type,
    include_soil_treatment = include_soil_treatment,
    stages_keep = stages_keep,
    force_run = FALSE
  )
)
transition_effects <- dplyr::bind_rows(
  purrr::map(transition_results, "effects")
)
transition_model_plot <- plot_phenology_transition_effects(transition_effects) +
  patchwork::plot_annotation(
    title = "Transition-timing model estimates",
    subtitle = paste(
      "This is the notebook's transition-timing model, not the SEM:",
      "points are fixed-effect treatment shifts in DOY with 95% Wald intervals."
    )
  )

summary_plot <- (
  patchwork::wrap_elements(full = original_phenology_plot) /
    patchwork::wrap_elements(full = transition_model_plot) /
    patchwork::wrap_elements(full = combined_four_panel) /
    patchwork::wrap_elements(full = summary_graph)
) +
  patchwork::plot_layout(heights = c(1.05, 1.25, 1.25, 1.85)) +
  patchwork::plot_annotation(
    title = "Phenology treatment assessment and average-transition SEM",
    subtitle = paste(
      "Progression and transition-timing panels show the original phenology assessment.",
      "Heatmap and arrows show the average-transition SEM decomposition."
    ),
    tag_levels = "A"
  ) &
  ggplot2::theme(plot.margin = ggplot2::margin(t = 10, r = 12, b = 10, l = 12))

summary_plot_file <- file.path(output_dir, "phenology-average-transition-sem-summary.png")
ggsave(
  summary_plot_file,
  summary_plot,
  width = 18,
  height = 32,
  dpi = 300,
  limitsize = FALSE
)

print(
  sem_effects %>%
    dplyr::select(
      "species",
      "treatment_label",
      "direct",
      "p_direct",
      "indirect",
      "p_indirect",
      "total",
      "p_total",
      "n_trees",
      "n_boxes"
    ),
  n = Inf
)
message("Saved phenology average-transition SEM effects: ", effects_file)
message("Saved phenology average-transition SEM matrix data: ", matrix_file)
message("Saved SEM graph and heatmap PNGs in: ", output_dir)
