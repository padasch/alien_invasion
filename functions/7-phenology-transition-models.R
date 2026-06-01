suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(lme4)
})

phenology_transition_treatment_info <- function(include_soil_treatment = NULL,
                                                soil_filter = NULL) {
  alinv_treatment_config(
    include_soil_treatment = include_soil_treatment,
    soil_filter = soil_filter
  ) %>%
    dplyr::filter(.data$effect != "extreme_event") %>%
    dplyr::transmute(
      effect = .data$effect,
      baseline_level = .data$baseline_level,
      treatment_level = .data$treatment_level,
      treatment_label = .data$short_label,
      contrast_label = .data$contrast_label,
      plot_order = .data$plot_order
    )
}

phenology_transition_cache_path <- function(species_keep,
                                            soil_type = "both",
                                            include_soil_treatment = NULL) {
  include_soil_treatment <- alinv_resolve_include_soil_treatment(
    include_soil_treatment = include_soil_treatment,
    soil_filter = soil_type
  )

  out_dir <- alinv_data_path("model-phenology-transition", create_dir = TRUE)
  soil_mode_tag <- alinv_soil_mode_tag(
    soil_filter = soil_type,
    include_soil_treatment = include_soil_treatment
  )

  file.path(
    out_dir,
    paste0("phenology-transition-", species_keep, "-", soil_mode_tag, ".rds")
  )
}

write_phenology_transition_csv_bundle <- function(result,
                                                  export_stem) {
  result <- result %||% list()

  readr::write_csv(result$data %||% tibble::tibble(), paste0(export_stem, "-data.csv"))
  readr::write_csv(result$effects %||% tibble::tibble(), paste0(export_stem, "-effects.csv"))
}

prepare_phenology_transition_data <- function(species_keep,
                                              soil_type = "both",
                                              include_soil_treatment = NULL,
                                              stages_keep = 2:4) {
  include_soil_treatment <- alinv_resolve_include_soil_treatment(
    include_soil_treatment = include_soil_treatment,
    soil_filter = soil_type
  )

  get_data("tree", "phenology_transitions") %>%
    dplyr::filter(.data$species == species_keep) %>%
    alinv_apply_soil_context(
      soil_filter = soil_type,
      include_soil_treatment = include_soil_treatment
    ) %>%
    dplyr::filter(
      .data$stage %in% stages_keep,
      !is.na(.data$doy),
      !is.na(.data$stage_date)
    ) %>%
    dplyr::mutate(
      tree_id = factor(.data$tree_id),
      boxlabel = factor(.data$boxlabel),
      stage = as.integer(.data$stage),
      stage_label = factor(
        paste0("Stage ", .data$stage),
        levels = paste0("Stage ", stages_keep)
      ),
      precipitation = factor(.data$precipitation, levels = alinv_factor_levels("precipitation")),
      robinia = factor(.data$robinia, levels = alinv_factor_levels("robinia")),
      culture = factor(.data$culture, levels = alinv_factor_levels("culture")),
      soiltype = alinv_relevel_soiltype(.data$soiltype)
    ) %>%
    droplevels()
}

phenology_transition_terms <- function(df,
                                       include_soil_treatment = NULL,
                                       soil_filter = NULL) {
  include_soil_treatment <- alinv_resolve_include_soil_treatment(
    include_soil_treatment = include_soil_treatment,
    soil_filter = soil_filter
  )

  terms <- c("precipitation", "robinia", "culture")
  if (isTRUE(include_soil_treatment)) {
    terms <- c(terms, "soiltype")
  }

  terms[vapply(
    terms,
    function(term_i) {
      term_i %in% names(df) && dplyr::n_distinct(df[[term_i]], na.rm = TRUE) > 1
    },
    logical(1)
  )]
}

fit_phenology_transition_stage_model <- function(df_stage,
                                                 include_soil_treatment = NULL,
                                                 soil_filter = NULL) {
  rhs_terms <- phenology_transition_terms(
    df_stage,
    include_soil_treatment = include_soil_treatment,
    soil_filter = soil_filter
  )

  rhs <- if (length(rhs_terms)) paste(rhs_terms, collapse = " + ") else "1"
  fml <- stats::as.formula(paste("doy ~", rhs, "+ (1 | boxlabel)"))

  lme4::lmer(fml, data = df_stage, REML = TRUE)
}

fit_phenology_transition_pooled_model <- function(df_all,
                                                  include_soil_treatment = NULL,
                                                  soil_filter = NULL) {
  rhs_terms <- phenology_transition_terms(
    df_all,
    include_soil_treatment = include_soil_treatment,
    soil_filter = soil_filter
  )

  rhs <- c("stage_label", rhs_terms)
  fml <- stats::as.formula(
    paste("doy ~", paste(rhs, collapse = " + "), "+ (1 | boxlabel) + (1 | tree_id)")
  )

  lme4::lmer(fml, data = df_all, REML = TRUE)
}

extract_phenology_transition_effects <- function(fit,
                                                 stage_label,
                                                 species_keep,
                                                 include_soil_treatment = NULL,
                                                 soil_filter = NULL,
                                                 model_type = c("stage_specific", "average")) {
  model_type <- match.arg(model_type)
  trt_info <- phenology_transition_treatment_info(
    include_soil_treatment = include_soil_treatment,
    soil_filter = soil_filter
  )
  beta <- lme4::fixef(fit)
  vc <- as.matrix(stats::vcov(fit))

  out <- lapply(seq_len(nrow(trt_info)), function(i) {
    info_i <- trt_info[i, ]
    effect_i <- info_i$effect[[1]]
    coef_name <- paste0(effect_i, info_i$treatment_level[[1]])

    if (!coef_name %in% names(beta)) {
      return(NULL)
    }
    est <- unname(beta[[coef_name]])
    se <- sqrt(unname(vc[coef_name, coef_name]))
    lower <- est - 1.96 * se
    upper <- est + 1.96 * se
    stat <- est / se
    p_value <- 2 * stats::pnorm(-abs(stat))

    tibble::tibble(
      species = species_keep,
      model_type = model_type,
      stage_label = stage_label,
      effect = effect_i,
      treatment_label = info_i$treatment_label[[1]],
      contrast_label = info_i$contrast_label[[1]],
      estimate = est,
      se = se,
      lower = lower,
      upper = upper,
      p_value = p_value,
      earlier_later = dplyr::case_when(
        .data$estimate < 0 ~ "earlier",
        .data$estimate > 0 ~ "later",
        TRUE ~ "no change"
      ),
      abs_days = abs(.data$estimate),
      plot_order = info_i$plot_order[[1]]
    )
  })

  dplyr::bind_rows(out)
}

run_phenology_transition_analysis <- function(species_keep,
                                              soil_type = "both",
                                              include_soil_treatment = NULL,
                                              stages_keep = 2:4,
                                              force_run = FALSE) {
  include_soil_treatment <- alinv_resolve_include_soil_treatment(
    include_soil_treatment = include_soil_treatment,
    soil_filter = soil_type
  )

  cache_path <- phenology_transition_cache_path(
    species_keep = species_keep,
    soil_type = soil_type,
    include_soil_treatment = include_soil_treatment
  )
  export_stem <- tools::file_path_sans_ext(cache_path)

  if (file.exists(cache_path) && !force_run) {
    result <- readRDS(cache_path)
    write_phenology_transition_csv_bundle(result, export_stem = export_stem)
    return(result)
  }

  df <- prepare_phenology_transition_data(
    species_keep = species_keep,
    soil_type = soil_type,
    include_soil_treatment = include_soil_treatment,
    stages_keep = stages_keep
  )

  stage_models <- lapply(stages_keep, function(stage_i) {
    df_stage <- df %>% dplyr::filter(.data$stage == stage_i)
    fit_phenology_transition_stage_model(
      df_stage,
      include_soil_treatment = include_soil_treatment,
      soil_filter = soil_type
    )
  })
  names(stage_models) <- paste0("Stage ", stages_keep)

  pooled_model <- fit_phenology_transition_pooled_model(
    df,
    include_soil_treatment = include_soil_treatment,
    soil_filter = soil_type
  )

  effects_stage <- dplyr::bind_rows(
    lapply(seq_along(stage_models), function(i) {
      extract_phenology_transition_effects(
        fit = stage_models[[i]],
        stage_label = names(stage_models)[[i]],
        species_keep = species_keep,
        include_soil_treatment = include_soil_treatment,
        soil_filter = soil_type,
        model_type = "stage_specific"
      )
    })
  )

  effects_avg <- extract_phenology_transition_effects(
    fit = pooled_model,
    stage_label = "Average (Stages 2-4)",
    species_keep = species_keep,
    include_soil_treatment = include_soil_treatment,
    soil_filter = soil_type,
    model_type = "average"
  )

  result <- list(
    data = df,
    stage_models = stage_models,
    pooled_model = pooled_model,
    effects = dplyr::bind_rows(effects_stage, effects_avg)
  )

  saveRDS(result, cache_path)
  write_phenology_transition_csv_bundle(result, export_stem = export_stem)

  result
}

plot_phenology_transition_effects <- function(effects_df) {
  if (is.null(effects_df) || !nrow(effects_df)) {
    return(
      ggplot() +
        theme_void() +
        ggtitle("Phenology transition timing") +
        labs(subtitle = "No phenology transition effects were available.")
    )
  }

  effects_plot <- effects_df %>%
    dplyr::mutate(
      stage_label = factor(
        .data$stage_label,
        levels = c("Stage 2", "Stage 3", "Stage 4", "Average (Stages 2-4)")
      ),
      treatment_label = factor(
        .data$treatment_label,
        levels = c("Drought", "With robinia", "Mixed culture", "Drier soil")
      ),
      species = factor(.data$species, levels = c("fagus", "quercus"))
    ) %>%
    dplyr::arrange(.data$plot_order)

  if (!requireNamespace("patchwork", quietly = TRUE)) {
    warning("Package 'patchwork' not available; using a single faceted phenology transition plot.")
    return(
      ggplot(effects_plot, aes(x = .data$estimate, y = .data$treatment_label, color = .data$treatment_label)) +
        geom_vline(xintercept = 0, linetype = 2, color = "grey50", linewidth = 0.5) +
        geom_segment(
          aes(
            x = .data$lower,
            xend = .data$upper,
            y = .data$treatment_label,
            yend = .data$treatment_label
          ),
          linewidth = 0.85,
          alpha = 0.8
        ) +
        geom_point(size = 2.4) +
        facet_grid(
          species ~ stage_label,
          labeller = ggplot2::labeller(
            species = c(fagus = "Fagus", quercus = "Quercus")
          )
        ) +
        scale_color_brewer(palette = "Dark2", drop = FALSE) +
        labs(
          title = "Phenology transition timing by treatment",
          subtitle = "Response = DOY of reaching each stage. Negative values mean the named treatment reaches the stage earlier; positive values mean later.",
          x = "Shift in transition DOY under treatment (days)",
          y = NULL
        ) +
        theme_classic(base_size = 11) +
        theme(
          legend.position = "none",
          strip.background = ggplot2::element_rect(fill = "black", color = "black"),
          strip.text = ggplot2::element_text(color = "white", face = "bold")
        )
    )
  }

  x_limits <- range(c(effects_plot$lower, effects_plot$upper, 0), na.rm = TRUE)
  x_pad <- diff(x_limits) * 0.08
  if (!is.finite(x_pad) || x_pad <= 0) {
    x_pad <- 0.5
  }
  x_limits <- c(x_limits[[1]] - x_pad, x_limits[[2]] + x_pad)

  species_labels <- c(fagus = "Fagus", quercus = "Quercus")

  build_species_panel <- function(species_key) {
    ggplot(
      dplyr::filter(effects_plot, .data$species == species_key),
      aes(x = .data$estimate, y = .data$treatment_label, color = .data$treatment_label)
    ) +
      geom_vline(xintercept = 0, linetype = 2, color = "grey50", linewidth = 0.5) +
      geom_segment(
        aes(
          x = .data$lower,
          xend = .data$upper,
          y = .data$treatment_label,
          yend = .data$treatment_label
        ),
        linewidth = 0.85,
        alpha = 0.8
      ) +
      geom_point(size = 2.4) +
      facet_wrap(~stage_label, nrow = 1) +
      scale_color_brewer(palette = "Dark2", drop = FALSE) +
      coord_cartesian(xlim = x_limits) +
      labs(
        title = species_labels[[species_key]] %||% species_key,
        x = "Shift in transition DOY under treatment (days)",
        y = NULL
      ) +
      theme_classic(base_size = 11) +
      theme(
        legend.position = "none",
        plot.title = ggplot2::element_text(face = "bold", hjust = 0),
        strip.background = ggplot2::element_rect(fill = "black", color = "black"),
        strip.text = ggplot2::element_text(color = "white", face = "bold")
      )
  }

  species_present <- levels(droplevels(effects_plot$species))
  species_plots <- lapply(species_present, build_species_panel)

  patchwork::wrap_plots(species_plots, ncol = 1) +
    patchwork::plot_annotation(
      title = "Phenology transition timing by treatment",
      subtitle = "Response = DOY of reaching each stage. Negative values mean the named treatment reaches the stage earlier; positive values mean later."
    )
}
