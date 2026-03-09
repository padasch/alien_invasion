# 0-biomass.R
# Biomass wrangling, plotting, and GLMM functions

wrangle_tree_biomass <- function(fp, sheet = "Biomass") {
  meta_tree <- get_meta("tree") %>%
    mutate(tree_id = as.character(tree_id)) %>%
    select(treelabel, tree_id, boxlabel, soiltype, culture, robinia) %>%
    distinct(treelabel, .keep_all = TRUE)

  read_excel(fp, sheet = sheet) %>%
    remove_empty(which = "cols") %>%
    clean_names() %>%
    mutate(across(where(is.character), ~ tolower(trimws(.x)))) %>%
    mutate(across(where(is.character), ~ gsub("_", "-", .x))) %>%
    mutate(
      treelabel = sub("^(fagus|quercus)-", "", label),
      root_biomass = parse_number(as.character(root_biomass)),
      shoot_biomass = parse_number(as.character(shoot_biomass))
    ) %>%
    mutate(
      root_shoot_biomass = if_else(
        !is.na(shoot_biomass) & shoot_biomass > 0,
        root_biomass / shoot_biomass,
        NA_real_
      )
    ) %>%
    left_join(meta_tree, by = "treelabel") %>%
    select(
      tree_id, treelabel, label, species,
      precipitation, soiltype, culture, robinia, boxlabel,
      compartment, root_bag, shoot_bag,
      root_biomass_with_bag, shoot_biomass_with_bag,
      root_biomass, shoot_biomass, root_shoot_biomass,
      comment
    ) %>%
    arrange(species, soiltype, robinia, precipitation, culture, tree_id)
}

plot_tree_biomass_treatments <- function(df_biomass, species_keep = "fagus") {
  stopifnot(all(c("root_biomass", "shoot_biomass", "root_shoot_biomass") %in% names(df_biomass)))

  df_base <- df_biomass %>%
    filter(species == species_keep) %>%
    filter(!is.na(precipitation), !is.na(soiltype), !is.na(culture), !is.na(robinia)) %>%
    mutate(
      precipitation = factor(precipitation, levels = c("control", "drought")),
      culture = factor(culture, levels = c("mono", "mixed")),
      robinia = factor(robinia, levels = c("without-robinia", "with-robinia")),
      soiltype = factor(soiltype, levels = c("inoc-beech", "inoc-robinia"))
    )

  robinia_labels <- c(
    `without-robinia` = "Without robinia",
    `with-robinia` = "With robinia"
  )

  soiltype_labels <- c(
    `inoc-beech` = "Beech Soil (was drier)",
    `inoc-robinia` = "Robinia Soil (was wetter)"
  )

  metric_info <- list(
    root_biomass       = list(col = "root_biomass",       title = "Root biomass (g)"),
    shoot_biomass      = list(col = "shoot_biomass",      title = "Shoot biomass (g)"),
    root_shoot_biomass = list(col = "root_shoot_biomass", title = "Root:shoot biomass (-)")
  )

  make_panel <- function(info) {
    df_m <- df_base %>%
      filter(!is.na(.data[[info$col]])) %>%
      rename(value = !!info$col)

    ggplot2::ggplot(
      df_m,
      ggplot2::aes(
        x = culture, y = value,
        fill = precipitation, color = precipitation, alpha = culture
      )
    ) +
      ggplot2::geom_boxplot(
        outlier.shape = NA, width = 0.65,
        position = ggplot2::position_dodge2(width = 0.75, preserve = "single")
      ) +
      ggplot2::geom_point(
        position = ggplot2::position_jitterdodge(jitter.width = 0.12, dodge.width = 0.75),
        size = 1.6, stroke = 0
      ) +
      ggplot2::facet_grid(
        soiltype ~ robinia,
        labeller = ggplot2::labeller(
          robinia = robinia_labels,
          soiltype = soiltype_labels
        )
      ) +
      ggplot2::scale_fill_manual(values = c(control = "grey50", drought = "indianred"), drop = FALSE) +
      ggplot2::scale_color_manual(values = c(control = "grey50", drought = "indianred"), drop = FALSE) +
      ggplot2::scale_alpha_manual(values = c(mono = 0.55, mixed = 0.9), guide = "none") +
      ggplot2::labs(x = "Culture", y = info$title, fill = "Precipitation", color = "Precipitation") +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(
        strip.background = ggplot2::element_rect(fill = "black", color = "black"),
        strip.text = ggplot2::element_text(color = "white", face = "bold"),
        legend.position = "bottom"
      )
  }

  panels <- lapply(metric_info, make_panel)

  combined <- patchwork::wrap_plots(panels, ncol = 1) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(
      title = paste("Tree biomass by treatments -", species_keep),
      theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
    ) &
    ggplot2::theme(legend.position = "bottom")

  combined
}

plot_tree_biomass_treatments_all_species <- function(df_biomass, species_vec = c("fagus", "quercus")) {
  species_vec <- unique(species_vec)
  setNames(
    lapply(species_vec, function(sp) plot_tree_biomass_treatments(df_biomass, species_keep = sp)),
    species_vec
  )
}

fit_biomass_glmm <- function(df_biomass, species_keep = "fagus", metric = "root_biomass") {
  df_m <- df_biomass %>%
    filter(species == species_keep, !is.na(.data[[metric]])) %>%
    mutate(
      precipitation = factor(precipitation, levels = c("control", "drought")),
      culture       = factor(culture,       levels = c("mono", "mixed")),
      robinia       = factor(robinia,       levels = c("without-robinia", "with-robinia")),
      soiltype_f    = factor(soiltype,      levels = c("inoc-beech", "inoc-robinia"))
    ) %>%
    rename(y = !!metric) %>%
    mutate(y_z = as.numeric(scale(y)))

  mod <- lme4::lmer(
    y_z ~ precipitation + culture + robinia + soiltype_f + (1 | boxlabel),
    data = df_m,
    REML = TRUE
  )

  list(model = mod, data = df_m, metric = metric, species = species_keep)
}

extract_biomass_effects <- function(fit_obj) {
  mod <- fit_obj$model
  cf  <- summary(mod)$coefficients
  df_eff <- tibble::tibble(
    term     = rownames(cf),
    estimate = cf[, "Estimate"],
    se       = cf[, "Std. Error"],
    t_value  = cf[, "t value"]
  ) %>%
    filter(term != "(Intercept)") %>%
    mutate(
      ci_lo = estimate - 1.96 * se,
      ci_hi = estimate + 1.96 * se,
      significant = !(ci_lo <= 0 & ci_hi >= 0)
    )

  df_eff$metric  <- fit_obj$metric
  df_eff$species <- fit_obj$species
  df_eff
}

# Readable labels for GLMM coefficient names
.rename_glmm_terms <- function(terms) {
  lookup <- c(
    "precipitationdrought"  = "Precipitation: Control \u2192 Drought",
    "culturemixed"          = "Culture: Mono \u2192 Mixed",
    "robiniawith-robinia"   = "Robinia: Without \u2192 With",
    "soiltype_finoc-robinia" = "Soil: Beech Soil \u2192 Robinia Soil"
  )
  ifelse(terms %in% names(lookup), lookup[terms], terms)
}

extract_model_performance <- function(fit_obj) {

  mod <- fit_obj$model
  r2  <- tryCatch(MuMIn::r.squaredGLMM(mod), error = function(e) matrix(c(NA, NA), nrow = 1, dimnames = list(NULL, c("R2m", "R2c"))))
  tibble::tibble(
    metric   = fit_obj$metric,
    species  = fit_obj$species,
    AIC      = stats::AIC(mod),
    BIC      = stats::BIC(mod),
    R2_marginal    = r2[1, "R2m"],
    R2_conditional = r2[1, "R2c"],
    n_obs    = nrow(fit_obj$data)
  )
}

plot_biomass_effects <- function(df_effects) {
  metric_labels <- c(
    root_biomass       = "Root biomass",
    shoot_biomass      = "Shoot biomass",
    root_shoot_biomass = "Root:shoot biomass"
  )

  df_effects <- df_effects %>%
    mutate(
      term_label   = .rename_glmm_terms(term),
      metric_label = metric_labels[metric]
    )

  ggplot2::ggplot(df_effects, ggplot2::aes(
    x = term_label, y = estimate, color = significant
  )) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_pointrange(
      ggplot2::aes(ymin = ci_lo, ymax = ci_hi),
      size = 0.6, linewidth = 0.7, fatten = 2.5
    ) +
    ggplot2::facet_wrap(~ metric_label, ncol = 1, scales = "free_x") +
    ggplot2::coord_flip() +
    ggplot2::scale_color_manual(
      values = c(`TRUE` = "indianred", `FALSE` = "grey50"),
      labels = c(`TRUE` = "95% CI excludes 0", `FALSE` = "95% CI includes 0"),
      name = NULL
    ) +
    ggplot2::labs(
      title = paste("Treatment effects on biomass \u2013", df_effects$species[1]),
      x = NULL, y = "Standardized effect (SD units, \u00b195% CI)"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "black", color = "black"),
      strip.text = ggplot2::element_text(color = "white", face = "bold"),
      legend.position = "bottom"
    )
}

run_biomass_glmm_species <- function(df_biomass, species_keep = "fagus",
                                     metrics = c("root_biomass", "shoot_biomass", "root_shoot_biomass")) {
  fits <- lapply(metrics, function(m) fit_biomass_glmm(df_biomass, species_keep, m))
  effects <- do.call(rbind, lapply(fits, extract_biomass_effects))
  perf    <- do.call(rbind, lapply(fits, extract_model_performance))
  p <- plot_biomass_effects(effects)
  list(fits = fits, effects = effects, performance = perf, plot = p)
}
