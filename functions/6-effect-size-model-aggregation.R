library(patchwork)

collect_all_effects <- function(
    var_grid      = get_data_var_grid(),
    species_vec   = c("fagus", "quercus"),
    soil_type     = "both",
    add_covars    = FALSE,
    covars_fun    = NULL
) {
  # needs: dplyr, purrr, tibble loaded
  purrr::pmap_dfr(
    var_grid,
    function(type, data, resp_var) {
      purrr::map_dfr(
        species_vec,
        function(sp) {
          res <- tryCatch(
            make_effect_figure_generic(
              type           = type,
              data_name      = data,
              resp_var       = resp_var,
              target_species = sp,
              soil_type      = soil_type,
              add_covars     = add_covars,
              covars_fun     = covars_fun
            ),
            error = function(e) NULL
          )
          
          if (is.null(res) || is.null(res$effects) || !nrow(res$effects)) {
            return(NULL)
          }
          
          res$effects %>%
            dplyr::mutate(
              type     = type,
              data     = data,
              resp_var = resp_var,
              species  = sp
            )
        }
      )
    }
  )
}

add_importance_metric <- function(
    effects_df,
    use_stat      = TRUE,   # if FALSE, uses |estimate| instead
    scale_within  = TRUE    # scale within type–data–resp_var–species–effect
) {
  # decide globally whether to use stat
  use_stat_now <- use_stat &&
    ("stat" %in% names(effects_df)) &&
    !all(is.na(effects_df$stat))
  
  df <- if (use_stat_now) {
    effects_df %>%
      dplyr::mutate(importance_raw = abs(.data$stat))
  } else {
    effects_df %>%
      dplyr::mutate(importance_raw = abs(.data$estimate))
  }
  
  if (!scale_within) {
    return(df %>% dplyr::mutate(importance = importance_raw))
  }
  
  df %>%
    dplyr::group_by(type, data, resp_var, species, effect) %>%
    dplyr::mutate(
      sd_imp     = stats::sd(importance_raw, na.rm = TRUE),
      sd_imp     = dplyr::if_else(is.na(sd_imp) | sd_imp == 0, 1, sd_imp),
      importance = importance_raw / sd_imp
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-sd_imp)
}

plot_treatment_importance_phases <- function(
    effects_imp_phase,
    which_effects = c("robinia", "precipitation", "culture", "soiltype"),
    stat          = c(
      "boxplot",
      "mean_sd", "mean_se", "mean_ci",
      "pointrange_sd", "pointrange_se", "pointrange_ci"
    ),
    # ylab          = "Standardized importance (|test statistic|)"
    ylab          = "Importance (Standardized Effect Size)"
) {
  stat <- match.arg(stat)
  
  df_plot <- effects_imp_phase %>%
    dplyr::filter(.data$effect %in% which_effects) %>%
    dplyr::mutate(effect = factor(.data$effect, levels = which_effects))
  
  if (!nrow(df_plot)) stop("No effects to plot for the selected treatments.")
  
  # ---------------------------------------------------------------
  # Phase-level info
  # ---------------------------------------------------------------
  phase_counts <- df_plot %>%
    dplyr::group_by(phase) %>%
    dplyr::summarise(
      n_estimates = dplyr::n(),
      .groups = "drop"
    )
  
  phase_data_info <- df_plot %>%
    dplyr::group_by(phase) %>%
    dplyr::summarise(
      n_sources = dplyr::n_distinct(paste(type, data, data, sep = "_")),
      sources   = paste(sort(unique(data)), collapse = ", "),
      .groups   = "drop"
    )
  
  phase_caption <- dplyr::left_join(phase_counts, phase_data_info, by = "phase") %>%
    dplyr::mutate(
      lbl = paste0(
        phase, " (n = ", n_estimates,
        ", data: ", n_sources,
        " [", sources, "])"
      )
    ) %>%
    dplyr::pull(lbl) %>%
    paste(collapse = "\n")
  
  # ============================================================
  # BOXplot option
  # ============================================================
  if (stat == "boxplot") {
    return(
      ggplot2::ggplot(
        df_plot,
        ggplot2::aes(x = phase, y = importance, fill = effect, color = effect)
      ) +
        ggplot2::geom_hline(
          yintercept = c(1, 2),
          linetype = "dotted",
          alpha=0.5,
          color = "grey60",
          linewidth=0.5
        ) +
        ggplot2::geom_boxplot(
          position = ggplot2::position_dodge(width = 0.8),
          width = 0.7,
          outlier.size = 0.6
        ) +
        ggplot2::geom_hline(
          yintercept = c(1, 2),
          linetype = "dotted",
          alpha=0.5,
          color = "grey60",
          linewidth=0.5
        ) +
        ggplot2::scale_fill_brewer(palette = "Dark2") +
        ggplot2::scale_color_brewer(palette = "Dark2") +
        ggplot2::labs(
          x = NULL,
          y = ylab,
          fill = "Treatment",
          color = "Treatment",
          caption = paste(phase_caption)
        ) +
        ggplot2::theme_classic(base_size = 12) +
        ggplot2::theme(legend.position = "bottom")
    )
  }
  
  # ============================================================
  # BAR or POINTRANGE summary stats (SD / SE / CI)
  # ============================================================
  summary_type <- dplyr::case_when(
    stat %in% c("mean_sd",      "pointrange_sd") ~ "sd",
    stat %in% c("mean_se",      "pointrange_se") ~ "se",
    stat %in% c("mean_ci",      "pointrange_ci") ~ "ci"
  )
  
  df_sum <- df_plot %>%
    dplyr::group_by(phase, effect) %>%
    dplyr::summarise(
      mean_imp = mean(importance, na.rm = TRUE),
      sd_imp   = stats::sd(importance, na.rm = TRUE),
      n        = dplyr::n(),
      .groups  = "drop"
    ) %>%
    dplyr::mutate(
      se_imp = sd_imp / sqrt(n),
      ci_low = mean_imp - 1.96 * se_imp,
      ci_up  = mean_imp + 1.96 * se_imp,
      ymin = dplyr::case_when(
        .env$summary_type == "sd" ~ mean_imp - sd_imp,
        .env$summary_type == "se" ~ mean_imp - se_imp,
        .env$summary_type == "ci" ~ ci_low
      ),
      ymax = dplyr::case_when(
        .env$summary_type == "sd" ~ mean_imp + sd_imp,
        .env$summary_type == "se" ~ mean_imp + se_imp,
        .env$summary_type == "ci" ~ ci_up
      )
    )
  
  # ============================================================
  # MEAN ± SD/SE/CI BARS
  # ============================================================
  if (stat %in% c("mean_sd", "mean_se", "mean_ci")) {
    return(
      ggplot2::ggplot(
        df_sum,
        ggplot2::aes(
          x = phase, y = mean_imp,
          ymin = ymin, ymax = ymax,
          fill = effect, color = effect
        )
      ) +
        ggplot2::geom_hline(
          yintercept = c(1, 2),
          linetype = "dotted",
          alpha=0.5,
          color = "grey60",
          linewidth=0.5
        ) +
        ggplot2::geom_col(
          position = ggplot2::position_dodge(width = 0.8),
          width = 0.65,
          alpha = 0.7
        ) +
        ggplot2::geom_errorbar(
          position = ggplot2::position_dodge(width = 0.8),
          width = 0.2,
          linewidth = 0.6
        ) +
        ggplot2::scale_fill_brewer(palette = "Dark2") +
        ggplot2::scale_color_brewer(palette = "Dark2") +
        ggplot2::labs(
          x = NULL,
          y = ylab,
          fill = "Treatment",
          color = "Treatment",
          caption = paste(phase_caption)
        ) +
        ggplot2::theme_classic(base_size = 12) +
        ggplot2::theme(legend.position = "bottom")
    )
  }
}

add_phase_to_effects <- function(effects_imp) {
  effects_imp %>%
    dplyr::mutate(
      month = lubridate::month(.data$date),
      phase = dplyr::case_when(
        month <= 6 ~ "until June",      # Jan–Jun
        month <= 8 ~ "July–August",     # Jul–Aug
        TRUE       ~ "September+"       # Sep–Dec
      ),
      phase = factor(
        phase,
        levels = c("until June", "July–August", "September+")
      )
    ) %>%
    dplyr::select(-month)
}

# ============================================================
# EXAMPLE CODE: Uncomment and run to generate effect summaries
# ============================================================
# all_eff <- collect_all_effects()
# all_eff_imp <- all_eff %>% add_importance_metric(use_stat = FALSE, scale_within = TRUE)
# all_eff_imp_phase <- add_phase_to_effects(all_eff_imp)
#
# mystat <- "mean_se"
# p_all <- plot_treatment_importance_phases(all_eff_imp_phase, stat = mystat)
# p_fag <- plot_treatment_importance_phases(all_eff_imp_phase |> filter(species=="fagus"), stat = mystat)
# p_que <- plot_treatment_importance_phases(all_eff_imp_phase |> filter(species=="quercus"), stat = mystat)
#
# ymin <- 0
# ymax <- 3
#
# p_all_fix <- p_all + coord_cartesian(ylim = c(ymin, ymax)) + ggtitle("All species")
# p_fag_fix <- p_fag + coord_cartesian(ylim = c(ymin, ymax)) + ggtitle("Fagus")
# p_que_fix <- p_que + coord_cartesian(ylim = c(ymin, ymax)) + ggtitle("Quercus")
#
# p_3panel <- (p_all_fix | p_fag_fix | p_que_fix) +
#   plot_layout(guides = "collect") +
#   plot_annotation(tag_levels = "A") &
#   theme(legend.position = "bottom",
#         plot.tag = element_text(face = "bold"))
#
# p_3panel
