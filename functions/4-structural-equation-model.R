# Packages
library(dplyr)

prepare_sem_data <- function(df_prepared,
                             scale_all_numeric = TRUE) {
  df <- df_prepared %>%
    dplyr::filter(
      !is.na(y),
      !is.na(swc),
      !is.na(date)
    ) %>%
    dplyr::mutate(
      # make sure 'date' is a real Date, not factor
      date = as.Date(as.character(date)),

      # factors / treatments
      tree_id = factor(tree_id),
      boxlabel = factor(boxlabel),
      robinia = factor(robinia, levels = c("without-robinia", "with-robinia")),
      extreme_event = factor(extreme_event, levels = c("no", "yes")),
      precipitation = factor(precipitation, levels = c("control", "drought")),
      culture = factor(culture, levels = c("mono", "mixed")),
      soiltype = alinv_relevel_soiltype(soiltype),

      # time components (raw + centered)
      doy = as.numeric(format(date, "%j")),
      doy_c_raw = scale(doy, center = TRUE, scale = FALSE)[, 1],

      # keep originals of continuous variables
      swc_org = swc,
      y_org = y
    )

  if (isTRUE(scale_all_numeric)) {
    # explicitly scale the predictors we actually use
    df <- df %>%
      dplyr::mutate(
        doy_c_org = doy_c_raw,
        swc       = as.numeric(scale(swc_org)),
        y         = as.numeric(scale(y_org)),
        doy_c     = as.numeric(scale(doy_c_raw)),
        doy_c2    = doy_c^2
      )
  } else {
    # no scaling – just use centered doy
    df <- df %>%
      dplyr::mutate(
        doy_c  = doy_c_raw,
        doy_c2 = doy_c_raw^2
      )
  }

  df %>%
    dplyr::select(
      tree_id, boxlabel, date,
      swc, swc_org,
      y, y_org,
      robinia, extreme_event, precipitation, culture, soiltype,
      doy, doy_c, doy_c2
    ) %>%
    tidyr::drop_na()
}

build_sem_rhs <- function(df, include_interaction = TRUE) {
  base_terms <- c(
    "doy_c", "doy_c2",
    "precipitation", "robinia", "soiltype", "culture", "extreme_event"
  )

  # only add interaction if requested AND both factors have >1 level
  if (isTRUE(include_interaction) &&
    length(unique(df$robinia)) > 1 &&
    length(unique(df$extreme_event)) > 1) {
    base_terms <- c(base_terms, "extreme_event:robinia")
  }

  list(
    rhs_swc  = paste(base_terms, collapse = " + "),
    rhs_resp = paste(c("swc", base_terms), collapse = " + ")
  )
}

fit_sem_models <- function(df_sem_ready,
                           include_interaction = TRUE) {
  rhs <- build_sem_rhs(df_sem_ready, include_interaction = include_interaction)

  fml_swc <- as.formula(
    paste("swc ~", rhs$rhs_swc, "+ (1|boxlabel)")
  )

  fml_resp <- as.formula(
    paste("y ~", rhs$rhs_resp, "+ (1|boxlabel) + (1|tree_id)")
  )

  mod_swc <- lme4::lmer(fml_swc, data = df_sem_ready)
  mod_resp <- lme4::lmer(fml_resp, data = df_sem_ready)

  list(
    mod_swc  = mod_swc,
    mod_resp = mod_resp,
    fml_swc  = fml_swc,
    fml_resp = fml_resp
  )
}

extract_sem_effects_with_se <- function(mod_swc,
                                        mod_resp,
                                        factors = c("robinia", "precipitation", "culture", "soiltype", "extreme_event")) {
  beta_swc <- lme4::fixef(mod_swc)
  beta_resp <- lme4::fixef(mod_resp)
  vcov_swc <- as.matrix(vcov(mod_swc))
  vcov_resp <- as.matrix(vcov(mod_resp))

  # SWC → response
  b <- beta_resp["swc"]
  var_b <- vcov_resp["swc", "swc"]
  se_b <- sqrt(var_b)
  p_b <- 2 * pnorm(-abs(b / se_b))

  get_coef_var <- function(beta, V, pattern) {
    nm <- names(beta)
    idx <- grep(pattern, nm, fixed = TRUE)
    if (!length(idx)) {
      return(list(est = NA_real_, se = NA_real_, p = NA_real_))
    }
    est <- beta[idx[1]]
    var <- V[idx[1], idx[1]]
    se <- sqrt(var)
    p <- 2 * pnorm(-abs(est / se))
    list(est = est, se = se, p = p)
  }

  out <- lapply(factors, function(fac) {
    a_info <- get_coef_var(beta_swc, vcov_swc, fac)
    c_info <- get_coef_var(beta_resp, vcov_resp, fac)

    a <- a_info$est
    se_a <- a_info$se
    p_a <- a_info$p

    c_dir <- c_info$est
    se_c <- c_info$se
    p_c <- c_info$p

    indirect <- a * b
    var_ind <- (b^2) * (se_a^2) + (a^2) * (se_b^2)
    se_ind <- sqrt(var_ind)
    p_ind <- 2 * pnorm(-abs(indirect / se_ind))

    total <- c_dir + indirect
    var_total <- (se_c^2) + var_ind
    se_total <- sqrt(var_total)
    p_total <- 2 * pnorm(-abs(total / se_total))

    tibble::tibble(
      factor   = fac,
      a        = a,
      se_a     = se_a,
      p_a      = p_a,
      b        = b,
      se_b     = se_b,
      p_b      = p_b,
      c_direct = c_dir,
      se_c     = se_c,
      p_c      = p_c,
      indirect = indirect,
      se_ind   = se_ind,
      p_ind    = p_ind,
      total    = total,
      se_tot   = se_total,
      p_tot    = p_total
    )
  })

  dplyr::bind_rows(out)
}

extract_interaction_effect <- function(mod_swc,
                                       mod_resp,
                                       treat_factors = c("robinia", "precipitation", "culture", "soiltype", "extreme_event")) {
  beta_swc <- lme4::fixef(mod_swc)
  beta_resp <- lme4::fixef(mod_resp)
  vcov_swc <- as.matrix(vcov(mod_swc))
  vcov_resp <- as.matrix(vcov(mod_resp))

  # SWC → response
  b <- beta_resp["swc"]
  var_b <- vcov_resp["swc", "swc"]
  se_b <- sqrt(var_b)

  out <- list()

  if (length(treat_factors) < 2) {
    return(tibble::tibble(
      factor   = character(),
      a        = double(), se_a   = double(), p_a   = double(),
      b        = double(), se_b   = double(), p_b   = double(),
      c_direct = double(), se_c   = double(), p_c   = double(),
      indirect = double(), se_ind = double(), p_ind = double(),
      total    = double(), se_tot = double(), p_tot = double()
    ))
  }

  for (pair in combn(treat_factors, 2, simplify = FALSE)) {
    f1 <- pair[1]
    f2 <- pair[2]

    # find the unique coefficient that contains both factor names
    idx_swc <- which(grepl(f1, names(beta_swc)) & grepl(f2, names(beta_swc)))
    idx_resp <- which(grepl(f1, names(beta_resp)) & grepl(f2, names(beta_resp)))

    # if we don't have exactly one interaction term in both models, skip this pair
    if (length(idx_swc) != 1L || length(idx_resp) != 1L) next

    a_int <- beta_swc[idx_swc]
    var_a <- vcov_swc[idx_swc, idx_swc]
    se_a <- sqrt(var_a)

    c_int <- beta_resp[idx_resp]
    var_c <- vcov_resp[idx_resp, idx_resp]
    se_c <- sqrt(var_c)

    p_a <- 2 * pnorm(-abs(a_int / se_a))
    p_c <- 2 * pnorm(-abs(c_int / se_c))

    ind <- a_int * b
    var_ind <- (b^2) * var_a + (a_int^2) * var_b
    se_ind <- sqrt(var_ind)
    p_ind <- 2 * pnorm(-abs(ind / se_ind))

    tot <- c_int + ind
    var_tot <- var_c + var_ind
    se_tot <- sqrt(var_tot)
    p_tot <- 2 * pnorm(-abs(tot / se_tot))

    out[[length(out) + 1L]] <- tibble::tibble(
      factor   = paste(f1, f2, sep = ":"),
      a        = as.numeric(a_int),
      se_a     = as.numeric(se_a),
      p_a      = as.numeric(p_a),
      b        = as.numeric(b),
      se_b     = as.numeric(se_b),
      p_b      = NA_real_, # not testing SWC→response here
      c_direct = as.numeric(c_int),
      se_c     = as.numeric(se_c),
      p_c      = as.numeric(p_c),
      indirect = as.numeric(ind),
      se_ind   = as.numeric(se_ind),
      p_ind    = as.numeric(p_ind),
      total    = as.numeric(tot),
      se_tot   = as.numeric(se_tot),
      p_tot    = as.numeric(p_tot)
    )
  }

  if (!length(out)) {
    return(tibble::tibble(
      factor   = character(),
      a        = double(), se_a   = double(), p_a   = double(),
      b        = double(), se_b   = double(), p_b   = double(),
      c_direct = double(), se_c   = double(), p_c   = double(),
      indirect = double(), se_ind = double(), p_ind = double(),
      total    = double(), se_tot = double(), p_tot = double()
    ))
  }

  dplyr::bind_rows(out)
}

build_sem_matrix_data <- function(effects_main,
                                  effects_int = NULL,
                                  resp_var,
                                  species,
                                  include_interaction,
                                  swc_source = "measured",
                                  p_sig = 0.05,
                                  phase = NA_character_) {
  required_cols <- c("factor", "a", "p_a", "b", "p_b", "c_direct", "p_c", "indirect", "p_ind", "total", "p_tot")
  missing_cols <- setdiff(required_cols, names(effects_main))
  if (length(missing_cols)) {
    stop("effects_main is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  make_panel <- function(df, effect_class = "main") {
    p_swc <- df %>%
      dplyr::transmute(
        treatment = factor,
        response_var = "swc",
        path_type = "treatment_to_swc",
        estimate = a,
        p_value = p_a,
        effect_class = effect_class
      )

    p_dir <- df %>%
      dplyr::transmute(
        treatment = factor,
        response_var = resp_var,
        path_type = "direct",
        estimate = c_direct,
        p_value = p_c,
        effect_class = effect_class
      )

    p_ind <- df %>%
      dplyr::transmute(
        treatment = factor,
        response_var = resp_var,
        path_type = "indirect",
        estimate = indirect,
        p_value = p_ind,
        effect_class = effect_class
      )

    p_tot <- df %>%
      dplyr::transmute(
        treatment = factor,
        response_var = resp_var,
        path_type = "total",
        estimate = total,
        p_value = p_tot,
        effect_class = effect_class
      )

    p_swc_target <- tibble::tibble(
      treatment = "swc",
      response_var = resp_var,
      path_type = "swc_to_target",
      estimate = df$b[[1]],
      p_value = df$p_b[[1]],
      effect_class = effect_class
    )

    dplyr::bind_rows(p_swc, p_dir, p_ind, p_tot, p_swc_target)
  }

  out <- make_panel(effects_main, effect_class = "main")

  if (!is.null(effects_int) && nrow(effects_int) > 0L) {
    out <- dplyr::bind_rows(out, make_panel(effects_int, effect_class = "interaction"))
  }

  out %>%
    dplyr::mutate(
      significant = !is.na(p_value) & p_value < p_sig,
      estimate_sig = dplyr::if_else(significant, estimate, NA_real_),
      species = species,
      include_interaction = include_interaction,
      swc_source = swc_source,
      phase = phase
    ) %>%
    alinv_apply_response_orientation(
      resp_col = "response_var",
      estimate_col = "estimate",
      estimate_sig_col = "estimate_sig"
    )
}

add_sem_phase <- function(df, phase_sel = "all") {
  df2 <- df %>%
    dplyr::mutate(
      phase_window = dplyr::case_when(
        lubridate::month(as.Date(date)) <= 6 ~ "until June",
        lubridate::month(as.Date(date)) <= 8 ~ "July-August",
        TRUE ~ "September-end"
      )
    )

  if (!identical(phase_sel, "all")) {
    df2 <- df2 %>% dplyr::filter(.data$phase_window == phase_sel)
  }

  df2
}

summarize_sem_matrix_significant_mean <- function(matrix_df) {
  needed <- c("species", "include_interaction", "swc_source", "phase", "response_var", "treatment", "path_type", "estimate", "significant")
  missing <- setdiff(needed, names(matrix_df))
  if (length(missing)) {
    stop("matrix_df is missing required columns: ", paste(missing, collapse = ", "))
  }

  matrix_df %>%
    dplyr::filter(significant) %>%
    dplyr::group_by(
      species, include_interaction, swc_source, phase,
      response_var, treatment, path_type
    ) %>%
    dplyr::summarise(
      mean_estimate = mean(estimate, na.rm = TRUE),
      n_significant = dplyr::n(),
      .groups = "drop"
    )
}

plot_overarching_sem_effect_matrices <- function(summary_df,
                                                 target_path_type = c("direct", "indirect", "total"),
                                                 title = NULL) {
  target_path_type <- match.arg(target_path_type)

  df_target <- summary_df %>%
    dplyr::filter(path_type == target_path_type, treatment != "swc") %>%
    dplyr::mutate(panel = paste0("Treatment -> target (", target_path_type, ")"))

  df_swc <- summary_df %>%
    dplyr::filter(path_type == "treatment_to_swc", treatment != "swc") %>%
    dplyr::mutate(panel = "Treatment -> SWC")

  df_swc_target <- summary_df %>%
    dplyr::filter(path_type == "swc_to_target") %>%
    dplyr::mutate(panel = "SWC -> target")

  df_plot <- dplyr::bind_rows(df_target, df_swc, df_swc_target)
  if (!nrow(df_plot)) {
    stop("No rows available for overarching SEM matrix plotting.")
  }

  max_abs <- max(abs(df_plot$mean_estimate), na.rm = TRUE)
  if (!is.finite(max_abs) || max_abs == 0) {
    max_abs <- 1
  }

  df_plot <- df_plot %>%
    dplyr::mutate(
      panel = factor(panel, levels = c(
        paste0("Treatment -> target (", target_path_type, ")"),
        "Treatment -> SWC",
        "SWC -> target"
      )),
      alpha_val = pmin(1, abs(mean_estimate) / max_abs)
    )

  ggplot(df_plot, aes(x = treatment, y = response_var)) +
    geom_tile(aes(fill = mean_estimate, alpha = alpha_val), color = "grey90", linewidth = 0.2) +
    scale_fill_gradient2(
      low = "indianred3",
      mid = "white",
      high = "steelblue4",
      midpoint = 0,
      na.value = "white",
      name = "Mean oriented std. effect"
    ) +
    scale_alpha(range = c(0, 1), guide = "none") +
    facet_wrap(~panel, ncol = 3, scales = "free_x") +
    labs(x = NULL, y = NULL, title = title) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.grid = element_blank())
}

plot_sem_effect_matrices <- function(matrix_df,
                                     target_path_type = c("direct", "indirect", "total"),
                                     title = NULL) {
  target_path_type <- match.arg(target_path_type)

  df_target <- matrix_df %>%
    dplyr::filter(path_type == target_path_type, treatment != "swc") %>%
    dplyr::mutate(panel = paste0("Treatment -> target (", target_path_type, ")"))

  df_swc <- matrix_df %>%
    dplyr::filter(path_type == "treatment_to_swc", treatment != "swc") %>%
    dplyr::mutate(panel = "Treatment -> SWC")

  df_swc_target <- matrix_df %>%
    dplyr::filter(path_type == "swc_to_target") %>%
    dplyr::mutate(panel = "SWC -> target")

  df_plot <- dplyr::bind_rows(df_target, df_swc, df_swc_target) %>%
    dplyr::mutate(
      panel = factor(panel, levels = c(
        paste0("Treatment -> target (", target_path_type, ")"),
        "Treatment -> SWC",
        "SWC -> target"
      )),
      alpha_val = dplyr::if_else(
        is.na(estimate_sig),
        0,
        pmin(1, abs(estimate_sig) / max(abs(estimate_sig), na.rm = TRUE))
      )
    )

  if (!nrow(df_plot)) {
    stop("No SEM matrix rows available for plotting.")
  }

  ggplot(df_plot, aes(x = treatment, y = response_var)) +
    geom_tile(aes(fill = estimate_sig, alpha = alpha_val), color = "grey90", linewidth = 0.2) +
    scale_fill_gradient2(
      low = "indianred3",
      mid = "white",
      high = "steelblue4",
      midpoint = 0,
      na.value = "white",
      name = "Oriented std. effect"
    ) +
    scale_alpha(range = c(0, 1), guide = "none") +
    facet_wrap(~panel, ncol = 3, scales = "free_x") +
    labs(
      x = NULL,
      y = NULL,
      title = title
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1),
      panel.grid = element_blank()
    )
}

sem_heatmap_specs <- function() {
  tibble::tribble(
    ~path_type, ~file_stub, ~panel_title,
    "direct", "direct-effect_of_treatment_on_target", "Direct effect of treatment on target",
    "indirect", "indirect-effect_of_treatment_on_target", "Indirect effect of treatment on target",
    "total", "total-effect_of_treatment_on_target", "Total effect of treatment on target",
    "treatment_to_swc", "effect_of_treatment_on_swc", "Effect of treatment on SWC",
    "swc_to_target", "effect_of_swc_on_target", "Effect of SWC on target"
  )
}

sem_heatmap_treatment_labels <- function() {
  c(
    robinia = "Robinia: without -> with",
    precipitation = "Precipitation: control -> drought",
    culture = "Culture: mono -> mixed",
    soiltype = "Soil: wetter soil (robinia soil) -> drier soil (beech soil)",
    extreme_event = "Extreme event: no -> yes",
    swc = "SWC"
  )
}

compute_shared_sem_heatmap_limit <- function(matrix_df) {
  vals <- matrix_df$estimate_sig
  vals <- vals[is.finite(vals)]

  if (!length(vals)) {
    vals <- matrix_df$estimate
    vals <- vals[is.finite(vals)]
  }

  if (!length(vals)) {
    return(1)
  }

  max_abs <- max(abs(vals), na.rm = TRUE)
  if (!is.finite(max_abs) || max_abs <= 0) {
    1
  } else {
    max_abs
  }
}

prepare_sem_heatmap_panel <- function(matrix_df,
                                      path_type,
                                      resp_labels = NULL,
                                      treat_labels = sem_heatmap_treatment_labels()) {
  resp_labels <- resp_labels %||% c()
  treat_labels <- treat_labels %||% sem_heatmap_treatment_labels()

  df_panel <- matrix_df %>%
    dplyr::filter(.data$path_type == .env$path_type)

  if (!nrow(df_panel)) {
    return(tibble::tibble())
  }

  if (!"estimate_raw" %in% names(df_panel)) {
    df_panel$estimate_raw <- df_panel$estimate
  }
  if (!"estimate_sig_raw" %in% names(df_panel)) {
    df_panel$estimate_sig_raw <- df_panel$estimate_sig
  }

  if (identical(path_type, "treatment_to_swc")) {
    df_panel <- df_panel %>%
      dplyr::filter(.data$treatment != "swc") %>%
      dplyr::mutate(
        row_label = dplyr::recode(.data$treatment, !!!treat_labels),
        col_label = "SWC"
      )
    row_order <- unname(treat_labels[names(treat_labels) %in% unique(df_panel$treatment)])
    row_order <- row_order[!is.na(row_order)]
    col_order <- "SWC"
  } else if (identical(path_type, "swc_to_target")) {
    df_panel <- df_panel %>%
      dplyr::mutate(
        row_label = "SWC",
        col_label = dplyr::recode(.data$response_var, !!!resp_labels, .default = .data$response_var)
      )
    row_order <- "SWC"
    resp_levels <- unique(df_panel$response_var)
    col_order <- unname(c(resp_labels[resp_levels], stats::setNames(resp_levels, resp_levels)))
    col_order <- col_order[!is.na(col_order) & !duplicated(col_order)]
  } else {
    df_panel <- df_panel %>%
      dplyr::filter(.data$treatment != "swc") %>%
      dplyr::mutate(
        row_label = dplyr::recode(.data$treatment, !!!treat_labels),
        col_label = dplyr::recode(.data$response_var, !!!resp_labels, .default = .data$response_var)
      )
    row_order <- unname(treat_labels[names(treat_labels) %in% unique(df_panel$treatment)])
    row_order <- row_order[!is.na(row_order)]
    resp_levels <- unique(df_panel$response_var)
    col_order <- unname(c(resp_labels[resp_levels], stats::setNames(resp_levels, resp_levels)))
    col_order <- col_order[!is.na(col_order) & !duplicated(col_order)]
  }

  df_panel %>%
    dplyr::transmute(
      row_label = .data$row_label,
      col_label = .data$col_label,
      value = .data$estimate_sig,
      value_raw = dplyr::coalesce(.data$estimate_sig_raw, .data$estimate_raw, .data$estimate)
    ) %>%
    dplyr::mutate(
      row_label = factor(.data$row_label, levels = row_order),
      col_label = factor(.data$col_label, levels = col_order)
    ) %>%
    dplyr::arrange(.data$row_label, .data$col_label)
}

sem_heatmap_panel_to_wide <- function(panel_df) {
  if (!nrow(panel_df)) {
    return(tibble::tibble(row_label = character()))
  }

  row_order <- levels(panel_df$row_label)
  col_order <- levels(panel_df$col_label)

  tidyr::expand_grid(
    row_label = factor(row_order, levels = row_order),
    col_label = factor(col_order, levels = col_order)
  ) %>%
    dplyr::left_join(panel_df, by = c("row_label", "col_label")) %>%
    dplyr::select(.data$row_label, .data$col_label, .data$value) %>%
    dplyr::mutate(
      row_label = as.character(.data$row_label),
      col_label = as.character(.data$col_label)
    ) %>%
    tidyr::pivot_wider(names_from = .data$col_label, values_from = .data$value) %>%
    dplyr::rename(row = .data$row_label)
}

write_sem_heatmap_csvs <- function(matrix_df,
                                   species,
                                   export_dir,
                                   resp_labels = NULL,
                                   treat_labels = sem_heatmap_treatment_labels()) {
  dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)
  specs <- sem_heatmap_specs()

  if (is.null(matrix_df) || !"species" %in% names(matrix_df) || !nrow(matrix_df)) {
    purrr::walk(seq_len(nrow(specs)), function(i) {
      spec <- specs[i, ]
      out_file <- file.path(export_dir, paste0(species, "-", spec$file_stub, ".csv"))
      readr::write_csv(tibble::tibble(row = character()), out_file)
    })
    return(invisible(NULL))
  }

  purrr::walk(seq_len(nrow(specs)), function(i) {
    spec <- specs[i, ]
    panel_df <- prepare_sem_heatmap_panel(
      matrix_df = matrix_df %>% dplyr::filter(.data$species == .env$species),
      path_type = spec$path_type,
      resp_labels = resp_labels,
      treat_labels = treat_labels
    )
    wide_df <- sem_heatmap_panel_to_wide(panel_df)
    out_file <- file.path(export_dir, paste0(species, "-", spec$file_stub, ".csv"))
    readr::write_csv(wide_df, out_file)
  })
}

plot_sem_heatmap_panel <- function(panel_df,
                                   limit,
                                   title = NULL,
                                   annotate_values = FALSE,
                                   value_digits = 2,
                                   value_text_size = 3.1) {
  if (!nrow(panel_df)) {
    return(
      ggplot() +
        theme_void() +
        ggtitle(title %||% "No SEM heatmap data")
    )
  }

  label_df <- if (isTRUE(annotate_values)) {
    panel_df %>%
      dplyr::filter(is.finite(.data$value)) %>%
      dplyr::mutate(
        value_label = formatC(.data$value, format = "f", digits = value_digits),
        value_color = dplyr::if_else(abs(.data$value) >= 0.55 * limit, "white", "black")
      )
  } else {
    tibble::tibble()
  }

  p <- ggplot(
    panel_df,
    aes(x = .data$col_label, y = .data$row_label, fill = .data$value)
  ) +
    geom_tile(color = "grey90", linewidth = 0.2) +
    scale_fill_gradient2(
      low = "indianred3",
      mid = "white",
      high = "steelblue4",
      midpoint = 0,
      limits = c(-limit, limit),
      na.value = "white",
      name = "Oriented std. effect"
    ) +
    scale_color_identity() +
    labs(
      x = NULL,
      y = NULL,
      title = title
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1),
      panel.grid = element_blank()
    )

  if (nrow(label_df)) {
    p <- p +
      geom_text(
        data = label_df,
        mapping = aes(
          x = .data$col_label,
          y = .data$row_label,
          label = .data$value_label,
          color = .data$value_color
        ),
        inherit.aes = FALSE,
        size = value_text_size,
        fontface = "bold"
      )
  }

  p
}


build_sem_graph_components <- function(effects_main,
                                       effects_int = NULL,
                                       sem_mod,
                                       mod_swc = NULL,
                                       mod_resp = NULL,
                                       resp_var,
                                       species,
                                       soil_type,
                                       include_interaction = TRUE,
                                       modeled_factors = NULL,
                                       p_sig = 0.1) {
  display_sign <- alinv_response_display_sign(resp_var)[[1]]

  if (!is.null(sem_mod) && requireNamespace("piecewiseSEM", quietly = TRUE)) {
    r2_tbl <- piecewiseSEM::rsquared(sem_mod)
    r2_swc_marg <- r2_tbl$Marginal[r2_tbl$Response == "swc"]
    r2_swc_cond <- r2_tbl$Conditional[r2_tbl$Response == "swc"]
    r2_resp_marg <- r2_tbl$Marginal[r2_tbl$Response == "y"]
    r2_resp_cond <- r2_tbl$Conditional[r2_tbl$Response == "y"]

    r2_label <- sprintf(
      "SWC:    R²(marg) = %.2f, R²(cond) = %.2f\n%s: R²(marg) = %.2f, R²(cond) = %.2f",
      r2_swc_marg, r2_swc_cond,
      resp_var, r2_resp_marg, r2_resp_cond
    )
  } else {
    r2_swc <- alinv_lmm_r2(mod_swc)
    r2_resp <- alinv_lmm_r2(mod_resp)
    r2_swc_marg <- r2_swc$r2_marginal[[1]]
    r2_swc_cond <- r2_swc$r2_conditional[[1]]
    r2_resp_marg <- r2_resp$r2_marginal[[1]]
    r2_resp_cond <- r2_resp$r2_conditional[[1]]
    r2_label <- if (all(is.na(c(r2_swc_marg, r2_swc_cond, r2_resp_marg, r2_resp_cond)))) {
      "R² unavailable in this R session."
    } else {
      sprintf(
        "SWC:    R²(marg) = %.2f, R²(cond) = %.2f\n%s: R²(marg) = %.2f, R²(cond) = %.2f",
        r2_swc_marg, r2_swc_cond,
        resp_var, r2_resp_marg, r2_resp_cond
      )
    }
  }

  resp_node <- resp_var

  # automatic node layout
  # fixed order of main predictors on the left
  modeled_factors <- modeled_factors %||% c(
    "extreme_event", "robinia", "precipitation", "soiltype", "culture"
  )
  main_nodes_order <- modeled_factors

  nodes_main <- tibble::tibble(
    node = main_nodes_order,
    x    = 0.5,
    y    = seq(from = length(main_nodes_order), to = 1, by = -1)
  )

  # interaction nodes in a second column, slightly to the right
  if (!is.null(effects_int) && nrow(effects_int) > 0L && isTRUE(include_interaction)) {
    int_labels <- unique(effects_int$factor)

    nodes_int <- lapply(int_labels, function(lbl) {
      parts <- strsplit(lbl, ":", fixed = TRUE)[[1]]
      fac1 <- parts[1]
      fac2 <- parts[2]

      y1 <- nodes_main$y[nodes_main$node == fac1]
      y2 <- nodes_main$y[nodes_main$node == fac2]

      # place interaction between its two parents; fall back to center if needed
      y_int <- if (length(y1) && length(y2)) {
        mean(c(y1, y2))
      } else if (length(y1)) {
        y1
      } else if (length(y2)) {
        y2
      } else {
        mean(nodes_main$y)
      }

      tibble::tibble(
        node = lbl,
        x    = 7,
        y    = y_int
      )
    }) %>%
      dplyr::bind_rows()
  } else {
    nodes_int <- tibble::tibble(
      node = character(), x = double(), y = double()
    )
  }

  # SWC and response nodes further to the right
  # y_swc  <- max(nodes_main$y) + 0.25
  # y_resp <- mean(range(nodes_main$y))

  y_swc <- mean(range(nodes_main$y)) * 1.15
  y_resp <- mean(range(nodes_main$y)) * 0.85

  node_swc <- tibble::tibble(
    node = "swc",
    x    = 4,
    y    = y_swc
  )

  node_resp <- tibble::tibble(
    node = resp_node,
    x    = 4,
    y    = y_resp
  )

  nodes <- dplyr::bind_rows(nodes_main, nodes_int, node_swc, node_resp)

  edges_a <- effects_main %>%
    dplyr::transmute(
      from = factor,
      to   = "swc",
      est_raw = a,
      est  = a,
      se   = se_a,
      p    = p_a,
      path = "a (via swc)"
    )

  edges_c <- effects_main %>%
    dplyr::transmute(
      from = factor,
      to   = resp_node,
      est_raw = c_direct,
      est  = c_direct * display_sign,
      se   = se_c,
      p    = p_c,
      path = "c' (direct)"
    )

  b_row <- effects_main %>% dplyr::slice(1)
  edge_b <- b_row %>%
    dplyr::transmute(
      from = "swc",
      to   = resp_node,
      est_raw = b,
      est  = b * display_sign,
      se   = se_b,
      p    = p_b,
      path = "b (swc → response)"
    )

  edges_int_a <- NULL
  edges_int_c <- NULL

  if (!is.null(effects_int) && nrow(effects_int) > 0L && isTRUE(include_interaction)) {
    edges_int_a <- effects_int %>%
      dplyr::transmute(
        from = factor,
        to   = "swc",
        est_raw = a,
        est  = a,
        se   = se_a,
        p    = p_a,
        path = "interaction via swc"
      )

    edges_int_c <- effects_int %>%
      dplyr::transmute(
        from = factor,
        to   = resp_node,
        est_raw = c_direct,
        est  = c_direct * display_sign,
        se   = se_c,
        p    = p_c,
        path = "interaction direct"
      )
  }

  edges <- dplyr::bind_rows(edges_a, edges_c, edge_b, edges_int_a, edges_int_c) %>%
    dplyr::mutate(
      alpha = dplyr::case_when(
        # p < 0.01 ~ 1.00,
        # p < 0.05 ~ 0.5,
        # p < 0.10 ~ 0,
        # TRUE ~ 0
        p < 0.01 ~ 1.00,
        p < 0.05 ~ 0.75,
        p < 0.10 ~ 0.50,
        TRUE ~ 0.25
      ),
      col_dir = dplyr::if_else(est > 0, "positive", "negative"),
      weight = abs(est)
    )

  # Keep only significant edges
  edges <- edges %>%
    dplyr::filter(!is.na(p) & p < p_sig)

  # Identify nodes that appear in at least one significant edge
  used_nodes <- unique(c(edges$from, edges$to))

  # Main nodes we always keep (even if they have no sig edges)
  main_nodes_order <- c(modeled_factors, "swc", resp_node)

  # Keep:
  # - all main nodes
  # - only those interaction nodes that have at least one sig edge
  nodes <- nodes %>%
    dplyr::filter(node %in% main_nodes_order | node %in% used_nodes)

  edges_plot <- edges %>%
    dplyr::left_join(nodes, by = c("from" = "node")) %>%
    dplyr::rename(x_from = x, y_from = y) %>%
    dplyr::left_join(nodes, by = c("to" = "node")) %>%
    dplyr::rename(x_to = x, y_to = y) %>%
    dplyr::mutate(
      x_end = x_to,
      y_end = y_to,
      lw    = scales::rescale(weight, to = c(0.3, 2.5))
    )

  title_txt <- paste0(
    "SEM: ", resp_var,
    " (species: ", species,
    ", soil: ", soil_type,
    if (include_interaction) ", with robinia:extreme" else ", no robinia:extreme",
    ")"
  )

  metrics <- tibble::tibble(
    resp_var = resp_var,
    species = species,
    soil_filter = soil_type,
    include_interaction = isTRUE(include_interaction),
    p_sig = p_sig,
    r2_swc_marg = r2_swc_marg,
    r2_swc_cond = r2_swc_cond,
    r2_resp_marg = r2_resp_marg,
    r2_resp_cond = r2_resp_cond,
    r2_label = r2_label,
    title = title_txt
  )

  list(
    nodes = nodes,
    edges = edges_plot,
    metrics = metrics
  )
}

plot_sem_graph_components <- function(nodes_df,
                                      edges_df,
                                      metrics_df = NULL) {
  metric_value <- function(col, default = NA) {
    if (is.null(metrics_df) || !nrow(metrics_df) || !col %in% names(metrics_df)) {
      return(default)
    }
    metrics_df[[col]][[1]]
  }

  title_txt <- metric_value("title", "SEM graph")
  r2_label <- metric_value("r2_label", "")

  if (is.null(nodes_df) || !nrow(nodes_df)) {
    return(ggplot() + theme_void() + ggtitle(title_txt))
  }

  max_x <- max(nodes_df$x) + 0.5
  max_y <- max(nodes_df$y) + 0.7
  min_y <- min(nodes_df$y) - 0.7

  ggplot() +
    geom_segment(
      data = edges_df,
      aes(
        x = x_from, y = y_from,
        xend = x_end, yend = y_end,
        size = lw,
        colour = col_dir,
        alpha = alpha
      ),
      arrow = grid::arrow(length = grid::unit(0.5, "cm"), type = "open"),
      lineend = "round"
    ) +
    geom_label(
      data = nodes_df,
      aes(x = x, y = y, label = node),
      size = 3.5,
      label.size = 0.3,
      label.r = grid::unit(0.15, "lines"),
      fill = "white",
      label.padding = grid::unit(0.15, "lines")
    ) +
    annotate(
      "text",
      x = max_x,
      y = max_y,
      label = r2_label,
      hjust = 1,
      size = 3
    ) +
    scale_size_identity() +
    scale_alpha_identity() +
    scale_colour_manual(
      values = c(positive = "steelblue4", negative = "indianred3"),
      name   = "Effect direction"
    ) +
    coord_equal(
      xlim   = c(0, max_x + 0.2),
      ylim   = c(min_y, max_y + 0.2),
      expand = FALSE
    ) +
    labs(title = title_txt) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position      = c(0.5, 0.06),
      legend.justification = c(0.5, 0),
      axis.text            = element_blank(),
      axis.title           = element_blank(),
      panel.grid           = element_blank(),
      plot.margin          = margin(5, 15, 5, 15)
    )
}

plot_sem_graph <- function(effects_main,
                           effects_int = NULL,
                           sem_mod,
                           mod_swc = NULL,
                           mod_resp = NULL,
                           resp_var,
                           species,
                           soil_type,
                           include_interaction = TRUE,
                           modeled_factors = NULL) {
  graph_bits <- build_sem_graph_components(
    effects_main = effects_main,
    effects_int = effects_int,
    sem_mod = sem_mod,
    mod_swc = mod_swc,
    mod_resp = mod_resp,
    resp_var = resp_var,
    species = species,
    soil_type = soil_type,
    include_interaction = include_interaction,
    modeled_factors = modeled_factors
  )

  plot_sem_graph_components(
    nodes_df = graph_bits$nodes,
    edges_df = graph_bits$edges,
    metrics_df = graph_bits$metrics
  )
}

sem_cache_path <- function(type = "tree",
                           data_name,
                           resp_var,
                           species,
                           soil_type = "both",
                           include_soil_treatment = NULL,
                           phase_window = "all",
                           include_interaction = TRUE,
                           scale_all_numeric = TRUE,
                           do_rfe = FALSE,
                           aic_improve = 2,
                           swc_source = "measured") {
  include_soil_treatment <- alinv_resolve_include_soil_treatment(
    include_soil_treatment = include_soil_treatment,
    soil_filter = soil_type
  )

  out_dir <- alinv_data_path("model-sem", create_dir = TRUE)
  resp_tag <- if (!is.null(resp_var)) resp_var else "default"
  int_tag <- if (isTRUE(include_interaction)) "int" else "noInt"
  scale_tag <- if (isTRUE(scale_all_numeric)) "scaled" else "unscaled"
  rfe_tag <- if (isTRUE(do_rfe)) paste0("rfeAIC", aic_improve) else "noRFE"
  swc_tag <- if (swc_source == "measured") "swcMeas" else "swcImputed"
  phase_tag <- gsub("[^a-zA-Z0-9]+", "", tolower(phase_window))
  soil_mode_tag <- alinv_soil_mode_tag(
    soil_filter = soil_type,
    include_soil_treatment = include_soil_treatment
  )

  file_name <- paste0(
    "sem-",
    type, "-",
    data_name, "-",
    resp_tag, "-",
    species, "-",
    soil_mode_tag, "-",
    int_tag, "-",
    scale_tag, "-",
    rfe_tag, "-",
    phase_tag, "-",
    swc_tag,
    ".rds"
  )

  file.path(out_dir, file_name)
}

augment_sem_result_for_exports <- function(result,
                                           resp_var,
                                           species,
                                           soil_type,
                                           include_interaction = TRUE,
                                           modeled_factors = NULL) {
  if (is.null(result)) {
    return(result)
  }

  modeled_factors <- modeled_factors %||% result$modeled_factors
  if (is.null(modeled_factors) && !is.null(result$effects) && nrow(result$effects)) {
    modeled_factors <- unique(as.character(result$effects$factor))
  }

  if (is.null(result$effects) || !nrow(result$effects)) {
    result$modeled_factors <- modeled_factors
    result$graph_nodes <- tibble::tibble()
    result$graph_edges <- tibble::tibble()
    result$graph_metrics <- tibble::tibble(
      resp_var = resp_var,
      species = species,
      soil_filter = soil_type,
      include_interaction = isTRUE(include_interaction),
      p_sig = 0.1,
      r2_swc_marg = NA_real_,
      r2_swc_cond = NA_real_,
      r2_resp_marg = NA_real_,
      r2_resp_cond = NA_real_,
      r2_label = "",
      title = paste0("SEM: ", resp_var, " (species: ", species, ", soil: ", soil_type, ")")
    )
    return(result)
  }

  display_sign <- alinv_response_display_sign(resp_var)[[1]]
  orient_effect_table <- function(df_effects) {
    if (is.null(df_effects) || !nrow(df_effects)) {
      return(df_effects)
    }

    df_effects %>%
      dplyr::mutate(
        display_sign = display_sign,
        a_display = .data$a,
        b_display = .data$b * .data$display_sign,
        c_direct_display = .data$c_direct * .data$display_sign,
        indirect_display = .data$indirect * .data$display_sign,
        total_display = .data$total * .data$display_sign
      )
  }

  graph_bits <- build_sem_graph_components(
    effects_main = result$effects %||% tibble::tibble(),
    effects_int = result$effects_int %||% tibble::tibble(),
    sem_mod = result$sem,
    mod_swc = result$mod_swc %||% NULL,
    mod_resp = result$mod_resp %||% NULL,
    resp_var = resp_var,
    species = species,
    soil_type = soil_type,
    include_interaction = include_interaction,
    modeled_factors = modeled_factors
  )

  result$effects <- orient_effect_table(result$effects)
  result$effects_int <- orient_effect_table(result$effects_int)
  if (!is.null(result$matrix_data) && nrow(result$matrix_data)) {
    if (!"estimate_raw" %in% names(result$matrix_data)) {
      result$matrix_data <- result$matrix_data %>%
        alinv_apply_response_orientation(
          resp_col = "response_var",
          estimate_col = "estimate",
          estimate_sig_col = "estimate_sig"
        )
    }
    result$matrix_overarching <- summarize_sem_matrix_significant_mean(result$matrix_data)
  }
  result$modeled_factors <- modeled_factors
  result$graph_nodes <- graph_bits$nodes
  result$graph_edges <- graph_bits$edges
  result$graph_metrics <- graph_bits$metrics
  result$plot <- plot_sem_graph_components(
    nodes_df = graph_bits$nodes,
    edges_df = graph_bits$edges,
    metrics_df = graph_bits$metrics
  )
  result
}

extract_sem_model_performance <- function(result,
                                          data_name,
                                          resp_var,
                                          species,
                                          soil_type,
                                          include_soil_treatment = NULL,
                                          swc_source = "measured") {
  include_soil_treatment <- alinv_resolve_include_soil_treatment(
    include_soil_treatment = include_soil_treatment,
    soil_filter = soil_type
  )

  metric_value <- function(col, default = NA_real_) {
    metrics_df <- result$graph_metrics %||% tibble::tibble()
    if (!nrow(metrics_df) || !col %in% names(metrics_df)) {
      return(default)
    }
    metrics_df[[col]][[1]]
  }

  data_df <- result$data %||% tibble::tibble()
  mod_swc <- result$mod_swc %||% NULL
  mod_resp <- result$mod_resp %||% NULL
  r2_swc <- alinv_lmm_r2(mod_swc)
  r2_resp <- alinv_lmm_r2(mod_resp)

  safe_stat <- function(expr) {
    tryCatch(expr, error = function(e) NA_real_)
  }

  tibble::tibble(
    model_scope = "SEM",
    data_name = data_name,
    resp_var = resp_var,
    species = species,
    soil_filter = soil_type,
    include_soil_treatment = isTRUE(include_soil_treatment),
    swc_source = swc_source,
    modeled_factors = paste(result$modeled_factors %||% character(), collapse = ", "),
    n_obs = nrow(data_df),
    n_boxes = if ("boxlabel" %in% names(data_df)) dplyr::n_distinct(data_df$boxlabel) else NA_integer_,
    n_trees = if ("tree_id" %in% names(data_df)) dplyr::n_distinct(data_df$tree_id) else NA_integer_,
    aic_swc = safe_stat(stats::AIC(mod_swc)),
    bic_swc = safe_stat(stats::BIC(mod_swc)),
    aic_response = safe_stat(stats::AIC(mod_resp)),
    bic_response = safe_stat(stats::BIC(mod_resp)),
    r2_swc_marginal = dplyr::coalesce(metric_value("r2_swc_marg"), r2_swc$r2_marginal[[1]]),
    r2_swc_conditional = dplyr::coalesce(metric_value("r2_swc_cond"), r2_swc$r2_conditional[[1]]),
    r2_response_marginal = dplyr::coalesce(metric_value("r2_resp_marg"), r2_resp$r2_marginal[[1]]),
    r2_response_conditional = dplyr::coalesce(metric_value("r2_resp_cond"), r2_resp$r2_conditional[[1]]),
    piecewise_r2_available = !all(is.na(c(
      dplyr::coalesce(metric_value("r2_swc_marg"), r2_swc$r2_marginal[[1]]),
      dplyr::coalesce(metric_value("r2_swc_cond"), r2_swc$r2_conditional[[1]]),
      dplyr::coalesce(metric_value("r2_resp_marg"), r2_resp$r2_marginal[[1]]),
      dplyr::coalesce(metric_value("r2_resp_cond"), r2_resp$r2_conditional[[1]])
    ))),
    note = result$note %||% ""
  )
}

write_sem_csv_bundle <- function(result,
                                 export_stem) {
  result <- result %||% list()

  readr::write_csv(result$data %||% tibble::tibble(), paste0(export_stem, "-data.csv"))
  readr::write_csv(result$effects %||% tibble::tibble(), paste0(export_stem, "-effects_main.csv"))
  readr::write_csv(result$effects_int %||% tibble::tibble(), paste0(export_stem, "-effects_interactions.csv"))
  readr::write_csv(result$matrix_data %||% tibble::tibble(), paste0(export_stem, "-matrix_data.csv"))
  readr::write_csv(result$matrix_overarching %||% tibble::tibble(), paste0(export_stem, "-matrix_overarching.csv"))
  readr::write_csv(result$graph_nodes %||% tibble::tibble(), paste0(export_stem, "-graph_nodes.csv"))
  readr::write_csv(result$graph_edges %||% tibble::tibble(), paste0(export_stem, "-graph_edges.csv"))
  readr::write_csv(result$graph_metrics %||% tibble::tibble(), paste0(export_stem, "-graph_metrics.csv"))
}

run_sem_for_trait <- function(type = "tree",
                              data_name,
                              resp_var,
                              species,
                              soil_type = "both",
                              include_soil_treatment = NULL,
                              phase_window = "all",
                              include_interaction = TRUE,
                              scale_all_numeric = TRUE,
                              do_rfe = FALSE,
                              aic_improve = 2,
                              swc_source = "measured",
                              force_run = FALSE) {
  include_soil_treatment <- alinv_resolve_include_soil_treatment(
    include_soil_treatment = include_soil_treatment,
    soil_filter = soil_type
  )

  # --- 0) caching: save / load SEM run for today ---
  cache_path <- sem_cache_path(
    type = type,
    data_name = data_name,
    resp_var = resp_var,
    species = species,
    soil_type = soil_type,
    include_soil_treatment = include_soil_treatment,
    phase_window = phase_window,
    include_interaction = include_interaction,
    scale_all_numeric = scale_all_numeric,
    do_rfe = do_rfe,
    aic_improve = aic_improve,
    swc_source = swc_source
  )
  export_stem <- tools::file_path_sans_ext(cache_path)

  if (file.exists(cache_path) && !force_run) {
    message("Loading cached SEM results from: ", cache_path)
    cached_result <- readRDS(cache_path)
    cache_has_error <- !is.null(cached_result$error) && nzchar(cached_result$error)

    if (cache_has_error) {
      message("Cached SEM result contained an error. Re-running: ", cache_path)
    } else {
      cached_result <- augment_sem_result_for_exports(
        result = cached_result,
        resp_var = resp_var,
        species = species,
        soil_type = soil_type,
        include_interaction = include_interaction
      )
      write_sem_csv_bundle(cached_result, export_stem = export_stem)
      return(cached_result)
    }
  }

  # 1) Load analysis-ready data in your project’s canonical way
  df_pre <- prepare_df_generic(
    type         = type,
    data_name    = data_name,
    resp_var     = resp_var,
    species_keep = species,
    standardize_response = isTRUE(scale_all_numeric),
    add_covars   = FALSE, # SEM doesn’t use the extra temporal covars
    covars_fun   = NULL,
    soil_type    = soil_type,
    include_soil_treatment = include_soil_treatment,
    swc_source   = swc_source
  )

  if (!nrow(df_pre)) {
    result <- list(
      data = tibble::tibble(),
      effects = tibble::tibble(),
      effects_int = tibble::tibble(),
      matrix_data = tibble::tibble(),
      matrix_plots = list(),
      matrix_overarching = tibble::tibble(),
      matrix_overarching_plots = list(),
      modeled_factors = character(),
      note = "No SEM data available after filtering.",
      plot = alinv_empty_plot(
        title = "No SEM data",
        subtitle = paste0("No model rows available for ", species, " / ", data_name, " / ", resp_var)
      )
    )
    result <- augment_sem_result_for_exports(
      result = result,
      resp_var = resp_var,
      species = species,
      soil_type = soil_type,
      include_interaction = include_interaction
    )
    saveRDS(result, cache_path)
    write_sem_csv_bundle(result, export_stem = export_stem)
    return(result)
  }

  # optional exploratory SEM by seasonal phase
  if (!identical(phase_window, "all")) {
    df_pre <- df_pre %>%
      dplyr::mutate(
        phase_window = dplyr::case_when(
          lubridate::month(as.Date(date)) <= 6 ~ "until June",
          lubridate::month(as.Date(date)) <= 8 ~ "July-August",
          TRUE ~ "September-end"
        )
      ) %>%
      dplyr::filter(.data$phase_window == phase_window)
  }

  # 2) Prepare specifically for SEM (generic y, scaling, etc.)
  df_sem_ready <- prepare_sem_data(
    df_prepared       = df_pre,
    scale_all_numeric = scale_all_numeric
  )

  if (!nrow(df_sem_ready)) {
    result <- list(
      data = df_sem_ready,
      effects = tibble::tibble(),
      effects_int = tibble::tibble(),
      matrix_data = tibble::tibble(),
      matrix_plots = list(),
      matrix_overarching = tibble::tibble(),
      matrix_overarching_plots = list(),
      modeled_factors = character(),
      note = "No SEM rows remained after preparation.",
      plot = alinv_empty_plot(
        title = "No SEM data",
        subtitle = paste0("No analyzable rows remained for ", species, " / ", data_name, " / ", resp_var)
      )
    )
    result <- augment_sem_result_for_exports(
      result = result,
      resp_var = resp_var,
      species = species,
      soil_type = soil_type,
      include_interaction = include_interaction
    )
    saveRDS(result, cache_path)
    write_sem_csv_bundle(result, export_stem = export_stem)
    return(result)
  }

  # 3) Fit SEM submodels
  mods <- fit_sem_models(
    df_sem_ready        = df_sem_ready,
    include_interaction = include_interaction,
    include_soil_treatment = include_soil_treatment,
    do_rfe              = do_rfe,
    aic_improve         = aic_improve
  )

  mod_swc <- mods$mod_swc
  mod_resp <- mods$mod_resp

  # 4) Combine into piecewise SEM
  sem_note <- NULL
  sem_mod <- if (requireNamespace("piecewiseSEM", quietly = TRUE)) {
    piecewiseSEM::psem(
      mod_swc,
      mod_resp,
      data = df_sem_ready
    )
  } else {
    sem_note <- "Package 'piecewiseSEM' unavailable; SEM path coefficients were computed from the fitted submodels and submodel R² values were estimated directly from the fitted mixed models."
    message(sem_note)
    NULL
  }

  # 5) Extract effects
  effects_main <- extract_sem_effects_with_se(
    mod_swc,
    mod_resp,
    factors = mods$used_factors
  )
  effects_int <- if (include_interaction) {
    extract_interaction_effect(
      mod_swc,
      mod_resp,
      treat_factors = mods$used_factors
    )
  } else {
    NULL
  }

  matrix_data <- build_sem_matrix_data(
    effects_main = effects_main,
    effects_int = effects_int,
    resp_var = resp_var,
    species = species,
    include_interaction = include_interaction,
    swc_source = swc_source,
    phase = phase_window
  )

  matrix_overarching <- summarize_sem_matrix_significant_mean(matrix_data)

  matrix_overarching_plots <- if (nrow(matrix_overarching)) {
    list(
      direct = plot_overarching_sem_effect_matrices(matrix_overarching, target_path_type = "direct", title = paste0("Overarching mean SEM matrix (direct): ", species, " - ", resp_var)),
      indirect = plot_overarching_sem_effect_matrices(matrix_overarching, target_path_type = "indirect", title = paste0("Overarching mean SEM matrix (indirect): ", species, " - ", resp_var)),
      total = plot_overarching_sem_effect_matrices(matrix_overarching, target_path_type = "total", title = paste0("Overarching mean SEM matrix (total): ", species, " - ", resp_var))
    )
  } else {
    list()
  }

  matrix_plots <- list(
    direct = plot_sem_effect_matrices(matrix_data, target_path_type = "direct", title = paste0("SEM matrix (direct): ", species, " - ", resp_var)),
    indirect = plot_sem_effect_matrices(matrix_data, target_path_type = "indirect", title = paste0("SEM matrix (indirect): ", species, " - ", resp_var)),
    total = plot_sem_effect_matrices(matrix_data, target_path_type = "total", title = paste0("SEM matrix (total): ", species, " - ", resp_var))
  )

  # 6) Plot graph with full metadata (species, soil_type, resp_var)
  p <- plot_sem_graph(
    effects_main = effects_main,
    effects_int = effects_int,
    sem_mod = sem_mod,
    mod_swc = mod_swc,
    mod_resp = mod_resp,
    resp_var = resp_var,
    species = species,
    soil_type = soil_type,
    include_interaction = include_interaction,
    modeled_factors = mods$used_factors
  )

  graph_bits <- build_sem_graph_components(
    effects_main = effects_main,
    effects_int = effects_int,
    sem_mod = sem_mod,
    mod_swc = mod_swc,
    mod_resp = mod_resp,
    resp_var = resp_var,
    species = species,
    soil_type = soil_type,
    include_interaction = include_interaction,
    modeled_factors = mods$used_factors
  )

  result <- list(
    data        = df_sem_ready,
    mod_swc     = mod_swc,
    mod_resp    = mod_resp,
    sem         = sem_mod,
    effects     = effects_main,
    effects_int = effects_int,
    matrix_data = matrix_data,
    matrix_plots = matrix_plots,
    matrix_overarching = matrix_overarching,
    matrix_overarching_plots = matrix_overarching_plots,
    modeled_factors = mods$used_factors,
    graph_nodes = graph_bits$nodes,
    graph_edges = graph_bits$edges,
    graph_metrics = graph_bits$metrics,
    phase_window = phase_window,
    phase_sem_exploratory = !identical(phase_window, "all"),
    n_rows_phase = nrow(df_sem_ready),
    include_soil_treatment = include_soil_treatment,
    note        = sem_note,
    plot        = p
  )

  saveRDS(result, cache_path)
  write_sem_csv_bundle(result, export_stem = export_stem)
  message("Saved SEM results to: ", cache_path)

  result
}
# ....................................................------------------------

build_sem_terms <- function(df,
                            include_interaction = TRUE,
                            include_soil_treatment = NULL) {
  include_soil_treatment <- alinv_resolve_include_soil_treatment(
    include_soil_treatment = include_soil_treatment,
    soil_filter = if ("soiltype" %in% names(df) &&
      dplyr::n_distinct(df$soiltype, na.rm = TRUE) > 1) "both" else "single"
  )

  # all potential treatment factors
  treat_factors_all <- c(
    "precipitation",
    "robinia",
    "culture",
    "extreme_event"
  )
  if (isTRUE(include_soil_treatment)) {
    treat_factors_all <- c("precipitation", "robinia", "soiltype", "culture", "extreme_event")
  }

  # keep only those that:
  #  - exist in df
  #  - have at least 2 distinct (non-NA) levels
  treat_factors <- treat_factors_all[
    treat_factors_all %in% names(df) &
      vapply(
        treat_factors_all,
        function(v) {
          if (!v %in% names(df)) {
            return(FALSE)
          }
          dplyr::n_distinct(df[[v]], na.rm = TRUE) > 1
        },
        logical(1)
      )
  ]

  # main effects used in SEM selection
  main_terms <- c("doy_c", "doy_c2", treat_factors)

  # all two-way interactions among treatments that actually vary
  if (isTRUE(include_interaction) && length(treat_factors) >= 2) {
    int_mat <- combn(treat_factors, 2)
    int_terms <- apply(int_mat, 2, function(x) paste(x, collapse = ":"))
  } else {
    int_terms <- character(0)
  }

  list(
    terms_swc = c(main_terms, int_terms), # fixed effects in swc model
    terms_resp = c(main_terms, int_terms), # fixed effects in response model (swc added later)
    used_factors = treat_factors # optional: to inspect which treatments survived
  )
}

backward_select_lmer_terms <- function(
  response,
  fixed_terms,
  random_part,
  data,
  always_keep = character(),
  aic_improve = 2
) {
  # ensure unique terms and remove always_keep from candidate pool
  fixed_terms <- unique(fixed_terms)
  curr_terms <- fixed_terms
  curr_terms <- union(curr_terms, always_keep) # just in case
  history <- list()

  # helper to fit a model and get AIC
  fit_and_aic <- function(terms) {
    rhs <- paste(terms, collapse = " + ")
    fml <- as.formula(
      paste(response, "~", rhs, "+", random_part)
    )
    mod <- lme4::lmer(fml, data = data, REML = FALSE)
    list(model = mod, aic = AIC(mod), formula = fml)
  }

  # initial full model
  best <- fit_and_aic(curr_terms)
  history[[length(history) + 1]] <- data.frame(
    step = 0L,
    dropped = NA_character_,
    k_terms = length(curr_terms),
    AIC = best$aic,
    stringsAsFactors = FALSE
  )

  repeat {
    # candidate droppable terms: not always_keep
    droppable <- setdiff(curr_terms, always_keep)
    if (!length(droppable)) break

    # evaluate dropping each term (respect hierarchy)
    candidates <- lapply(droppable, function(term_drop) {
      # if main effect: drop it and any interactions involving it
      if (!grepl(":", term_drop, fixed = TRUE)) {
        pattern <- paste0("(^|:)", term_drop, "(:|$)")
        to_drop <- grep(pattern, curr_terms, value = TRUE)
      } else {
        # interaction: drop only that interaction
        to_drop <- term_drop
      }
      new_terms <- setdiff(curr_terms, to_drop)
      # ensure always_keep are present
      new_terms <- union(new_terms, always_keep)
      if (!length(new_terms)) {
        return(NULL)
      }

      fit <- fit_and_aic(new_terms)
      data.frame(
        term_dropped = term_drop,
        aic = fit$aic,
        k_terms = length(new_terms),
        stringsAsFactors = FALSE
      )
    })

    candidates <- dplyr::bind_rows(candidates)
    # if nothing could be evaluated, stop
    if (nrow(candidates) == 0L) break

    # choose best candidate (lowest AIC)
    best_candidate <- candidates[which.min(candidates$aic), ]

    # check improvement
    delta_aic <- best_candidate$aic - best$aic
    if (delta_aic > -aic_improve) {
      # no substantial improvement
      break
    }

    # accept this drop
    term_drop <- best_candidate$term_dropped
    if (!grepl(":", term_drop, fixed = TRUE)) {
      pattern <- paste0("(^|:)", term_drop, "(:|$)")
      to_drop <- grep(pattern, curr_terms, value = TRUE)
    } else {
      to_drop <- term_drop
    }
    curr_terms <- setdiff(curr_terms, to_drop)
    curr_terms <- union(curr_terms, always_keep)

    # refit best model with updated terms
    best <- fit_and_aic(curr_terms)

    history[[length(history) + 1]] <- data.frame(
      step = length(history),
      dropped = paste(to_drop, collapse = ","),
      k_terms = length(curr_terms),
      AIC = best$aic,
      stringsAsFactors = FALSE
    )
  }

  history_df <- dplyr::bind_rows(history)
  list(
    selected_terms = curr_terms,
    history        = history_df,
    best_model     = best$model,
    best_formula   = best$formula,
    best_AIC       = best$aic
  )
}

fit_sem_models <- function(
  df_sem_ready,
  include_interaction = TRUE,
  include_soil_treatment = NULL,
  do_rfe = FALSE,
  aic_improve = 2
) {
  # build full set of main and interaction terms
  terms_list <- build_sem_terms(
    df_sem_ready,
    include_interaction = include_interaction,
    include_soil_treatment = include_soil_treatment
  )

  if (isTRUE(do_rfe)) {
    # SWC equation: swc ~ terms_swc + (1|boxlabel)
    sel_swc <- backward_select_lmer_terms(
      response    = "swc",
      fixed_terms = terms_list$terms_swc,
      random_part = "(1|boxlabel)",
      data        = df_sem_ready,
      always_keep = c("doy_c", "doy_c2"),
      aic_improve = aic_improve
    )

    rhs_swc_terms <- sel_swc$selected_terms

    # Response equation: y ~ swc + terms_resp + (1|boxlabel) + (1|tree_id)
    sel_resp <- backward_select_lmer_terms(
      response    = "y",
      fixed_terms = c("swc", terms_list$terms_resp),
      random_part = "(1|boxlabel) + (1|tree_id)",
      data        = df_sem_ready,
      always_keep = c("swc", "doy_c", "doy_c2"),
      aic_improve = aic_improve
    )

    rhs_resp_terms <- setdiff(sel_resp$selected_terms, "swc") # swc is explicit in formula below
  } else {
    # old behavior: keep all terms (no selection)
    rhs_swc_terms <- terms_list$terms_swc
    rhs_resp_terms <- terms_list$terms_resp
  }

  # turn term vectors into RHS strings
  rhs_swc <- paste(rhs_swc_terms, collapse = " + ")
  rhs_resp <- paste(c("swc", rhs_resp_terms), collapse = " + ")

  # build final formulas (you can use REML = TRUE here if you want)
  fml_swc <- as.formula(
    paste("swc ~", rhs_swc, "+ (1|boxlabel)")
  )

  fml_resp <- as.formula(
    paste("y ~", rhs_resp, "+ (1|boxlabel) + (1|tree_id)")
  )

  mod_swc <- lme4::lmer(fml_swc, data = df_sem_ready)
  mod_resp <- lme4::lmer(fml_resp, data = df_sem_ready)

  list(
    mod_swc = mod_swc,
    mod_resp = mod_resp,
    fml_swc = fml_swc,
    fml_resp = fml_resp,
    rhs_swc_terms = rhs_swc_terms,
    rhs_resp_terms = rhs_resp_terms,
    used_factors = terms_list$used_factors
  )
}

# # ....................................................------------------------
# # Temporal? -----------------------------------------------------------------------------------
#
# plot_swc_event_windows <- function(df, robinia_levels = c("without-robinia", "with-robinia"),
#                                    event_windows = list(
#                                      c(as.Date("2025-06-20"), as.Date("2025-07-02")),
#                                      c(as.Date("2025-08-12"), as.Date("2025-08-20"))
#                                    )) {
#
#   # Assign period label
#   df2 <- df %>%
#     mutate(period = case_when(
#       date <  event_windows[[1]][1] ~ "pre",
#       between(date, event_windows[[1]][1], event_windows[[1]][2]) ~ "during1",
#       between(date, event_windows[[2]][1], event_windows[[2]][2]) ~ "during2",
#       date >  event_windows[[2]][2] ~ "post",
#       TRUE ~ NA_character_
#     )) %>%
#     filter(!is.na(period))
#
#   # Summaries
#   df_sum <- df2 %>%
#     group_by(robinia, period) %>%
#     summarise(
#       mean_swc = mean(swc, na.rm = TRUE),
#       se_swc   = sd(swc, na.rm = TRUE) / sqrt(n()),
#       .groups = "drop"
#     ) %>%
#     mutate(period = factor(period, levels = c("pre", "post")))
#
#   ggplot(df_sum, aes(period, mean_swc, color = robinia, group = robinia)) +
#     geom_point(size = 3, position = position_dodge(width = 0.3)) +
#     geom_errorbar(
#       aes(ymin = mean_swc - se_swc, ymax = mean_swc + se_swc),
#       width = 0.2, position = position_dodge(width = 0.3)
#     ) +
#     geom_line(position = position_dodge(width = 0.3)) +
#     scale_color_manual(values = c("without-robinia" = "black", "with-robinia" = "firebrick")) +
#     theme_bw(base_size = 12) +
#     labs(
#       title = "SWC response to extreme events, grouped by robinia presence",
#       x = "Period",
#       y = "Mean SWC (%)",
#       color = "Robinia"
#     )
# }
#
# res <- run_sem_for_trait(
#   type               = "tree",
#   data_name          = "growth",
#   resp_var           = "diameter",
#   species            = "fagus",
#   soil_type          = "both",
#   include_interaction = TRUE,
#   scale_all_numeric  = FALSE
# )
#
# plot_swc_event_windows(res$data)
#
#
#
# # Both droughts -------------------------------------------------------------------------------
#
# plot_swc_pre_post_by_window <- function(
#     df,
#     event_windows,
#     swc_var = "swc",
#     robinia_var = "robinia"
# ) {
#   # deps
#   require(dplyr)
#   require(ggplot2)
#
#   stopifnot("date" %in% names(df))
#   stopifnot(swc_var %in% names(df))
#   stopifnot(robinia_var %in% names(df))
#
#   # basic prep
#   df <- df %>%
#     mutate(date = as.Date(date)) %>%
#     arrange(date)
#
#   df$swc_col     <- df[[swc_var]]
#   df$robinia_col <- df[[robinia_var]]
#
#   all_dates <- sort(unique(df$date))
#
#   # ---------- build table of pre/post dates per window ----------
#   # event_windows is a list of c(start, end)
#   ew_tbl <- lapply(seq_along(event_windows), function(i) {
#     win <- event_windows[[i]]
#     start_i <- as.Date(win[1])
#     end_i   <- as.Date(win[2])
#
#     pre_dates  <- all_dates[all_dates <  start_i]
#     post_dates <- all_dates[all_dates >  end_i]
#
#     if (!length(pre_dates) || !length(post_dates)) {
#       warning("Event ", i, " has no valid pre or post date - skipping.")
#       return(NULL)
#     }
#
#     tibble::tibble(
#       event_id   = i,
#       pre_date   = max(pre_dates),
#       post_date  = min(post_dates),
#       pre_label  = paste0("pre_",  i),
#       post_label = paste0("post_", i)
#     )
#   }) %>%
#     bind_rows()
#
#   if (nrow(ew_tbl) == 0L) {
#     stop("No valid pre/post dates could be derived from event_windows.")
#   }
#
#   # ---------- tag all rows on those pre/post dates ----------
#   df$event_period <- NA_character_
#
#   for (i in seq_len(nrow(ew_tbl))) {
#     pre_d  <- ew_tbl$pre_date[i]
#     post_d <- ew_tbl$post_date[i]
#
#     df$event_period[df$date == pre_d]  <- ew_tbl$pre_label[i]
#     df$event_period[df$date == post_d] <- ew_tbl$post_label[i]
#   }
#
#   # ---------- summarise ----------
#   df_sum <- df %>%
#     filter(!is.na(event_period)) %>%
#     group_by(robinia_col, event_period) %>%
#     summarise(
#       mean_swc = mean(swc_col, na.rm = TRUE),
#       se_swc   = sd(swc_col,  na.rm = TRUE) / sqrt(dplyr::n()),
#       .groups  = "drop"
#     )
#
#   # custom order: pre_1, post_1, pre_2, post_2, ...
#   event_ids    <- sort(unique(ew_tbl$event_id))
#   event_levels <- as.vector(rbind(
#     paste0("pre_",  event_ids),
#     paste0("post_", event_ids)
#   ))
#
#   df_sum <- df_sum %>%
#     mutate(
#       event_period = factor(event_period, levels = event_levels),
#       robinia_col  = factor(
#         robinia_col,
#         levels = c("without-robinia", "with-robinia")
#       )
#     )
#
#   # ---------- plot ----------
#   p <- ggplot(
#     df_sum,
#     aes(x = event_period,
#         y = mean_swc,
#         colour = robinia_col,
#         group = robinia_col)
#   ) +
#     geom_line(linewidth = 0.8, position = position_dodge(width = 0.1)) +
#     geom_point(size = 2.6, position = position_dodge(width = 0.1)) +
#     geom_errorbar(
#       aes(ymin = mean_swc - se_swc,
#           ymax = mean_swc + se_swc),
#       width = 0.15,
#       position = position_dodge(width = 0.1)
#     ) +
#     scale_color_manual(
#       values = c("without-robinia" = "black",
#                  "with-robinia"    = "firebrick"),
#       name   = "Robinia"
#     ) +
#     labs(
#       x = "Event stage",
#       y = "Mean SWC",
#       title    = "SWC response to extreme events",
#       subtitle = "Pre/post for each event window, separated by robinia presence"
#     ) +
#     theme_minimal(base_size = 13) +
#     theme(
#       panel.grid.minor = element_blank()
#     )
#
#   attr(p, "summary_table") <- df_sum
#   p
# }
#
# event_windows <- list(
#   c(as.Date("2025-06-20"), as.Date("2025-07-02")),
#   c(as.Date("2025-08-12"), as.Date("2025-08-20"))
# )
#
# p_events <- plot_swc_pre_post_by_window(
#   df            = df_sem_ready,   # or whatever df you use
#   event_windows = event_windows,
#   swc_var       = "swc",
#   robinia_var   = "robinia"
# )
#
# p_events
# attr(p_events, "summary_table")
#
# # MARGINAL EFFECT OF ROB ON SWC ---------------------------------------------------------------
#
# plot_marginal_swc_effects <- function(
#     mod_swc,
#     df_sem_ready,
#     event_windows,
#     robinia_var = "robinia"
# ) {
#   require(emmeans)
#   require(dplyr)
#   require(ggplot2)
#
#   # -----------------------------------------------
#   # 1. Determine representative DOY_c values
#   # -----------------------------------------------
#   # pick midpoints of pre, during, post windows
#   periods <- lapply(event_windows, function(w) {
#     w <- as.Date(w)
#     tibble(
#       pre_doy    = as.numeric(format(w[1] - 5, "%j")),
#       during_doy = as.numeric(format((w[1] + w[2]) / 2, "%j")),
#       post_doy   = as.numeric(format(w[2] + 5, "%j"))
#     )
#   }) |> bind_rows()
#
#   # take unique time points
#   doy_values <- unique(unlist(periods))
#
#   # center for model
#   doy_center <- mean(df_sem_ready$doy)
#   doy_c_vals <- doy_values - doy_center
#
#   # -----------------------------------------------
#   # 2. Build prediction grid
#   # -----------------------------------------------
#   pred_grid <- expand.grid(
#     robinia       = levels(df_sem_ready[[robinia_var]]),
#     extreme_event = c("no", "yes"),
#     precipitation = levels(df_sem_ready$precipitation),
#     soiltype      = levels(df_sem_ready$soiltype),
#     culture       = levels(df_sem_ready$culture),
#     doy_c         = doy_c_vals
#   ) %>%
#     mutate(
#       doy_c2 = doy_c^2,
#       period = factor(rep(c("pre", "during", "post"), each = n() / length(doy_c_vals)),
#                       levels = c("pre", "during", "post"))
#     )
#
#   # -----------------------------------------------
#   # 3. Get marginal predictions
#   # -----------------------------------------------
#   emm <- emmeans(
#     mod_swc,
#     specs = ~ robinia * extreme_event * period,
#     at = list(
#       doy_c  = doy_c_vals,
#       doy_c2 = doy_c_vals^2
#     ),
#     data = pred_grid,
#     type = "response"
#   )
#
#   emm_df <- as.data.frame(emm)
#
#   # -----------------------------------------------
#   # 4. Plot
#   # -----------------------------------------------
#   p <- ggplot(
#     emm_df,
#     aes(
#       x = period,
#       y = emmean,
#       colour = robinia,
#       group = robinia
#     )
#   ) +
#     geom_line(linewidth = 1) +
#     geom_point(size = 3) +
#     geom_errorbar(
#       aes(ymin = lower.CL, ymax = upper.CL),
#       width = 0.15
#     ) +
#     facet_wrap(~ extreme_event, ncol = 2) +
#     scale_color_manual(values = c("black", "firebrick")) +
#     labs(
#       title = "Marginal SWC estimates",
#       subtitle = "Adjusted for drought, soil type, culture, and temporal trend",
#       y = "Predicted SWC (marginal)",
#       x = "Event period",
#       colour = "Robinia"
#     ) +
#     theme_minimal(base_size = 13)
#
#   attr(p, "emm_table") <- emm_df
#   p
# }
#
#
# res <- run_sem_for_trait(
#   type               = "tree",
#   data_name          = "growth",
#   resp_var           = "diameter",
#   species            = "fagus",
#   soil_type          = "both",
#   include_interaction = TRUE,
#   scale_all_numeric  = TRUE
# )
#
# p_swc <- plot_marginal_swc_effects(
#   mod_swc      = res$mod_swc,
#   df_sem_ready = res$data,
#   event_windows = list(
#     c(as.Date("2025-06-20"), as.Date("2025-07-02")),
#     c(as.Date("2025-08-12"), as.Date("2025-08-20"))
#   )
# )
#
# p_swc
#
#
#
# ## New approach --------------------------------------------------------------------------------
# # 1. Start from your SEM data
# df <- res$data %>%
#   mutate(date = as.Date(date))
#
# # 2. Build a period variable using your event windows
# event_windows <- list(
#   c(as.Date("2025-06-20"), as.Date("2025-07-02")+7),
#   c(as.Date("2025-08-12"), as.Date("2025-08-20")+7)
# )
#
# df$period <- NA_character_
#
# all_dates <- sort(unique(df$date))
#
# for (i in seq_along(event_windows)) {
#   win <- event_windows[[i]]
#   start_i <- win[1]
#   end_i   <- win[2]
#
#   pre_date  <- max(all_dates[all_dates <  start_i])
#   post_date <- min(all_dates[all_dates >  end_i])
#   during_dates <- all_dates[all_dates >= start_i & all_dates <= end_i]
#
#   df$period[df$date == pre_date]         <- paste0("pre_", i)
#   df$period[df$date %in% during_dates]   <- paste0("during_", i)
#   df$period[df$date == post_date]        <- paste0("post_", i)
# }
#
# df <- df %>% filter(!is.na(period))
# df$period <- factor(df$period,
#                     levels = sort(unique(df$period)))  # e.g. pre_1, during_1, post_1, ...
#
# # 3. Marginal SWC per robinia × extreme × period
# emm <- emmeans::emmeans(
#   mod_swc,
#   specs = ~ robinia * extreme_event | period,
#   type  = "response"
# )
#
# emm_df <- as.data.frame(emm)
#
# # 4. Plot
# ggplot(emm_df,
#        aes(x = period, y = emmean,
#            colour = robinia, group = robinia)) +
#   geom_line(linewidth = 1) +
#   geom_point(size = 3) +
#   geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
#                 width = 0.15) +
#   facet_wrap(~ extreme_event) +
#   theme_minimal()
#
#
#
# # 🚨 ANOTHER APPROACH ----------------------------------------------------------------------------
# library(dplyr)
# library(emmeans)
# library(ggplot2)
# library(tibble)
#
# # ------------------------------------------------------------------
# # 1) YOUR EXISTING SEM WRAPPER (unchanged, just for context)
# # ------------------------------------------------------------------
# # run_sem_for_trait <- function(type = "tree",
# #                               data_name,
# #                               resp_var,
# #                               species,
# #                               soil_type = "both",
# #                               include_interaction = TRUE,
# #                               scale_all_numeric = TRUE) {
# #   df_pre <- prepare_df_generic(
# #     type         = type,
# #     data_name    = data_name,
# #     resp_var     = resp_var,
# #     species_keep = species,
# #     add_covars   = FALSE,
# #     covars_fun   = NULL,
# #     soil_type    = soil_type
# #   )
# #
# #   df_sem_ready <- prepare_sem_data(
# #     df_prepared        = df_pre,
# #     scale_all_numeric  = scale_all_numeric
# #   )
# #
# #   mods <- fit_sem_models(
# #     df_sem_ready       = df_sem_ready,
# #     include_interaction = include_interaction
# #   )
# #
# #   mod_swc  <- mods$mod_swc
# #   mod_resp <- mods$mod_resp
# #
# #   sem_mod <- piecewiseSEM::psem(
# #     mod_swc,
# #     mod_resp,
# #     data = df_sem_ready
# #   )
# #
# #   effects_main <- extract_sem_effects_with_se(mod_swc, mod_resp)
# #   effects_int  <- if (include_interaction) {
# #     extract_interaction_effect(mod_swc, mod_resp)
# #   } else {
# #     NULL
# #   }
# #
# #   p <- plot_sem_graph(
# #     effects_main       = effects_main,
# #     effects_int        = effects_int,
# #     sem_mod            = sem_mod,
# #     resp_var           = resp_var,
# #     species            = species,
# #     soil_type          = soil_type,
# #     include_interaction = include_interaction
# #   )
# #
# #   list(
# #     data        = df_sem_ready,
# #     mod_swc     = mod_swc,
# #     mod_resp    = mod_resp,
# #     sem         = sem_mod,
# #     effects     = effects_main,
# #     effects_int = effects_int,
# #     plot        = p
# #   )
# # }
#
# # ------------------------------------------------------------------
# # 2) TAG OBSERVED DATES AS pre/during/post FOR EACH EXTREME EVENT
# # ------------------------------------------------------------------
#
# add_event_periods <- function(df, event_windows) {
#   df <- df %>% mutate(date = as.Date(date))
#
#   all_dates <- sort(unique(df$date))
#   df$period <- NA_character_
#
#   for (i in seq_along(event_windows)) {
#     win     <- as.Date(event_windows[[i]])
#     start_i <- win[1]
#     end_i   <- win[2]
#
#     pre_date     <- max(all_dates[all_dates <  start_i], na.rm = TRUE)
#     post_date    <- min(all_dates[all_dates >  end_i],   na.rm = TRUE)
#     during_dates <- all_dates[all_dates >= start_i & all_dates <= end_i]
#
#     df$period[df$date == pre_date]       <- paste0("pre_",    i)
#     df$period[df$date %in% during_dates] <- paste0("during_", i)
#     df$period[df$date == post_date]      <- paste0("post_",   i)
#   }
#
#   df <- df %>% filter(!is.na(period))
#
#   # pre_1, during_1, post_1, pre_2, during_2, post_2, ...
#   k <- length(event_windows)
#   levs <- unlist(lapply(seq_len(k), function(i) {
#     paste0(c("pre_", "during_", "post_"), i)
#   }))
#   levs <- levs[levs %in% unique(df$period)]
#
#   df %>%
#     mutate(period = factor(period, levels = levs))
# }
#
# # ------------------------------------------------------------------
# # 3) MARGINAL SWC EFFECTS PER PERIOD (safe against 1-level factors)
# # ------------------------------------------------------------------
#
# get_swc_marginals_by_period <- function(mod_swc, df_with_period) {
#   df_use  <- df_with_period %>% filter(!is.na(period))
#   periods <- levels(df_use$period)
#
#   out <- lapply(periods, function(per) {
#     df_sub <- df_use %>%
#       filter(period == per) %>%
#       droplevels()
#
#     # if robinia or extreme_event collapses to 1 level, skip
#     if (nlevels(df_sub$robinia) < 2L || nlevels(df_sub$extreme_event) < 2L) {
#       return(tibble())
#     }
#
#     emm <- tryCatch(
#       emmeans::emmeans(
#         mod_swc,
#         specs = ~ robinia * extreme_event,
#         data  = df_sub,
#         type  = "response"
#       ),
#       error = function(e) NULL
#     )
#
#     if (is.null(emm)) {
#       return(tibble())
#     }
#
#     as.data.frame(emm) %>%
#       mutate(period = per)
#   })
#
#   res <- bind_rows(out)
#
#   if (!nrow(res)) {
#     warning("No valid periods with ≥2 levels in robinia and extreme_event.")
#     return(res)
#   }
#
#   res$period <- factor(res$period, levels = periods)
#   res
# }
#
# # ------------------------------------------------------------------
# # 4) PLOT MARGINAL SWC OVER pre/during/post × robinia × extreme_event
# # ------------------------------------------------------------------
#
# plot_swc_period_effects <- function(emm_df,
#                                     title = "Marginal SWC (model-adjusted) across pre/during/post periods") {
#   if (!nrow(emm_df)) {
#     stop("emm_df is empty - nothing to plot.")
#   }
#
#   ggplot(
#     emm_df,
#     aes(
#       x      = period,
#       y      = emmean,
#       group  = robinia,
#       colour = robinia
#     )
#   ) +
#     geom_line(linewidth = 1) +
#     geom_point(size = 3) +
#     geom_errorbar(
#       aes(ymin = lower.CL, ymax = upper.CL),
#       width = 0.15
#     ) +
#     facet_wrap(~ extreme_event) +
#     labs(
#       title  = title,
#       y      = "Predicted SWC",
#       x      = NULL,
#       colour = "Robinia"
#     ) +
#     theme_minimal(base_size = 12)
# }
#
# # ------------------------------------------------------------------
# # 5) END-TO-END EXAMPLE USING run_sem_for_trait()
# # ------------------------------------------------------------------
#
# # User input --------------------------------------------------------
# # (change these when you switch trait/species etc.)
#
# type_in        <- "tree"
# data_name_in   <- "growth"
# resp_var_in    <- "diameter_inc_t0"
# species_in     <- "quercus"
# soil_type_in   <- "both"
# include_int_in <- TRUE
# scale_num_in   <- TRUE
#
# # define your event windows (your preferred version)
# event_windows <- list(
#   c(as.Date("2025-06-20"), as.Date("2025-07-02")),
#   c(as.Date("2025-08-12"), as.Date("2025-08-20"))
# )
#
# # Run SEM as usual
# sem_out <- run_sem_for_trait(
#   type               = type_in,
#   data_name          = data_name_in,
#   resp_var           = resp_var_in,
#   species            = species_in,
#   soil_type          = soil_type_in,
#   include_interaction = include_int_in,
#   scale_all_numeric  = scale_num_in
# )
#
# # Extract SEM data + SWC submodel
# df_sem_ready <- sem_out$data
# mod_swc      <- sem_out$mod_swc
#
# # Tag pre/during/post periods based on *observed dates*
# df_sem_periods <- add_event_periods(df_sem_ready, event_windows)
#
# # Compute marginal SWC per period × robinia × extreme_event
# emm_df <- get_swc_marginals_by_period(mod_swc, df_sem_periods)
#
# # Plot
# p_swc_period <- plot_swc_period_effects(
#   emm_df,
#   title = paste0(
#     "SWC response to extreme events (",
#     species_in, ", ", soil_type_in, ", resp: ", resp_var_in, ")"
#   )
# )
#
# p_swc_period
