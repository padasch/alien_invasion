# Packages
library(lme4)
library(piecewiseSEM)
library(dplyr)
library(emmeans)
library(multcomp)

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
      soiltype = factor(soiltype, levels = c("inoc-beech", "inoc-robinia")),

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

extract_sem_effects_with_se <- function(mod_swc, mod_resp) {
  beta_swc <- fixef(mod_swc)
  beta_resp <- fixef(mod_resp)
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

  factors <- c("robinia", "precipitation", "culture", "soiltype", "extreme_event")

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

extract_interaction_effect <- function(mod_swc, mod_resp) {
  beta_swc <- fixef(mod_swc)
  beta_resp <- fixef(mod_resp)
  vcov_swc <- as.matrix(vcov(mod_swc))
  vcov_resp <- as.matrix(vcov(mod_resp))

  # SWC → response
  b <- beta_resp["swc"]
  var_b <- vcov_resp["swc", "swc"]
  se_b <- sqrt(var_b)

  # treatment factors we allow to interact
  treat_factors <- c("robinia", "precipitation", "culture", "soiltype", "extreme_event")

  out <- list()

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


plot_sem_graph <- function(effects_main,
                           effects_int = NULL,
                           sem_mod,
                           resp_var,
                           species,
                           soil_type,
                           include_interaction = TRUE) {
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

  resp_node <- resp_var

  # automatic node layout
  # fixed order of main predictors on the left
  main_nodes_order <- c("extreme_event", "robinia", "precipitation", "soiltype", "culture")

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
      est  = a,
      se   = se_a,
      p    = p_a,
      path = "a (via swc)"
    )

  edges_c <- effects_main %>%
    dplyr::transmute(
      from = factor,
      to   = resp_node,
      est  = c_direct,
      se   = se_c,
      p    = p_c,
      path = "c' (direct)"
    )

  b_row <- effects_main %>% dplyr::slice(1)
  edge_b <- b_row %>%
    dplyr::transmute(
      from = "swc",
      to   = resp_node,
      est  = b,
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
        est  = a,
        se   = se_a,
        p    = p_a,
        path = "interaction via swc"
      )

    edges_int_c <- effects_int %>%
      dplyr::transmute(
        from = factor,
        to   = resp_node,
        est  = c_direct,
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
  p_sig <- 0.1
  edges <- edges %>%
    dplyr::filter(!is.na(p) & p < p_sig)

  # Identify nodes that appear in at least one significant edge
  used_nodes <- unique(c(edges$from, edges$to))

  # Main nodes we always keep (even if they have no sig edges)
  main_nodes_order <- c(
    "extreme_event", "robinia", "precipitation",
    "soiltype", "culture", "swc", resp_node
  )

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

  max_x <- max(nodes$x) + 0.5
  max_y <- max(nodes$y) + 0.7
  min_y <- min(nodes$y) - 0.7

  ggplot() +
    geom_segment(
      data = edges_plot,
      aes(
        x = x_from, y = y_from,
        xend = x_end, yend = y_end,
        size = lw,
        colour = col_dir,
        alpha = alpha
      ),
      arrow = arrow(length = unit(0.5, "cm"), type = "open"),
      lineend = "round"
    ) +
    geom_label(
      data = nodes,
      aes(x = x, y = y, label = node),
      size = 3.5,
      label.size = 0.3,
      label.r = unit(0.15, "lines"),
      fill = "white",
      label.padding = unit(0.15, "lines")
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

run_sem_for_trait <- function(type = "tree",
                              data_name,
                              resp_var,
                              species,
                              soil_type = "both",
                              include_interaction = TRUE,
                              scale_all_numeric = TRUE,
                              do_rfe = FALSE,
                              aic_improve = 2,
                              force_run = FALSE) {
  # --- 0) caching: save / load SEM run for today ---
  today_str <- format(Sys.Date(), "%Y-%m-%d")

  out_dir <- file.path("output", today_str, "model-sem")
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  resp_tag <- if (!is.null(resp_var)) resp_var else "default"
  int_tag <- if (isTRUE(include_interaction)) "int" else "noInt"
  scale_tag <- if (isTRUE(scale_all_numeric)) "scaled" else "unscaled"
  rfe_tag <- if (isTRUE(do_rfe)) paste0("rfeAIC", aic_improve) else "noRFE"

  file_name <- paste0(
    "sem-",
    type, "-",
    data_name, "-",
    resp_tag, "-",
    species, "-",
    soil_type, "-",
    int_tag, "-",
    scale_tag, "-",
    rfe_tag,
    ".rds"
  )

  cache_path <- file.path(out_dir, file_name)

  if (file.exists(cache_path) && !force_run) {
    message("Loading cached SEM results from: ", cache_path)
    return(readRDS(cache_path))
  }

  # 1) Load analysis-ready data in your project’s canonical way
  df_pre <- prepare_df_generic(
    type         = type,
    data_name    = data_name,
    resp_var     = resp_var,
    species_keep = species,
    add_covars   = FALSE, # SEM doesn’t use the extra temporal covars
    covars_fun   = NULL,
    soil_type    = soil_type
  )

  # 2) Prepare specifically for SEM (generic y, scaling, etc.)
  df_sem_ready <- prepare_sem_data(
    df_prepared       = df_pre,
    scale_all_numeric = scale_all_numeric
  )

  # 3) Fit SEM submodels
  mods <- fit_sem_models(
    df_sem_ready        = df_sem_ready,
    include_interaction = include_interaction,
    do_rfe              = do_rfe,
    aic_improve         = aic_improve
  )

  mod_swc <- mods$mod_swc
  mod_resp <- mods$mod_resp

  # 4) Combine into piecewise SEM
  sem_mod <- piecewiseSEM::psem(
    mod_swc,
    mod_resp,
    data = df_sem_ready
  )

  # 5) Extract effects
  effects_main <- extract_sem_effects_with_se(mod_swc, mod_resp)
  effects_int <- if (include_interaction) {
    extract_interaction_effect(mod_swc, mod_resp)
  } else {
    NULL
  }

  # 6) Plot graph with full metadata (species, soil_type, resp_var)
  p <- plot_sem_graph(
    effects_main = effects_main,
    effects_int = effects_int,
    sem_mod = sem_mod,
    resp_var = resp_var,
    species = species,
    soil_type = soil_type,
    include_interaction = include_interaction
  )

  result <- list(
    data        = df_sem_ready,
    mod_swc     = mod_swc,
    mod_resp    = mod_resp,
    sem         = sem_mod,
    effects     = effects_main,
    effects_int = effects_int,
    plot        = p
  )

  saveRDS(result, cache_path)
  message("Saved SEM results to: ", cache_path)

  result
}
# ....................................................------------------------

build_sem_terms <- function(df, include_interaction = TRUE) {
  # all potential treatment factors
  treat_factors_all <- c(
    "precipitation",
    "robinia",
    "soiltype",
    "culture",
    "extreme_event"
  )

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
  do_rfe = FALSE,
  aic_improve = 2
) {
  # build full set of main and interaction terms
  terms_list <- build_sem_terms(df_sem_ready, include_interaction = include_interaction)

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
    rhs_resp_terms = rhs_resp_terms
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
