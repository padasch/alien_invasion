# Central analysis configuration for factor baselines, response labels, and
# treatment display wording used across notebooks and exported outputs.

ALINV_FACTOR_LEVELS <- list(
  precipitation = c("control", "drought"),
  culture = c("mono", "mixed"),
  robinia = c("without-robinia", "with-robinia"),
  soiltype = c("inoc-robinia", "inoc-beech"),
  extreme_event = c("no", "yes"),
  species = c("fagus", "quercus", "robinia")
)

ALINV_LEVEL_LABELS <- list(
  robinia = c(
    `without-robinia` = "without robinia",
    `with-robinia` = "with robinia"
  ),
  soiltype = c(
    `inoc-robinia` = "wetter soil (robinia soil)",
    `inoc-beech` = "drier soil (beech soil)",
    inoc_robinia = "wetter soil (robinia soil)",
    inoc_beech = "drier soil (beech soil)"
  )
)

ALINV_RESPONSE_LABELS <- c(
  chl = "Chlorophyll",
  condition = "Condition",
  height = "Height",
  diameter = "Diameter",
  volume = "Volume",
  height_inc_t0 = "Height inc. t0",
  diameter_inc_t0 = "Diameter inc. t0",
  volume_inc_t0 = "Volume inc. t0",
  height_inc_t0_rel = "Height rel. inc t0",
  diameter_inc_t0_rel = "Diameter rel. inc t0",
  volume_inc_t0_rel = "Volume rel. inc t0",
  height_inc_phase_abs = "Height inc. phase",
  diameter_inc_phase_abs = "Diameter inc. phase",
  volume_inc_phase_abs = "Volume inc. phase",
  height_inc_phase_rel = "Height rel. inc phase",
  diameter_inc_phase_rel = "Diameter rel. inc phase",
  volume_inc_phase_rel = "Volume rel. inc phase",
  qy = "Quantum yield",
  remaining_green = "Senescence",
  chlavg = "Senescence chlorophyll",
  stage = "Phenology stage"
)

ALINV_SCENARIO_GRID <- tibble::tribble(
  ~scenario_label, ~soil_filter, ~include_soil_treatment,
  "drier soil (beech soil)", "inoc-beech", FALSE,
  "wetter soil (robinia soil)", "inoc-robinia", FALSE,
  "both soils (with soil as treatment)", "both", TRUE,
  "both soils (without soil as treatment)", "both", FALSE
)

ALINV_TREATMENT_CONFIG <- tibble::tribble(
  ~effect, ~baseline_level, ~treatment_level, ~short_label, ~temporal_label, ~heatmap_label, ~contrast_label, ~plot_order,
  "precipitation", "control", "drought", "Drought", "Precipitation (control -> drought)", "Precipitation: control -> drought", "drought - control", 1L,
  "robinia", "without-robinia", "with-robinia", "With robinia", "Robinia (without -> with)", "Robinia: without -> with", "with robinia - without robinia", 2L,
  "culture", "mono", "mixed", "Mixed culture", "Culture (mono -> mixed)", "Culture: mono -> mixed", "mixed - mono", 3L,
  "soiltype", "inoc-robinia", "inoc-beech", "Drier soil", "Soil (wetter robinia soil -> drier beech soil)", "Soil: wetter soil (robinia soil) -> drier soil (beech soil)", "drier soil - wetter soil", 4L,
  "extreme_event", "no", "yes", "Extreme event", "Extreme event (no -> yes)", "Extreme event: no -> yes", "yes - no", 5L
)

alinv_factor_levels <- function(name) {
  ALINV_FACTOR_LEVELS[[as.character(name)]] %||% NULL
}

alinv_level_labels <- function(name) {
  ALINV_LEVEL_LABELS[[as.character(name)]] %||% NULL
}

alinv_response_labels <- function() {
  ALINV_RESPONSE_LABELS
}

alinv_response_label <- function(resp_var) {
  labels <- alinv_response_labels()
  out <- unname(labels[as.character(resp_var)])
  ifelse(is.na(out) | !nzchar(out), as.character(resp_var), out)
}

alinv_treatment_label_map <- function(style = c("temporal", "heatmap", "short", "contrast")) {
  style <- match.arg(style)
  col_name <- switch(
    style,
    temporal = "temporal_label",
    heatmap = "heatmap_label",
    short = "short_label",
    contrast = "contrast_label"
  )
  stats::setNames(ALINV_TREATMENT_CONFIG[[col_name]], ALINV_TREATMENT_CONFIG$effect)
}
