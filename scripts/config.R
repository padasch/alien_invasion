# Central analysis configuration for factor baselines, response labels, and
# treatment display wording used by the active analyses and figures.

# Central growth proxy options.
#
# Notes:
# - The public response name remains `volume` for compatibility.
# - The default method is the Annighöfer et al. (2016) Table 4 RCD²H allometry.
# - An alternative method from Zianis et al. (2005) is available and can be
#   selected by setting `ALINV_VOLUME_PROXY_METHOD=zianis_2005` when running
#   the data-cleaning pipeline.
# - The experimental oak is Quercus ilex biologically, but the repo treats the
#   oak group generically as `quercus`. We therefore use a configurable
#   surrogate species from Annighöfer et al.
# - The default surrogate for the oak proxy is Quercus robur because its published
#   calibration range covers the observed RCD²H values in this dataset.

ALINV_VOLUME_PROXY_METHODS <- c("annighofer", "zianis_2005")
ALINV_VOLUME_PROXY_DEFAULT <- "annighofer"
ALINV_VOLUME_PROXY_VERSION <- c(
  annighofer = "annighofer_table4_v1",
  zianis_2005 = "zianis_2005_v1"
)

ALINV_VOLUME_ALLOMETRY <- list(
  source = list(
    citation = paste(
      "Annighofer, P., Ameztegui, A., Ammer, C. et al.",
      "Species-specific and generic biomass equations for seedlings and saplings",
      "of European tree species. Eur J Forest Res 135, 313-329 (2016).",
      "https://doi.org/10.1007/s10342-016-0937-z"
    ),
    equation = "proxy = beta1 * (diameter_mm^2 * height_cm)^beta2",
    table = "Table 4"
  ),
  oak_surrogate = "quercus_robur",
  specs = tibble::tribble(
    ~repo_species, ~source_species_key, ~source_species_label, ~beta1, ~beta2, ~rcd2h_min, ~rcd2h_max, ~notes,
    "fagus", "fagus_sylvatica", "Fagus sylvatica", 0.62342, 0.87409, 0, 132559, "Species-specific Table 4 equation.",
    "quercus", "quercus_robur", "Quercus robur", 0.67311, 0.85202, 2, 65307, "Current default surrogate for repo `quercus` because the published range covers the observed data.",
    "quercus", "quercus_petraea", "Quercus petraea", 0.52740, 0.81213, 1, 16366, "Alternative oak surrogate kept switchable in config; would require extrapolation for current larger observations."
  )
)

ALINV_VOLUME_ZIANIS <- list(
  source = list(
    citation = paste(
      "Zianis, D., Muukkonen, P., Mäkipää, R. & Mencuccini, M. (2005).",
      "Biomass and stem volume equations for tree species in Europe.",
      "https://w1.aua.gr/dasologia/wp-content/uploads/sites/4/2023/01/2005_Zianis_SF.pdf"
    ),
    formula = "Species-specific equations from the linked document."
  ),
  specs = tibble::tribble(
    ~repo_species, ~source_species_key, ~source_species_label, ~equation_id, ~formula_type, ~a, ~b, ~c, ~d, ~e, ~f, ~g, ~source_page_formula, ~source_page_parameters, ~notes,
    "fagus", "fagus_sylvatica", "Fagus sylvatica", 49L, "poly6", -0.015572, 0.00290013, -7.0476e-6, 2.393e-6, -0.0013528, 3.9837e-5, NA_real_, 53, 54, "Equation 49 (Belgium): a + b*D + c*D² + d*D³ + e*H + f*D²*H (poly6).",
    "quercus", "quercus_ilex", "Quercus ilex", 201L, "d2h", 1.1909, 0.038639, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, 59, 60, "Equation 201 (Italy): a + b*D²*H. Using this as the default oak surrogate in this dataset.",
    "tilia", "tilia_cordata", "Tilia cordata", 225L, "log10_poly", 0.00004124, 1.9302, 0.0209, 0.129, -0.1903, NA_real_, NA_real_, 61, 62, "Equation 225 (Romania): a * 10^(b*log10(D) + c*log10(D)^2 + d*log10(H) + e*log10(H)^2)."
  )
)

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

alinv_volume_proxy_version <- function(method = NULL) {
  version <- ALINV_VOLUME_PROXY_VERSION[[alinv_volume_proxy_method(method)]]
  if (is.null(version)) {
    version <- ALINV_VOLUME_PROXY_VERSION[[ALINV_VOLUME_PROXY_DEFAULT]]
  }
  version
}

alinv_volume_allometry_catalog <- function() {
  ALINV_VOLUME_ALLOMETRY$specs
}

alinv_volume_allometry_source <- function() {
  ALINV_VOLUME_ALLOMETRY$source
}

alinv_volume_allometry_spec <- function(repo_species) {
  repo_species <- as.character(repo_species)
  specs <- alinv_volume_allometry_catalog()
  surrogate_key <- ALINV_VOLUME_ALLOMETRY$oak_surrogate

  if (identical(repo_species, "quercus")) {
    out <- specs[specs$repo_species == "quercus" & specs$source_species_key == surrogate_key, , drop = FALSE]
  } else {
    out <- specs[specs$repo_species == repo_species, , drop = FALSE]
    if (nrow(out) > 1L) {
      out <- out[1, , drop = FALSE]
    }
  }

  if (!nrow(out)) {
    stop("No configured volume allometry for repo species: ", repo_species, call. = FALSE)
  }

  tibble::as_tibble(out)
}
