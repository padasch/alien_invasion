# Phenology data preparation shared by the production bootstrap model.

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
      precipitation = factor(
        .data$precipitation,
        levels = alinv_factor_levels("precipitation")
      ),
      robinia = factor(.data$robinia, levels = alinv_factor_levels("robinia")),
      culture = factor(.data$culture, levels = alinv_factor_levels("culture")),
      soiltype = alinv_relevel_soiltype(.data$soiltype)
    ) %>%
    droplevels()
}
