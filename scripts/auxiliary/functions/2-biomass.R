# 2-biomass.R
# Biomass wrangling shared by publication and bootstrap workflows.

wrangle_tree_biomass <- function(fp, sheet = "Biomass") {
  meta_tree <- get_meta("tree") %>%
    mutate(tree_id = as.character(tree_id)) %>%
    select(treelabel, tree_id, boxlabel, soiltype, culture, robinia) %>%
    distinct(treelabel, .keep_all = TRUE)

  read_excel(fp, sheet = sheet) %>%
    alinv_drop_empty_cols() %>%
    alinv_clean_names_df() %>%
    mutate(across(where(is.character), ~ tolower(trimws(.x)))) %>%
    mutate(across(where(is.character), ~ gsub("_", "-", .x))) %>%
    mutate(
      treelabel = sub("^(fagus|quercus)-", "", label),
      root_biomass = parse_number(as.character(root_biomass), na = c("", "NA", "na")),
      shoot_biomass = parse_number(as.character(shoot_biomass), na = c("", "NA", "na"))
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
