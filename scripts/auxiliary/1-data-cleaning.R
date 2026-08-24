# 1-data-cleaning.R

# --- Setup --------------------------------------------------------------------
library(readxl)
library(tidyr)
library(readr)
library(dplyr)
library(purrr)

# Load helpers (if not already loaded)
if (!exists(".alinv_project_root", mode = "function")) {
  source_candidates <- c(
    "scripts/auxiliary/functions/_source.R",
    "_source.R"
  )
  source_path <- source_candidates[file.exists(source_candidates)][1]
  if (is.na(source_path)) {
    stop("Could not locate scripts/auxiliary/functions/_source.R for cleaning script setup.")
  }
  source(source_path)
}

# Set filepath
fp <- "./data/raw/data_20260309.xlsx"

# Print all sheet names
print(excel_sheets(fp))


# --- Meta Data -----------------------------------------------------------------

# Tree meta data
meta_tree <- read_excel(fp, sheet = "All_labels") %>%
  alinv_drop_empty_cols() %>%
  alinv_clean_names_df() %>%
  mutate(across(where(is.character), ~ tolower(trimws(.x)))) %>%
  mutate(across(where(is.character), ~ gsub("_", "-", .x))) %>%
  mutate(treelabel = sub("^[^-]+-", "", species_treelabel)) |>
  rename(tree_id = id_number) %>%
  # Remove robinia and NA in id
  filter(!is.na(species)) |> 
  # Select only tree-level information
  select(tree_id, treelabel, species, boxlabel)

# Box meta data
meta_box <- read_excel(fp, sheet = "Label_Compartment") %>%
  alinv_drop_empty_cols() %>%
  alinv_clean_names_df() %>%
  mutate(across(where(is.character), ~ tolower(trimws(.x)))) %>%
  mutate(across(where(is.character), ~ gsub("_", "-", .x)))

# Save meta
write_csv(meta_tree, "./data/interim/meta_tree.csv")
write_csv(meta_box, "./data/interim/meta_box.csv")

glimpse(meta_tree)
glimpse(meta_box)


# --- Tree-Level Data -----------------------------------------------------------

## QY ----
df <- read_excel(fp, sheet = "Fluoropen QY") %>%
  rename(boxlabel = BoxLabel) %>%
  pivot_longer(
    cols = -boxlabel,
    names_to = c("tree", "date"),
    names_sep = "_",
    values_to = "qy"
  ) %>%
  mutate(
    treelabel = paste(boxlabel, tree, sep = "-"),
    date = as.Date(date, format = "%d.%m.%Y"),
    qy = parse_number(as.character(qy))
  ) %>%
  select(treelabel, date, qy) %>%
  arrange(treelabel, date) %>%
  alinv_clean_names_df() %>%
  mutate(across(where(is.character), ~ tolower(trimws(.x))))

df <- df %>%
  left_join(meta_tree %>% select(treelabel, tree_id), by = "treelabel") %>%
  select(tree_id, date, qy) %>%
  arrange(tree_id, date) |>
  filter(!is.na(tree_id))

write_csv(df, "./data/interim/tree_quantum_yield.csv")
glimpse(df)


## Chlorophyll ----
df <- read_excel(fp, sheet = "Chlorophyll content") %>%
  rename(boxlabel = BoxLabel) %>%
  pivot_longer(
    cols = -boxlabel,
    names_to = c("tree", "date"),
    names_sep = "_",
    values_to = "chl"
  ) %>%
  mutate(
    treelabel = paste(boxlabel, tree, sep = "-"),
    date = as.Date(date, format = "%d.%m.%Y"),
    chl = parse_number(as.character(chl))
  ) %>%
  select(treelabel, date, chl) %>%
  arrange(treelabel, date) %>%
  alinv_clean_names_df() %>%
  mutate(across(where(is.character), ~ tolower(trimws(.x))))

df <- df %>%
  left_join(meta_tree %>% select(tree_id, treelabel), by = "treelabel") %>%
  select(tree_id, date, chl) %>%
  arrange(tree_id, date) |>
  filter(!is.na(tree_id))

write_csv(df, "./data/interim/tree_chlorophyll.csv")
glimpse(df)


## Leaf State ----
df <- read_excel(fp, sheet = "Tree Condition") %>%
  select(ID_number, starts_with("LeafState"), starts_with("Comment")) %>%
  mutate(across(starts_with("Comment"), as.character)) %>%
  mutate(across(starts_with("LeafState"), as.character)) %>%
  pivot_longer(
    cols = -ID_number,
    names_to = c("measure", "date"),
    names_sep = "_",
    values_to = "value"
  ) %>%
  mutate(
    date    = as.Date(date, format = "%d.%m.%Y"),
    measure = tolower(measure)
  ) %>%
  pivot_wider(
    names_from  = measure,
    values_from = value
  ) %>%
  mutate(
    tree_id = as.integer(ID_number),
    # Invert 1..4 so higher values indicate better condition (lower stress).
    condition = 5 - as.double(leafstate),
    comment = as.character(comment)
  ) %>%
  select(tree_id, date, condition, comment) %>%
  arrange(tree_id, date) |>
  filter(!is.na(tree_id))

write_csv(df, "./data/interim/tree_condition.csv")
glimpse(df)


## Senescence ----
df <- read_excel(fp, sheet = "Senescence") %>%
  select(
    ID_number,
    starts_with("%_"),
    starts_with("Chl1_"),
    starts_with("Chl2_"),
    starts_with("Comment_")
  ) %>%
  mutate(across(everything(), as.character)) %>%
  mutate(across(everything(), ~ na_if(.x, ""))) |> 
  pivot_longer(
    cols = -ID_number,
    names_to = c("measure", "date"),
    names_sep = "_",
    values_to = "value"
  ) %>%
  mutate(date = as.Date(date, format = "%d.%m.%Y")) %>%
  pivot_wider(
    names_from  = measure,
    values_from = value
  ) %>%
  mutate(
    tree_id = as.integer(ID_number),
    chl1 = parse_number(as.character(Chl1)),
    chl2 = parse_number(as.character(Chl2)),
    percent_senesced = parse_number(as.character(`%`)),
    # Derived analysis variable so larger values always indicate a more
    # favorable canopy state in downstream inferential models.
    remaining_green = 100 - percent_senesced,
    comment = as.character(Comment),
    chlavg = rowMeans(cbind(chl1, chl2), na.rm = TRUE)
  ) %>%
  select(tree_id, date, percent_senesced, remaining_green, chl1, chl2, chlavg, comment) %>%
  arrange(tree_id, date) |>
  filter(!is.na(tree_id))

write_csv(df, "./data/interim/tree_senescence.csv")
glimpse(df)


## Growth (Height and Diameter) ----

volume_proxy_method <- alinv_volume_proxy_method()
message("Using volume proxy method: ", volume_proxy_method)

df <- read_excel(fp, sheet = "Growth_Measurements_D_H") %>%
  dplyr::select(ID_number, starts_with("Diameter"), starts_with("Height")) %>%
  mutate(across(everything(), as.character)) %>%
  rename_with(tolower) %>%
  pivot_longer(
    cols = -id_number,
    names_to = c("metric", "date"),
    names_pattern = "([^_]+)_(\\d{2}\\.\\d{2}\\.\\d{4})",
    values_to = "value"
  ) %>%
  mutate(
    date   = as.Date(date, format = "%d.%m.%Y"),
    metric = tolower(metric),
    metric = gsub("\\[mm\\]", "_mm", metric) |> gsub("\\[cm\\]", "_cm", x = _)
  ) %>%
  filter(!is.na(id_number)) %>%
  pivot_wider(
    names_from  = metric,
    values_from = value
  ) %>%
  mutate(
    tree_id = as.integer(id_number),
    diameter = parse_number(diameter_mm),
    height = parse_number(height_cm)
  ) %>%
  left_join(
    meta_tree %>% transmute(tree_id = as.integer(tree_id), species),
    by = "tree_id"
  ) %>%
  mutate(
    # Public response name stays `volume` for continuity across the analysis
    # pipeline. Volume method is controlled by `ALINV_VOLUME_PROXY_METHOD`.
    volume = alinv_compute_volume_proxy(
      diameter_mm = diameter,
      height_cm = height,
      species = species,
      method = volume_proxy_method
    )
  ) %>%
  dplyr::select(tree_id, date, diameter, height, volume) %>%
  arrange(tree_id, date) |>
  filter(!is.na(tree_id)) |> 
  group_by(tree_id) |>
  mutate(
    phase = dplyr::case_when(
      lubridate::month(date) <= 6 ~ "until June",
      lubridate::month(date) <= 8 ~ "July-August",
      TRUE ~ "September+"
    ),
    phase = factor(phase, levels = c("until June", "July-August", "September+")),

    # first measurements per tree (baseline)
    first_diameter = first(diameter),
    first_height   = first(height),
    first_volume   = first(volume),
    
    # cumulative absolute increments from first measurement
    diameter_inc_t0 = diameter - first_diameter,
    height_inc_t0   = height   - first_height,
    volume_inc_t0   = volume   - first_volume,
    
    # relative size and relative cumulative increment
    diameter_rel        = diameter / first_diameter,
    height_rel          = height   / first_height,
    volume_rel          = volume   / first_volume,
    diameter_inc_t0_rel = diameter_inc_t0 / first_diameter,
    height_inc_t0_rel   = height_inc_t0   / first_height,
    volume_inc_t0_rel   = volume_inc_t0   / first_volume,

    # phase baselines: first measurement for phase 1,
    # last measurement of prior phase for phases 2 and 3
    phase_diameter_baseline = dplyr::case_when(
      phase == "until June"  ~ first_diameter,
      phase == "July-August" ~ dplyr::last(diameter[phase == "until June" & !is.na(diameter)], default = NA_real_),
      phase == "September+"  ~ dplyr::last(diameter[phase == "July-August" & !is.na(diameter)], default = NA_real_),
      TRUE ~ NA_real_
    ),
    phase_height_baseline = dplyr::case_when(
      phase == "until June"  ~ first_height,
      phase == "July-August" ~ dplyr::last(height[phase == "until June" & !is.na(height)], default = NA_real_),
      phase == "September+"  ~ dplyr::last(height[phase == "July-August" & !is.na(height)], default = NA_real_),
      TRUE ~ NA_real_
    ),
    phase_volume_baseline = dplyr::case_when(
      phase == "until June"  ~ first_volume,
      phase == "July-August" ~ dplyr::last(volume[phase == "until June" & !is.na(volume)], default = NA_real_),
      phase == "September+"  ~ dplyr::last(volume[phase == "July-August" & !is.na(volume)], default = NA_real_),
      TRUE ~ NA_real_
    ),

    # absolute increment within phase, chained to prior phase end baseline
    diameter_inc_phase_abs = diameter - phase_diameter_baseline,
    height_inc_phase_abs = height - phase_height_baseline,
    volume_inc_phase_abs = volume - phase_volume_baseline,

    # relative increment within phase, chained to prior phase end baseline
    diameter_inc_phase_rel = dplyr::if_else(
      !is.na(phase_diameter_baseline) & abs(phase_diameter_baseline) > .Machine$double.eps,
      diameter / phase_diameter_baseline - 1,
      NA_real_
    ),
    height_inc_phase_rel = dplyr::if_else(
      !is.na(phase_height_baseline) & abs(phase_height_baseline) > .Machine$double.eps,
      height / phase_height_baseline - 1,
      NA_real_
    ),
    volume_inc_phase_rel = dplyr::if_else(
      !is.na(phase_volume_baseline) & abs(phase_volume_baseline) > .Machine$double.eps,
      volume / phase_volume_baseline - 1,
      NA_real_
    ),
    
    # time between consecutive measurements [years]
    delta_t_years = as.numeric(date - dplyr::lag(date)) / 365.25,
    
    # increments between consecutive measurements
    diameter_inc_dt = diameter - dplyr::lag(diameter),
    height_inc_dt   = height   - dplyr::lag(height),
    volume_inc_dt   = volume   - dplyr::lag(volume),
    
    # absolute growth rate (AGR) between dates
    diameter_agr = diameter_inc_dt / delta_t_years,
    height_agr   = height_inc_dt   / delta_t_years,
    volume_agr   = volume_inc_dt   / delta_t_years,
    
    # relative growth rate (RGR) between dates (log-scale)
    diameter_rgr = (log(diameter) - log(dplyr::lag(diameter))) / delta_t_years,
    height_rgr   = (log(height)   - log(dplyr::lag(height)))   / delta_t_years,
    volume_rgr   = (log(volume)   - log(dplyr::lag(volume)))   / delta_t_years
  ) %>%
  ungroup()

volume_range_report <- df %>%
  left_join(
    meta_tree %>% transmute(tree_id = as.integer(tree_id), species),
    by = "tree_id"
  ) %>%
  filter(species %in% c("fagus", "quercus")) %>%
  select(species, diameter, height) %>%
  alinv_volume_allometry_range_report()

write_csv(df, "./data/interim/tree_growth.csv")
glimpse(df)
print(volume_range_report)


## Specific Leaf Area ----
df <- read_excel(fp, sheet = "SLA") %>%
  select(-Label) %>%
  rename_with(tolower) |>
  rename(tree_id = id_number) |>
  filter(!is.na(tree_id))

write_csv(df, "./data/interim/tree_sla.csv")
glimpse(df)


## Phenology ----
clean_phenology_transitions <- function(df) {
  stage_cols <- paste0("stage", 1:4)

  df_stage <- df %>%
    dplyr::filter(.data$discard == 0) %>%
    dplyr::mutate(
      tree_id = as.character(.data$tree_id),
      dplyr::across(dplyr::all_of(stage_cols), as.Date)
    )

  purrr::map_dfr(seq_len(nrow(df_stage)), function(i) {
    row_i <- df_stage[i, ]
    raw_dates <- as.Date(unlist(row_i[stage_cols], use.names = FALSE))
    clean_dates <- raw_dates
    nulled_conflict <- rep(FALSE, length(raw_dates))

    if (length(clean_dates) > 1) {
      for (stage_idx in seq(length(clean_dates) - 1, 1)) {
        later_dates <- clean_dates[(stage_idx + 1):length(clean_dates)]
        later_dates <- later_dates[!is.na(later_dates)]

        if (!length(later_dates) || is.na(clean_dates[stage_idx])) {
          next
        }

        later_min <- min(later_dates)
        if (clean_dates[stage_idx] > later_min) {
          clean_dates[stage_idx] <- as.Date(NA)
          nulled_conflict[stage_idx] <- TRUE
        }
      }
    }

    tibble(
      tree_id = row_i$tree_id,
      stage = seq_along(stage_cols),
      stage_label = paste("Stage", seq_along(stage_cols)),
      stage_date_raw = raw_dates,
      stage_date = clean_dates,
      doy = lubridate::yday(clean_dates),
      stage_order_conflict_nulled = nulled_conflict,
      stage_date_source = dplyr::case_when(
        is.na(raw_dates) ~ "missing_raw",
        nulled_conflict ~ "nulled_stage_order_conflict",
        !is.na(clean_dates) ~ "retained",
        TRUE ~ NA_character_
      )
    )
  }) %>%
    dplyr::arrange(.data$tree_id, .data$stage)
}

build_phenology_timeseries <- function(df_transitions) {
  meas_dates <- df_transitions %>%
    dplyr::filter(!is.na(.data$stage_date)) %>%
    dplyr::pull(.data$stage_date) %>%
    unique() %>%
    sort()

  tree_ids <- df_transitions %>%
    dplyr::pull(.data$tree_id) %>%
    unique()

  purrr::map_dfr(tree_ids, function(tree_id_i) {
    stage_tree <- df_transitions %>%
      dplyr::filter(.data$tree_id == tree_id_i, !is.na(.data$stage_date)) %>%
      dplyr::arrange(.data$stage_date, .data$stage)

    if (!nrow(stage_tree)) {
      return(tibble(
        tree_id = tree_id_i,
        date = meas_dates,
        stage = NA_integer_
      ))
    }

    idx <- findInterval(meas_dates, stage_tree$stage_date)
    stage_vec <- rep(NA_integer_, length(meas_dates))
    idx_valid <- idx > 0
    stage_vec[idx_valid] <- stage_tree$stage[idx[idx_valid]]

    tibble(
      tree_id = tree_id_i,
      date = meas_dates,
      stage = as.integer(stage_vec)
    )
  }) %>%
    dplyr::mutate(
      tree_id = as.character(.data$tree_id),
      date = as.Date(.data$date),
      doy = lubridate::yday(.data$date)
    ) %>%
    tidyr::drop_na(.data$stage)
}

df_phenology_raw <- read_excel(fp, sheet = "Phenology") %>%
  rename_with(tolower) %>%
  dplyr::select(id_number, starts_with("stage"), starts_with("doy"), discard, comments) %>%
  dplyr::rename(tree_id = id_number) %>%
  dplyr::filter(!is.na(.data$tree_id))

df_phenology_transitions <- clean_phenology_transitions(df_phenology_raw)
df_phenology <- build_phenology_timeseries(df_phenology_transitions)

write_csv(df_phenology_transitions, "./data/interim/tree_phenology_transitions.csv")
write_csv(df_phenology, "./data/interim/tree_phenology.csv")

glimpse(df_phenology_transitions)
glimpse(df_phenology)


# --- Box-Level Data ------------------------------------------------------------

## Soil Isotope CN ----
df <- read_excel(fp, sheet = "Soil isotope CN") %>%
  rename(boxlabel = BoxLabel) %>%
  select(-Soil, -Robinia, -Condition) %>%
  rename_with(tolower) %>%
  mutate(across(where(is.character), ~ tolower(trimws(.x)))) %>%
  mutate(across(where(is.character), ~ gsub("_", "-", .x)))

write_csv(df, "./data/interim/box_cn_isotopes.csv")
glimpse(df)


## Soil Water ----
df <- read_excel(fp, sheet = "Soil Water") %>%
  select(-Plot, -Bloc, -box, -Drought) %>%
  rename_with(tolower) %>%
  mutate(across(where(is.character), ~ tolower(trimws(.x)))) %>%
  mutate(across(where(is.character), ~ gsub("_", "-", .x))) %>%
  pivot_longer(
    cols = matches("^\\d{2}\\.\\d{2}\\.\\d{4}$"),
    names_to = "date",
    values_to = "swc"
  ) %>%
  mutate(date = as.Date(date, format = "%d.%m.%Y")) %>%
  arrange(boxlabel, date)

write_csv(df, "./data/interim/box_soilwater.csv")
glimpse(df)


## Soil Respiration ----
df <- read_excel(fp, sheet = "Soil Respiration") %>%
  rename(boxlabel = BoxLabel) %>%
  select(-Species, -Treatment) %>%
  rename_with(tolower) %>%
  mutate(across(where(is.character), ~ tolower(trimws(.x)))) %>%
  mutate(across(where(is.character), ~ gsub("_", "-", .x))) %>%
  pivot_longer(
    cols = starts_with("co2mean_"),
    names_to = "date",
    names_prefix = "co2mean_",
    values_to = "co2"
  ) %>%
  mutate(date = as.Date(date, format = "%d.%m.%Y")) %>%
  arrange(boxlabel, date)

write_csv(df, "./data/interim/box_respiration.csv")
glimpse(df)


# --- Examples for calling data easily ------------------------------------------
# get_meta() and get_data() are helpers to load data and join meta.
# Examples:
# tree_meta <- get_meta("tree")
# tree_growth <- get_data("tree", "growth", with_meta = TRUE, path = "./data/interim")
# box_soilwater <- get_data("box", "soilwater", with_meta = TRUE, path = "./data/interim")
