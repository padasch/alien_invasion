# 0-data_cleaning.R

# --- Setup --------------------------------------------------------------------
library(readxl)
library(tidyr)
library(readr)
library(janitor)
library(dplyr)
library(purrr)

# Load helpers (if any)
source("./functions/_source.R")

# Set filepath
fp <- "./data/raw/data_20251203.xlsx"

# Print all sheet names
print(excel_sheets(fp))


# --- Meta Data -----------------------------------------------------------------

# Tree meta data
meta_tree <- read_excel(fp, sheet = "All_labels") %>%
  remove_empty(which = "cols") %>%
  clean_names() %>%
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
  remove_empty(which = "cols") %>%
  clean_names() %>%
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
  clean_names() %>%
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
  clean_names() %>%
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
    condition = as.double(leafstate),
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
    comment = as.character(Comment),
    chlavg = rowMeans(cbind(chl1, chl2), na.rm = TRUE)
  ) %>%
  select(tree_id, date, percent_senesced, chl1, chl2, chlavg, comment) %>%
  arrange(tree_id, date) |>
  filter(!is.na(tree_id))

write_csv(df, "./data/interim/tree_senescence.csv")
glimpse(df)


## Growth (Height and Diameter) ----
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
  dplyr::select(tree_id, date, diameter, height) %>%
  arrange(tree_id, date) |>
  filter(!is.na(tree_id)) |> 
  group_by(tree_id) |>
  mutate(
    # first measurements per tree (baseline)
    first_diameter = first(diameter),
    first_height   = first(height),
    
    # cumulative absolute increments from first measurement
    diameter_inc_t0 = diameter - first_diameter,
    height_inc_t0   = height   - first_height,
    
    # relative size and relative cumulative increment
    diameter_rel        = diameter / first_diameter,
    height_rel          = height   / first_height,
    diameter_inc_t0_rel = diameter_inc_t0 / first_diameter,
    height_inc_t0_rel   = height_inc_t0   / first_height,
    
    # time between consecutive measurements [years]
    delta_t_years = as.numeric(date - dplyr::lag(date)) / 365.25,
    
    # increments between consecutive measurements
    diameter_inc_dt = diameter - dplyr::lag(diameter),
    height_inc_dt   = height   - dplyr::lag(height),
    
    # absolute growth rate (AGR) between dates
    diameter_agr = diameter_inc_dt / delta_t_years,
    height_agr   = height_inc_dt   / delta_t_years,
    
    # relative growth rate (RGR) between dates (log-scale)
    diameter_rgr = (log(diameter) - log(dplyr::lag(diameter))) / delta_t_years,
    height_rgr   = (log(height)   - log(dplyr::lag(height)))   / delta_t_years
  ) %>%
  ungroup()

write_csv(df, "./data/interim/tree_growth.csv")
glimpse(df)


## Specific Leaf Area ----
df <- read_excel(fp, sheet = "SLA") %>%
  select(-Label) %>%
  rename_with(tolower) |>
  rename(tree_id = id_number) |>
  filter(!is.na(tree_id))

write_csv(df, "./data/interim/tree_sla.csv")
glimpse(df)


## Phenology ----
build_phenology_timeseries <- function(df) {
  # 0) Derive measurement dates from all stage columns
  meas_dates <- df %>%
    select(starts_with("stage")) %>%
    mutate(across(everything(), as.Date)) %>%
    unlist(use.names = FALSE)
  
  meas_dates <- meas_dates[!is.na(meas_dates)]
  meas_dates <- sort(unique(as.Date(meas_dates)))
  
  # 1) Reshape stage1...stage4 to long format
  df_long <- df %>%
    filter(discard == 0) %>%                    # drop discarded trees if needed
    pivot_longer(
      cols      = starts_with("stage"),
      names_to  = "stage_name",
      values_to = "stage_date"
    ) %>%
    filter(!is.na(stage_date)) %>%              # keep only observed transitions
    mutate(
      stage_date = as.Date(stage_date),
      stage      = parse_number(stage_name)     # stage1 -> 1, stage2 -> 2, ...
    ) %>%
    arrange(tree_id, stage_date) %>%
    select(tree_id, stage, stage_date)
  
  # 2) Build time series per tree using a vectorized findInterval
  trees <- df_long %>%
    pull(tree_id) %>%
    unique()
  
  res <- map_dfr(trees, function(id) {
    stage_tree <- df_long %>%
      filter(tree_id == id) %>%
      arrange(stage_date)
    
    # If a tree has no stages at all (just in case)
    if (nrow(stage_tree) == 0) {
      return(
        tibble(
          tree_id = id,
          date    = meas_dates,
          stage   = NA_integer_
        )
      )
    }
    
    # For each meas_date: last stage_date <= meas_date
    idx <- findInterval(meas_dates, stage_tree$stage_date)
    
    stage_vec <- ifelse(
      idx == 0,
      NA_integer_,                # before first observed stage -> NA
      stage_tree$stage[idx]
    )
    
    tibble(
      tree_id = id,
      date    = meas_dates,
      stage   = stage_vec
    )
  })
  
  # Optional: make stage an ordered factor and add DOY
  res %>%
    mutate(
      tree_id = as.character(tree_id),
      stage = factor(stage,
                     levels = sort(unique(na.omit(stage))),
                     ordered = TRUE
      ),
      doy = yday(date)
    )
}

df <- read_excel(fp, sheet = "Phenology") %>%
  rename_with(tolower) %>%
  dplyr::select(id_number, starts_with("stage"), starts_with("doy"), discard, comments, ) |>
  rename(tree_id = id_number) |>
  filter(!is.na(tree_id)) |> 
  build_phenology_timeseries() |> 
  drop_na(stage)

write_csv(df, "./data/interim/tree_phenology.csv")
glimpse(df)


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

