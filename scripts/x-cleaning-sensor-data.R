# Libraries
library(tidyverse)
library(lubridate)


# Meteo Data ----------------------------------------------------------------------------------

# Load 10 minutes meteo data
meteo10 <- read.table(
  "./data/raw/sensor_data/meteo_10min.dat",
  sep = ",",
  header = FALSE,
  skip = 4,
  stringsAsFactors = FALSE,
  col.names = c(
    "datetime",
    "id",
    "air_temp",
    "rh",
    "patm",
    "vpd",
    "radiation"
  )
)

# Load 1 minutes meteo data
meteo1 <- read.table(
  "./data/raw/sensor_data/meteo_1min.dat",
  sep = ",",
  header = FALSE,
  skip = 4,
  stringsAsFactors = FALSE,
  col.names = c(
    "datetime",
    "id",
    "air_temp",
    "rh",
    "patm",
    "vpd",
    "radiation",
    "CS320_Raw",
    "CS320Temp"
  )
)

# Get daily means from 10 minutes data
meteo10_daily <- x %>%
  mutate(
    datetime = ymd_hms(datetime),
    date = as.Date(datetime)
  ) %>%
  mutate(
    across(c(air_temp, rh, vpd, radiation), as.numeric)
  ) %>%
  group_by(id, date) %>%
  summarise(
    across(
      c(air_temp, rh, vpd, radiation),
      mean,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


# Soil Data -----------------------------------------------------------------------------------

# Load 
soil10 <- read.table(
  "./data/raw/sensor_data/soil_10min.dat",
  sep = ",",
  header = TRUE,
  skip = 2,
  stringsAsFactors = FALSE) |> 
  slice(-1)

# Daily mean
soil10_daily <- soil10 %>%
  mutate(
    datetime = ymd_hms(TS),
    date = as.Date(datetime)
  ) %>%
  select(date, starts_with("kPa"), starts_with("DegC")) %>%
  mutate(
    across(
      -date,
      ~ as.numeric(na_if(.x, "NAN"))
    )
  ) %>%
  group_by(date) %>%
  summarise(
    across(
      everything(),
      mean,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

# Match sensor to boxlabel
sensor_lookup <- tibble(
  sensor_nr = 1:24,
  boxlabel = c(
    "b2-p6-c1", "b2-p6-c2", "b2-p6-c3", "b2-p6-c4",
    "b2-p5-c1", "b2-p5-c2", "b2-p5-c3", "b2-p5-c4",
    "b2-p4-c1", "b2-p4-c2", "b2-p4-c3", "b2-p4-c4",
    "b2-p3-c1", "b2-p3-c2", "b2-p3-c3", "b2-p3-c4",
    "b2-p2-c1", "b2-p2-c2", "b2-p2-c3", "b2-p2-c4",
    "b2-p1-c1", "b2-p1-c2", "b2-p1-c3", "b2-p1-c4"
  )
)

soil_mp_long <- soil10_daily %>%
  pivot_longer(
    cols = starts_with("kPa"),
    names_to = "sensor_raw",
    values_to = "soil_mp"
  ) %>%
  mutate(
    sensor_nr = if_else(
      sensor_raw == "kPa",
      1L,
      as.integer(str_remove(sensor_raw, "kPa\\.")) + 1L
    )
  ) %>%
  select(date, sensor_nr, soil_mp)

soil_temp_long <- soil10_daily %>%
  pivot_longer(
    cols = starts_with("DegC"),
    names_to = "sensor_raw",
    values_to = "soil_temp"
  ) %>%
  mutate(
    sensor_nr = if_else(
      sensor_raw == "DegC",
      1L,
      as.integer(str_remove(sensor_raw, "DegC\\.")) + 1L
    )
  ) %>%
  select(date, sensor_nr, soil_temp)

# Take daily average
soil_daily_long <- soil_mp_long %>%
  left_join(
    soil_temp_long,
    by = c("date", "sensor_nr")
  ) %>%
  left_join(
    sensor_lookup,
    by = "sensor_nr"
  )

# Attach more box info
box_meta <- get_data("box", "soilwater") |> select(boxlabel, date, swc, precipitation_soiltype_robinia, precipitation, robinia, soiltype, culture, species_mix)
soil_daily_long <- soil_daily_long |> left_join(box_meta |> select(boxlabel, date, swc))
soil_daily_long <- soil_daily_long |> left_join(box_meta |> select(-date, -swc) |> unique())


# --- Plot Soil

library(tidyverse)
library(lubridate)

set.seed(44)

random_month <- soil_daily_long %>%
  mutate(month = floor_date(date, "month")) %>%
  distinct(month) %>%
  slice_sample(n = 1) %>%
  pull(month)


# Random month
# random_month <- lubridate::date("2025-06-01")
# df_plot <- soil_daily_long %>%
#   filter(
#     date >= random_month,
#     date < random_month + months(3)
#   )

# Defined period
df_plot <- soil_daily_long %>%
  filter(
    date >= lubridate::date("2025-01-01"),
    date < lubridate::date("2025-12-31")
  )

# Plot

## Limits
ylim <- -5 # for ylimit
ybar <- ylim - ylim * 0.05 # for heatwave bar

# Heatwave dates
DROUGHT_PERIODS <- list(c("2025-06-20","2025-07-07"),
                        c("2025-08-12","2025-08-20"))

ggplot(df_plot, aes(x = date, y = -log(-soil_mp, base=10), color = precipitation)) +
  geom_smooth() +
  # geom_line() +
  facet_wrap(~ soiltype * species_mix, scales = "free_y", ncol=6) +
  ylim(ylim, 1) +
  scale_color_manual(
    values = c(
      "drought" = "indianred",
      "control" = "steelblue"
    )
  ) +
  annotate("segment",
    x = as.Date(DROUGHT_PERIODS[[1]][1]), xend = as.Date(DROUGHT_PERIODS[[1]][2]),
        y = ybar, yend = ybar, size = 2.5, color = "orange", lineend = "round", alpha = 0.9) +
  annotate("segment",
           x = as.Date(DROUGHT_PERIODS[[2]][1]), xend = as.Date(DROUGHT_PERIODS[[2]][2]),
           y = ybar, yend = ybar, size = 2.5, color = "orange", lineend = "round", alpha = 0.9)


 # New ----------------------------------------------------------------------------------------
source("./functions/0-summary_figures.R")

df_plot <- soil_daily_long %>%
  filter(
    date >= lubridate::date("2025-01-01"),
    date < lubridate::date("2025-12-31")
  )

df_plot2 <- df_plot %>%
  dplyr::filter(!is.na(soil_mp), soil_mp < 0) %>%
  dplyr::mutate(soil_mp_log = log10(-soil_mp)) |> 
  dplyr::mutate(
    robinia = factor(robinia, levels = c("without-robinia", "with-robinia"))
  )

p <- build_timeseries(
  df_plot2,
  soil_mp_log,
  facet = "soil",
  style = "band",
  ylab = expression(log[10](-psi~"(kPa)")),
  title = "Soil water potential over time"
) +
  scale_y_reverse(
    breaks = c(-3,-2,-1,0,1,2,3,4),
    labels = function(x) paste0("-", 10^x)
  ) +
  annotate("segment",
           x = as.Date(DROUGHT_PERIODS[[1]][1]), xend = as.Date(DROUGHT_PERIODS[[1]][2]),
           y = 5, yend = 5, size = 2.5, color = "orange", lineend = "round", alpha = 0.9) +
  annotate("segment",
           x = as.Date(DROUGHT_PERIODS[[2]][1]), xend = as.Date(DROUGHT_PERIODS[[2]][2]),
           y = 5, yend = 5, size = 2.5, color = "orange", lineend = "round", alpha = 0.9)

p
