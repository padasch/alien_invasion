#!/usr/bin/env Rscript

# Build the shared site-level daily climate product.
#
# Daily means are calculated from every climate variable in Meteo_10min.
# Precipitation is summed from 10-minute observations, using the existing
# MeteoSwiss series with the required LWF replacement on 2025-07-13--15.

source("./functions/_source.R")

.daily_mean <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

meteo_file <- .resolve_path("data/raw/sensor_data/meteo_10min.dat")
climate_file <- .resolve_path("data/interim/site_climate_daily.csv")

meteo_daily <- readr::read_csv(
  meteo_file,
  skip = 4,
  col_names = c(
    "timestamp", "record", "air_temp", "rh", "patm", "vpd", "radiation"
  ),
  col_types = readr::cols(
    timestamp = readr::col_character(),
    record = readr::col_double(),
    air_temp = readr::col_double(),
    rh = readr::col_double(),
    patm = readr::col_double(),
    vpd = readr::col_double(),
    radiation = readr::col_double()
  ),
  na = c("", "NA", "NAN"),
  show_col_types = FALSE
) %>%
  dplyr::mutate(
    timestamp = lubridate::ymd_hms(.data$timestamp, tz = "UTC", quiet = TRUE),
    date = as.Date(.data$timestamp)
  ) %>%
  dplyr::filter(!is.na(.data$date)) %>%
  dplyr::group_by(.data$date) %>%
  dplyr::summarise(
    air_temp = .daily_mean(.data$air_temp),
    rh = .daily_mean(.data$rh),
    patm = .daily_mean(.data$patm),
    vpd = .daily_mean(.data$vpd),
    radiation = .daily_mean(.data$radiation),
    .groups = "drop"
  )

precipitation_daily <- get_site_precipitation_daily()

climate_daily <- dplyr::full_join(
  meteo_daily,
  precipitation_daily,
  by = "date"
) %>%
  dplyr::arrange(.data$date) %>%
  dplyr::select(dplyr::all_of(c(
    "date", "air_temp", "rh", "patm", "vpd", "radiation", "precip_mm"
  )))

readr::write_csv(climate_daily, climate_file)

message(
  "Saved ", climate_file, " (", nrow(climate_daily), " days; ",
  as.character(min(climate_daily$date)), " to ",
  as.character(max(climate_daily$date)), ")."
)
