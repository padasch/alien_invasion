# A minimal agnostic SWC imputation pipeline using GAM.
#
# Goal:
# - Build daily SWC estimates for all boxes without treatment predictors.
# - Use only measured SWC, measured SWP (sensor subset), and the shared
#   site-level daily climate product.
#
# Output:
# - data/interim/box_soilwater_daily_gam_agnostic.csv

source("./scripts/auxiliary/functions/_source.R")

library(tidyverse)
library(lubridate)
library(mgcv)

# -----------------------------------------------------------------------------
# 1) Load daily climate (site-level)
# -----------------------------------------------------------------------------
climate_daily <- get_climate()

# -----------------------------------------------------------------------------
# 2) Load and wrangle daily SWP signal (site-level aggregate from sensor subset)
# -----------------------------------------------------------------------------
soil_file <- .resolve_path("data/raw/sensor_data/soil_10min.dat")
soil_raw <- read.table(
  soil_file,
  sep = ",",
  header = TRUE,
  skip = 2,
  stringsAsFactors = FALSE
) %>%
  slice(-1)

soil_daily <- soil_raw %>%
  mutate(
    datetime = ymd_hms(TS),
    date = as.Date(datetime)
  ) %>%
  select(date, starts_with("kPa")) %>%
  mutate(
    across(-date, ~ as.numeric(na_if(.x, "NAN")))
  ) %>%
  group_by(date) %>%
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE)), .groups = "drop")

swp_site_daily <- soil_daily %>%
  pivot_longer(
    cols = starts_with("kPa"),
    names_to = "sensor_raw",
    values_to = "soil_mp"
  ) %>%
  filter(!is.na(soil_mp), soil_mp < 0) %>%
  mutate(swp_log10 = log10(-soil_mp)) %>%
  group_by(date) %>%
  summarise(
    swp_site = median(swp_log10, na.rm = TRUE),
    swp_site_sd = sd(swp_log10, na.rm = TRUE),
    .groups = "drop"
  )

# -----------------------------------------------------------------------------
# 3) Load observed SWC and create full box x date panel
# -----------------------------------------------------------------------------
swc_obs <- readr::read_csv(
  .resolve_path("data/interim/box_soilwater.csv"),
  show_col_types = FALSE
) %>%
  transmute(
    boxlabel = as.character(boxlabel),
    date = as.Date(date),
    swc_obs = as.numeric(swc)
)

# Keep only the date window where all site-level predictors are available.
date_min <- max(min(climate_daily$date), min(swp_site_daily$date))
date_max <- min(max(climate_daily$date), max(swp_site_daily$date))

all_boxes <- swc_obs %>% distinct(boxlabel)
all_dates <- tibble(date = seq.Date(date_min, date_max, by = "day"))

panel_daily <- tidyr::expand_grid(all_boxes, all_dates) %>%
  left_join(swc_obs, by = c("boxlabel", "date")) %>%
  left_join(climate_daily, by = "date") %>%
  left_join(swp_site_daily %>% select(date, swp_site), by = "date") %>%
  mutate(
    date_num = as.numeric(date),
    boxlabel = factor(boxlabel)
  )

fit_data <- panel_daily %>%
  filter(!is.na(swc_obs), !is.na(swp_site), !is.na(air_temp), !is.na(vpd),
         !is.na(radiation), !is.na(precip_mm))

cat("Fitting rows:", nrow(fit_data), "\n")
cat("Boxes in fit:", n_distinct(fit_data$boxlabel), "\n")
cat("Date range in fit:", as.character(min(fit_data$date)), "to", as.character(max(fit_data$date)), "\n")

if (nrow(fit_data) < 200) {
  stop("Too few rows for stable GAM fit.", call. = FALSE)
}

# -----------------------------------------------------------------------------
# 4) Fit agnostic GAM (no treatment covariates)
# -----------------------------------------------------------------------------
gam_swc <- mgcv::gam(
  swc_obs ~
    s(date_num, k = 10) +
    s(swp_site, k = 8) +
    s(air_temp, k = 8) +
    s(vpd, k = 8) +
    s(precip_mm, k = 8) +
    s(radiation, k = 8) +
    s(boxlabel, bs = "re"),
  data = fit_data,
  method = "REML"
)

cat("\nGAM fitted. Adjusted R2:", round(summary(gam_swc)$r.sq, 3), "\n")

# -----------------------------------------------------------------------------
# 5) Predict daily SWC for all boxes and export
# -----------------------------------------------------------------------------
pred_data <- panel_daily %>%
  filter(!is.na(swp_site), !is.na(air_temp), !is.na(vpd), !is.na(radiation), !is.na(precip_mm)) %>%
  select(
    boxlabel, date, swc_obs, air_temp, rh, vpd, radiation,
    precip_mm, swp_site, date_num
  )

pred <- predict(gam_swc, newdata = pred_data, se.fit = TRUE, type = "response")

out <- pred_data %>%
  mutate(
    swc_hat = as.numeric(pred$fit),
    swc_hat_se = as.numeric(pred$se.fit),
    swc_hat_lo = swc_hat - 1.96 * swc_hat_se,
    swc_hat_hi = swc_hat + 1.96 * swc_hat_se,
    # Keep SWC in physically meaningful bounds.
    swc_hat = pmin(pmax(swc_hat, 0), 100),
    swc_hat_lo = pmin(pmax(swc_hat_lo, 0), 100),
    swc_hat_hi = pmin(pmax(swc_hat_hi, 0), 100),
    swc_source = if_else(!is.na(swc_obs), "measured", "imputed_gam")
  ) %>%
  arrange(boxlabel, date)

out_file <- .resolve_path("data/interim/box_soilwater_daily_gam_agnostic.csv")
readr::write_csv(out, out_file)

cat("\nSaved:", out_file, "\n")
cat("Rows:", nrow(out), "| Measured:", sum(out$swc_source == "measured"), "| Imputed:", sum(out$swc_source == "imputed_gam"), "\n")
