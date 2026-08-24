#!/usr/bin/env Rscript
SUPP_SCRIPT_FILE <- normalizePath(gsub("~\\+~", " ", sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), mustWork = TRUE)
source(file.path(dirname(SUPP_SCRIPT_FILE), "_v1-figure-helpers.R"))
supp_timeseries_figure("fig-s2-2-height-timeseries", "growth", "height_inc_t0", c("fagus", "quercus", "robinia"), "Height increment (cm)", 2026082202L, c(0, 60))
