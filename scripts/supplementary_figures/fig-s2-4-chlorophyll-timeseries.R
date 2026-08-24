#!/usr/bin/env Rscript
SUPP_SCRIPT_FILE <- normalizePath(gsub("~\\+~", " ", sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), mustWork = TRUE)
source(file.path(dirname(SUPP_SCRIPT_FILE), "_v1-figure-helpers.R"))
supp_timeseries_figure("fig-s2-4-chlorophyll-timeseries", "chlorophyll", "chl", c("fagus", "quercus"), "Chlorophyll content (SPAD)", 2026082204L, c(0, 35))
