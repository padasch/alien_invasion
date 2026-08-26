#!/usr/bin/env Rscript
SUPP_SCRIPT_FILE <- normalizePath(gsub("~\\+~", " ", sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), mustWork = TRUE)
source(file.path(dirname(SUPP_SCRIPT_FILE), "figure_helpers.R"))
supp_timeseries_figure("fig_s6", "chlorophyll", "chl", c("fagus", "quercus"), "Chlorophyll content (SPAD)", 2026082204L, c(0, 35))
