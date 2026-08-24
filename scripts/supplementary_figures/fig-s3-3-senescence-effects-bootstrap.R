#!/usr/bin/env Rscript
SUPP_SCRIPT_FILE <- normalizePath(gsub("~\\+~", " ", sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), mustWork = TRUE)
source(file.path(dirname(SUPP_SCRIPT_FILE), "_v1-figure-helpers.R"))
supp_temporal_effect_figure("fig-s3-3-senescence-effects-bootstrap", "senescence", "remaining_green", c("fagus"), "remaining green crown", 2026082500L)
