#!/usr/bin/env Rscript
SUPP_SCRIPT_FILE <- normalizePath(gsub("~\\+~", " ", sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), mustWork = TRUE)
source(file.path(dirname(SUPP_SCRIPT_FILE), "_v1-figure-helpers.R"))
supp_temporal_effect_figure("fig-s3-2-vitality-effects-bootstrap", "condition", "condition", c("fagus", "quercus"), "vitality", 2026082400L)
