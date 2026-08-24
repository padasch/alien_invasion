#!/usr/bin/env Rscript
SUPP_SCRIPT_FILE <- normalizePath(gsub("~\\+~", " ", sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), mustWork = TRUE)
source(file.path(dirname(SUPP_SCRIPT_FILE), "_v1-figure-helpers.R"))
supp_timeseries_figure("fig_s5", "quantum_yield", "qy", c("fagus", "quercus"), "Quantum yield (Fv/Fm)", 2026082203L, c(0, 1))
