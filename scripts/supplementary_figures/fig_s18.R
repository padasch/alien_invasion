#!/usr/bin/env Rscript

# Date-specific treatment effects on phase-wise relative volume increments.
# Reuse the Figure 6 response; exclude initial baseline zeros and never add
# the artificial phase-restart rows used only for drawing that time series.
# Plot contract: two species panels, treatment contrasts with bootstrap bands,
# existing supplementary styling, and lines broken where phase baselines reset.
SUPP_SCRIPT_FILE <- normalizePath(gsub("~\\+~", " ", sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), mustWork = TRUE)
source(file.path(dirname(SUPP_SCRIPT_FILE), "figure_helpers.R"))
supp_temporal_effect_figure(
  "fig_s18", "growth", "volume_inc_phase_rel", c("fagus", "quercus"),
  "phase-wise relative volume increment", 2026082700L
)
