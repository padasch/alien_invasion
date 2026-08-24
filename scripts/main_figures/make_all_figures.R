#!/usr/bin/env Rscript

script_file <- {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) {
    gsub("~[+]~", " ", sub("^--file=", "", file_arg[[1]]))
  } else {
    "scripts/main_figures/make_all_figures.R"
  }
}

script_dir <- dirname(normalizePath(script_file, winslash = "/", mustWork = TRUE))

source(file.path(script_dir, "00_setup.R"))
source(file.path(script_dir, "02_fig2_variation_timeseries.R"))
source(file.path(script_dir, "03_fig3_diameter_increment_effects.R"))
source(file.path(script_dir, "04_fig4_quantum_yield_effects.R"))
source(file.path(script_dir, "05_fig5_biomass_effects.R"))
source(file.path(script_dir, "06_fig6_volume_sem.R"))

message("Writing final figures to: ", ALINV_FINAL_OUTPUT_DIR)

outputs <- c(
  fig2 = make_fig2(),
  fig3 = make_fig3(),
  fig4 = make_fig4(),
  fig5 = make_fig5(),
  fig6 = make_fig6()
)

message("Created PDFs:")
for (path in outputs) {
  message(" - ", path)
}
