# Shared setup for the supplementary phenology scripts.
#
# The calling script must define SUPP_SCRIPT_FILE as its normalized path.

if (!exists("SUPP_SCRIPT_FILE", inherits = FALSE)) {
  stop("SUPP_SCRIPT_FILE must be defined before sourcing phenology_setup.R.", call. = FALSE)
}

SUPP_SCRIPT_DIR <- dirname(SUPP_SCRIPT_FILE)
PROJECT_ROOT <- normalizePath(
  file.path(SUPP_SCRIPT_DIR, "..", ".."),
  winslash = "/",
  mustWork = TRUE
)
SUPP_OUTPUT_DIR <- file.path(PROJECT_ROOT, "output", "supplementary")
dir.create(SUPP_OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
setwd(PROJECT_ROOT)

renv_lib <- Sys.glob(file.path(PROJECT_ROOT, "renv", "library", "*", "R-*", "*"))
if (length(renv_lib)) {
  .libPaths(c(normalizePath(renv_lib[[1]], winslash = "/", mustWork = TRUE), .libPaths()))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
})

suppressMessages(suppressPackageStartupMessages({
  function_dir <- file.path(PROJECT_ROOT, "scripts", "functions")
  source(file.path(function_dir, "bootstrap.R"))
  source(file.path(function_dir, "figure_data.R"))
}))

SPECIES <- c("fagus", "quercus")
SPECIES_LABELS <- c(fagus = "Fagus", quercus = "Quercus")
SPECIES_SCIENTIFIC <- c(fagus = "Fagus sylvatica", quercus = "Quercus ilex")
STAGES_INFERENCE <- 2:4
SOIL_FILTER <- "both"
INCLUDE_SOIL_TREATMENT <- FALSE

PRECIPITATION_COLORS <- c(control = "#4F6674", drought = "#D65F5F")
TREATMENT_COLORS <- c(
  precipitation = "#1B9E77",
  robinia = "#D95F02",
  culture = "#7570B3"
)

theme_supplement <- function(base_size = 8) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      text = ggplot2::element_text(color = "black"),
      axis.text = ggplot2::element_text(color = "black"),
      axis.line = ggplot2::element_line(linewidth = 0.25, color = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.25, color = "black"),
      panel.grid.major = ggplot2::element_line(color = "grey90", linewidth = 0.2),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey92", color = "grey75", linewidth = 0.25),
      strip.text = ggplot2::element_text(color = "black", face = "plain", size = base_size + 0.2),
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank(),
      plot.caption = ggplot2::element_blank(),
      legend.position = "bottom"
    )
}

save_supplementary_plot <- function(plot, stem, width_mm = 160, height_mm) {
  if (width_mm > 160) {
    stop("Supplementary figure width exceeds 160 mm.", call. = FALSE)
  }

  pdf_path <- file.path(SUPP_OUTPUT_DIR, paste0(stem, ".pdf"))
  ggplot2::ggsave(
    pdf_path,
    plot,
    width = width_mm,
    height = height_mm,
    units = "mm",
    bg = "white",
    limitsize = FALSE
  )

  invisible(pdf_path)
}
