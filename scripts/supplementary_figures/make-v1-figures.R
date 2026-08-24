#!/usr/bin/env Rscript

arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- normalizePath(gsub("~\\+~", " ", sub("^--file=", "", arg[[1]])), mustWork = TRUE)
script_dir <- dirname(script_file)
scripts <- list.files(script_dir, pattern = "^fig_s([1-9]|1[0-5])[.]R$", full.names = TRUE)
figure_numbers <- as.integer(gsub("[^0-9]", "", basename(scripts)))
scripts <- scripts[order(figure_numbers)]

for (script in scripts) {
  message("Running ", basename(script))
  status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(script)))
  if (!identical(status, 0L)) stop("Figure script failed: ", script, call. = FALSE)
}

message("Completed ", length(scripts), " supplementary figure scripts.")
