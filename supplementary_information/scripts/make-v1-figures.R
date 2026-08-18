#!/usr/bin/env Rscript

arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- normalizePath(gsub("~\\+~", " ", sub("^--file=", "", arg[[1]])), mustWork = TRUE)
script_dir <- dirname(script_file)
scripts <- list.files(script_dir, pattern = "^fig-s[1-5]-.*\\.R$", full.names = TRUE)

for (script in sort(scripts)) {
  message("Running ", basename(script))
  status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(script)))
  if (!identical(status, 0L)) stop("Figure script failed: ", script, call. = FALSE)
}

message("Completed ", length(scripts), " supplementary figure scripts.")
