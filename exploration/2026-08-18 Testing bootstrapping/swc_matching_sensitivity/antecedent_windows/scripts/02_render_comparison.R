#!/usr/bin/env Rscript

# Re-render the compact diagnostic from saved bootstrap summaries without
# re-running the 48,000 mixed-model fits.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- gsub("~[+]~", " ", sub("^--file=", "", script_arg[[1]]))
script_file <- normalizePath(script_path, winslash = "/", mustWork = TRUE)
analysis_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
output_dir <- file.path(analysis_dir, "output")

plot_data <- read_csv(
  file.path(output_dir, "all-method-bootstrap-effects.csv"),
  show_col_types = FALSE
) %>%
  filter(component %in% c("direct", "indirect", "total")) %>%
  mutate(
    method = factor(
      method,
      levels = c("fuzzy_baseline", "antecedent_7d", "antecedent_14d"),
      labels = c("Fuzzy measured SWC", "Modeled prior 7-day mean", "Modeled prior 14-day mean")
    ),
    treatment = factor(
      treatment,
      levels = rev(c("precipitation", "robinia", "culture", "extreme_event"))
    )
  )

pdf(file.path(output_dir, "antecedent-window-sem-effect-comparison.pdf"),
    width = 11, height = 10)
for (sp in c("fagus", "quercus")) {
  p <- plot_data %>%
    filter(species == sp) %>%
    ggplot(aes(estimate, treatment, colour = method)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
    geom_errorbarh(
      aes(xmin = lower, xmax = upper),
      position = position_dodge(width = 0.55), height = 0
    ) +
    geom_point(position = position_dodge(width = 0.55), size = 1.5) +
    facet_grid(response_label ~ component, scales = "free_x") +
    labs(
      title = paste("Antecedent-SWC SEM sensitivity:", tools::toTitleCase(sp)),
      x = "Oriented standardized effect (95% container-bootstrap CI)", y = NULL,
      colour = "Mediator definition"
    ) +
    theme_classic(base_size = 8) +
    theme(
      legend.position = "bottom", strip.text = element_text(size = 7),
      plot.title = element_text(hjust = 0.5)
    )
  print(p)
}
dev.off()
