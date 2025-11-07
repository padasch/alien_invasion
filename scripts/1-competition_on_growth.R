# Setup ---------------------------------------------------------------------------------------

source("./functions/_source.R")


# Phenology -----------------------------------------------------------------------------------


df_pheno <- get_data("tree", "phenology")
df_pheno |> glimpse()

df_long <- df_pheno %>%
  select(species, robinia, precipitation, culture,
         doy_s1, doy_s2, doy_s3, doy_s4) %>%
  pivot_longer(
    cols = starts_with("doy_"),
    names_to = "stage",
    values_to = "doy"
  ) %>%
  filter(!is.na(doy)) %>%
  mutate(
    stage = factor(stage,
                   levels = c("doy_s1","doy_s2","doy_s3","doy_s4"),
                   labels = c("Stage 1","Stage 2","Stage 3","Stage 4")),
    precipitation = factor(precipitation, levels = c("control","drought")),
    culture       = factor(culture, levels = c("mixed","mono")),
    robinia       = factor(robinia, levels = c("without-robinia","with-robinia"))
  )

# summarise to mean & SD for each treatment
df_summary <- df_long %>%
  group_by(species, robinia, precipitation, culture, stage) %>%
  summarise(
    mean_doy = mean(doy, na.rm = TRUE),
    sd_doy   = sd(doy, na.rm = TRUE),
    .groups  = "drop"
  )

ggplot(df_summary,
       aes(x = stage, y = mean_doy,
           color = precipitation,
           linetype = culture,
           group = interaction(precipitation, culture))) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_doy - sd_doy, ymax = mean_doy + sd_doy),
                width = 0.1, alpha = 0.6) +
  facet_grid(
    rows = vars(species),
    cols = vars(robinia),
    labeller = labeller(
      robinia = c(`without-robinia` = "without robinia",
                  `with-robinia`    = "with robinia"))
  ) +
  scale_color_manual(values = c(control = "black", drought = "#D62728"),
                     name = "Precipitation") +
  scale_linetype_manual(values = c(mixed = "solid", mono = "dashed"),
                        name = "Culture") +
  labs(
    x = "Phenological stage",
    y = "Mean DOY",
    title = "Progression of phenological stages by species, robinia, precipitation, and culture"
  ) + 
  coord_cartesian(ylim=c(85, 125))


# Specific Leaf Area --------------------------------------------------------------------------

df_sla <- get_data("tree", "sla")
df_sla |> glimpse()

library(ggplot2)
library(dplyr)

df_sla %>%
  mutate(
    precipitation = factor(precipitation, levels = c("control", "drought")),
    culture       = factor(culture, levels = c("mixed", "mono")),
    robinia       = factor(robinia, levels = c("without-robinia", "with-robinia")),
    species       = factor(species)
  ) %>%
  ggplot(
    aes(x = culture, y = sla_mm2_per_mg,
        fill = precipitation, color = precipitation,
        alpha = culture)
  ) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.65,
    outlier.shape = NA
  ) +
  facet_grid(
    rows = vars(species),
    cols = vars(robinia),
    labeller = labeller(
      robinia = c(`without-robinia` = "without robinia",
                  `with-robinia`    = "with robinia"))
  ) +
  scale_fill_manual(values = c(control = "grey70", drought = "#E64B35"),
                    name = "Precipitation") +
  scale_color_manual(values = c(control = "black", drought = "#D62728"),
                     name = "Precipitation") +
  scale_alpha_manual(values = c(mixed = 1, mono = 0.6),
                     name = "Culture") +
  labs(
    x = "Culture", y = expression(SLA~(mm^2~mg^-1)),
    title = "Specific Leaf Area (SLA) across species, precipitation, culture, and robinia"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    axis.line.x.bottom = element_line(),
    axis.line.y.left   = element_line()
  )

# Growth ------------------------------------------------------------------------------------

df_growth <- get_data("tree", "growth")
df_growth |> glimpse()

# Choose soiltype filter ("inoc_beech", "inoc_robinia", or "both")
# filter_soiltype <- "inoc-robinia"
filter_soiltype <- "inoc-beech"
# filter_soiltype <- "both"

## Timeseries ----
### Diameter ----
df_plot <- df_growth %>%
  {
    if (filter_soiltype != "both") {
      filter(., soiltype == filter_soiltype)
    } else {
      .
    }
  } %>%
  filter(!is.na(diameter)) %>%
  group_by(tree_id) %>%
  arrange(date, .by_group = TRUE) %>%
  mutate(diam_inc = diameter - first(diameter)) %>%
  ungroup() %>%
  group_by(species, robinia, precipitation, culture, date) %>%
  summarize(
    mean_inc = mean(diam_inc, na.rm = TRUE),
    sd_inc = sd(diam_inc, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    ymin = mean_inc - sd_inc,
    ymax = mean_inc + sd_inc,
    precipitation = factor(precipitation, levels = c("control", "drought")),
    robinia = factor(robinia, levels = c("without-robinia", "with-robinia"))
  )

ggplot(df_plot, aes(
  x = date, y = mean_inc,
  color = precipitation, fill = precipitation,
  linetype = culture,
  group = interaction(precipitation, culture)
)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.15, color = NA) +
  geom_line(size = 1) +
  facet_grid(
    rows = vars(species),
    cols = vars(robinia),
    scales = "fixed"
  ) +
  scale_color_manual(values = c(control = "black", drought = "#D62728"), name = "Precipitation") +
  scale_fill_manual(values = c(control = "black", drought = "#D62728"), name = "Precipitation") +
  scale_x_date(date_labels = "%m/%y") +
  coord_cartesian(ylim = c(0, 8)) + # Diameter
  labs(
    title = paste(
      "Diameter growth increment over time",
      if (filter_soiltype != "both") paste0("(soiltype = ", filter_soiltype, ")") else ""
    ),
    x = "Date", y = "Diameter growth increment (mm)",
    linetype = "Culture"
  )

### Height ----
df_plot <- df_growth %>%
  {
    if (filter_soiltype != "both") {
      filter(., soiltype == filter_soiltype)
    } else {
      .
    }
  } %>%
  filter(!is.na(height)) %>%
  group_by(tree_id) %>%
  arrange(date, .by_group = TRUE) %>%
  mutate(inc = height - first(height)) %>%
  ungroup() %>%
  group_by(species, robinia, precipitation, culture, date) %>%
  summarize(
    mean_inc = mean(inc, na.rm = TRUE),
    sd_inc = sd(inc, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    ymin = mean_inc - sd_inc,
    ymax = mean_inc + sd_inc,
    precipitation = factor(precipitation, levels = c("control", "drought")),
    robinia = factor(robinia, levels = c("without-robinia", "with-robinia"))
  )

ggplot(df_plot, aes(
  x = date, y = mean_inc,
  color = precipitation, fill = precipitation,
  linetype = culture,
  group = interaction(precipitation, culture)
)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.15, color = NA) +
  geom_line(size = 1) +
  facet_grid(
    rows = vars(species),
    cols = vars(robinia),
    scales = "fixed"
  ) +
  scale_color_manual(values = c(control = "black", drought = "#D62728"), name = "Precipitation") +
  scale_fill_manual(values = c(control = "black", drought = "#D62728"), name = "Precipitation") +
  scale_x_date(date_labels = "%m/%y") +
  coord_cartesian(ylim = c(0, 75)) + # Height
  labs(
    title = paste(
      "Height growth increment over time",
      if (filter_soiltype != "both") paste0("(soiltype = ", filter_soiltype, ")") else ""
    ),
    x = "Date", y = "Height growth increment (cm)",
    linetype = "Culture"
  )

## Pre-Extreme Event ----
### Diameter ----
df_growth %>%
  filter(!is.na(diameter)) %>%
  group_by(tree_id) %>%
  arrange(date, .by_group = TRUE) %>%
  mutate(diam_inc = diameter - first(diameter)) %>%
  ungroup() %>%
  filter(month(date) == 7, !is.na(diam_inc)) %>% # only July measurements
  mutate(
    precipitation = factor(precipitation, levels = c("control", "drought")),
    culture = factor(culture, levels = c("mixed", "mono")),
    robinia = factor(robinia, levels = c("without-robinia", "with-robinia"))
  ) %>%
  ggplot(aes(
    x = culture, y = diam_inc,
    fill = precipitation, color = precipitation,
    alpha = culture
  )) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.65, outlier.shape = NA) +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
    size = 1.2, alpha = 0.6, show.legend = FALSE
  ) +
  facet_grid(
    rows = vars(species),
    cols = vars(robinia),
    labeller = labeller(robinia = c(
      `without-robinia` = "without robinia",
      `with-robinia` = "with robinia"
    ))
  ) +
  scale_color_manual(values = c(control = "black", drought = "#D62728"), name = "Precipitation") +
  scale_fill_manual(values = c(control = "grey70", drought = "#E64B35"), name = "Precipitation") +
  scale_alpha_manual(values = c(mixed = 1, mono = 0.6), name = "Culture") +
  labs(
    x = "Culture", y = "Diameter increment (mm)",
    title = "July diameter increments by precipitation & culture",
    subtitle = "Facets: rows = species, columns = robinia"
  )

### Height ----
df_growth %>%
  {
    if (exists("filter_soiltype") && filter_soiltype != "both") {
      dplyr::filter(., soiltype == filter_soiltype)
    } else {
      .
    }
  } %>%
  filter(!is.na(height)) %>%
  group_by(tree_id) %>%
  arrange(date, .by_group = TRUE) %>%
  mutate(inc = height - first(height)) %>%
  ungroup() %>%
  filter(month(date) == 9, !is.na(inc)) %>%        # only July
  mutate(
    precipitation = factor(precipitation, levels = c("control","drought")),
    culture = factor(culture, levels = c("mono","mixed")),
    robinia = factor(robinia, levels = c("without-robinia","with-robinia"))
  ) %>%
  ggplot(aes(
    x = culture, y = inc,
    fill = precipitation,
    color = precipitation,
    alpha = culture
  )) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.65, outlier.shape = NA) +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
    size = 1.2, alpha = 0.6, show.legend = FALSE
  ) +
  facet_grid(
    rows = vars(species),
    cols = vars(robinia),
    labeller = labeller(robinia = c(
      `without-robinia` = "without robinia",
      `with-robinia` = "with robinia"
    ))
  ) +
  scale_color_manual(values = c(control = "black", drought = "#D62728"), name = "Precipitation") +
  scale_fill_manual(values = c(control = "grey70", drought = "#E64B35"), name = "Precipitation") +
  scale_alpha_manual(values = c(mixed = 1, mono = 0.1), name = "Culture") +
  labs(
    x = "Culture", y = "Diameter increment (mm)",
    title = "July height increments by precipitation & culture",
    subtitle = "Facets: rows = species, columns = robinia"
  ) +
  coord_cartesian(ylim = c(0, 75)) # Height


# Quantum Yield -------------------------------------------------------------------------------

# Load quantum yield data
df_qy <- get_data("tree", "quantym_yield")

# Choose soiltype filter
filter_soiltype <- "inoc-robinia"
# filter_soiltype <- "inoc-beech"
# filter_soiltype <- "both"

df_plot <- df_qy %>%
  {
    if (filter_soiltype != "both") {
      filter(., soiltype == filter_soiltype)
    } else {
      .
    }
  } %>%
  filter(!is.na(qy)) %>%
  group_by(species, robinia, precipitation, culture, date) %>%
  summarize(
    mean_qy = mean(qy, na.rm = TRUE),
    sd_qy   = sd(qy, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    ymin = mean_qy - sd_qy,
    ymax = mean_qy + sd_qy,
    precipitation = factor(precipitation, levels = c("control", "drought")),
    robinia = factor(robinia, levels = c("without-robinia", "with-robinia"))
  )

ggplot(df_plot, aes(
  x = date, y = mean_qy,
  color = precipitation, fill = precipitation,
  linetype = culture,
  group = interaction(precipitation, culture)
)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.1, color = NA) +
  # geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 5, alpha = 0.1) +
  geom_line(size = 1) +
  facet_grid(
    rows = vars(species),
    cols = vars(robinia),
    scales = "fixed"
  ) +
  scale_color_manual(values = c(control = "black", drought = "#D62728"), name = "Precipitation") +
  scale_fill_manual(values = c(control = "black", drought = "#D62728"), name = "Precipitation") +
  scale_x_date(date_labels = "%m/%y") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = paste(
      "Quantum yield (QY) over time",
      if (filter_soiltype != "both") paste0("(soiltype = ", filter_soiltype, ")") else ""
    ),
    x = "Date", y = "Quantum yield (QY)",
    linetype = "Culture"
  )


# Senescence ----------------------------------------------------------------------------------

df_sn <- get_data("tree", "senescence")

# Choose soiltype filter
filter_soiltype <- "inoc-beech"
# filter_soiltype <- "inoc-beech"
# filter_soiltype <- "both"

## Percentage
df_plot <- df_sn %>%
  {
    if (filter_soiltype != "both") {
      filter(., soiltype == filter_soiltype)
    } else {
      .
    }
  } %>%
  filter(!is.na(percent_senesced)) %>%
  group_by(species, robinia, precipitation, culture, date) %>%
  summarize(
    mean_sn = mean(percent_senesced),
    sd_sn   = sd(percent_senesced, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    ymin = mean_sn - sd_sn,
    ymax = mean_sn + sd_sn,
    precipitation = factor(precipitation, levels = c("control", "drought")),
    robinia = factor(robinia, levels = c("without-robinia", "with-robinia"))
  )

ggplot(df_plot, aes(
  x = date, y = mean_sn,
  color = precipitation, fill = precipitation,
  linetype = culture,
  group = interaction(precipitation, culture)
)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.1, color = NA) +
  # geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 5, alpha = 0.1) +
  geom_line(size = 1) +
  facet_grid(
    rows = vars(species),
    cols = vars(robinia),
    scales = "fixed"
  ) +
  scale_color_manual(values = c(control = "black", drought = "#D62728"), name = "Precipitation") +
  scale_fill_manual(values = c(control = "black", drought = "#D62728"), name = "Precipitation") +
  scale_x_date(date_labels = "%m/%y") +
  coord_cartesian(ylim = c(0, 100)) +
  labs(
    title = paste(
      "Senescence (%) over time",
      if (filter_soiltype != "both") paste0("(soiltype = ", filter_soiltype, ")") else ""
    ),
    x = "Date", y = "Senescence (%)",
    linetype = "Culture"
  )

## Chlorophyll
df_plot <- df_sn %>%
  {
    if (filter_soiltype != "both") {
      filter(., soiltype == filter_soiltype)
    } else {
      .
    }
  } %>%
  filter(!is.na(chlavg)) %>%
  group_by(species, robinia, precipitation, culture, date) %>%
  summarize(
    mean_sn = mean(chlavg),
    sd_sn   = sd(chlavg, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    ymin = mean_sn - sd_sn,
    ymax = mean_sn + sd_sn,
    precipitation = factor(precipitation, levels = c("control", "drought")),
    robinia = factor(robinia, levels = c("without-robinia", "with-robinia"))
  )

ggplot(df_plot, aes(
  x = date, y = mean_sn,
  color = precipitation, fill = precipitation,
  linetype = culture,
  group = interaction(precipitation, culture)
)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.1, color = NA) +
  # geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 5, alpha = 0.1) +
  geom_line(size = 1) +
  facet_grid(
    rows = vars(species),
    cols = vars(robinia),
    scales = "fixed"
  ) +
  scale_color_manual(values = c(control = "black", drought = "#D62728"), name = "Precipitation") +
  scale_fill_manual(values = c(control = "black", drought = "#D62728"), name = "Precipitation") +
  scale_x_date(date_labels = "%m/%y") +
  coord_cartesian(ylim = c(0, 22)) +
  labs(
    title = paste(
      "Chlorophyll (-) over time",
      if (filter_soiltype != "both") paste0("(soiltype = ", filter_soiltype, ")") else ""
    ),
    x = "Date", y = "Chlorophyll (-)",
    linetype = "Culture"
  )


# Tree Condition ------------------------------------------------------------------------------
df_tc <- get_data("tree", "condition")

# Choose soiltype filter
filter_soiltype <- "inoc-beech"
# filter_soiltype <- "inoc-beech"
# filter_soiltype <- "both"

## Percentage
df_plot <- df_tc %>%
  {
    if (filter_soiltype != "both") {
      filter(., soiltype == filter_soiltype)
    } else {
      .
    }
  } %>%
  filter(!is.na(condition)) %>%
  group_by(species, robinia, precipitation, culture, date) %>%
  summarize(
    mean_sn = mean(condition),
    sd_sn   = sd(condition, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    ymin = mean_sn - sd_sn,
    ymax = mean_sn + sd_sn,
    precipitation = factor(precipitation, levels = c("control", "drought")),
    robinia = factor(robinia, levels = c("without-robinia", "with-robinia"))
  )

ggplot(df_plot, aes(
  x = date, y = mean_sn,
  color = precipitation, fill = precipitation,
  linetype = culture,
  group = interaction(precipitation, culture)
)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.1, color = NA) +
  # geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 5, alpha = 0.1) +
  geom_line(size = 1) +
  facet_grid(
    rows = vars(species),
    cols = vars(robinia),
    scales = "fixed"
  ) +
  scale_color_manual(values = c(control = "black", drought = "#D62728"), name = "Precipitation") +
  scale_fill_manual(values = c(control = "black", drought = "#D62728"), name = "Precipitation") +
  scale_x_date(date_labels = "%m/%y") +
  coord_cartesian(ylim = c(0, 5)) +
  labs(
    title = paste(
      "Tree Stress (-) over time",
      if (filter_soiltype != "both") paste0("(soiltype = ", filter_soiltype, ")") else ""
    ),
    x = "Date", y = "Tree Stress (-)",
    linetype = "Culture"
  )


# Soil Water Content --------------------------------------------------------------------------

df <- get_data("box", "soilwater")
df |> glimpse()

df %>%
  mutate(
    soiltype = gsub("-", "_", soiltype),
    precipitation = factor(precipitation, levels = c("control", "drought")),
    robinia = factor(robinia, levels = c("without-robinia", "with-robinia")),
    culture = factor(culture, levels = c("mixed", "mono"))
  ) %>%
  group_by(robinia, culture, precipitation, soiltype, date) %>%
  summarize(
    mean_swc = mean(swc, na.rm = TRUE),
    sd_swc = sd(swc, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(ymin = mean_swc - sd_swc, ymax = mean_swc + sd_swc) %>%
  ggplot(aes(
    x = date, y = mean_swc,
    color = precipitation, fill = precipitation,
    linetype = soiltype,
    group = interaction(precipitation, soiltype)
  )) +
  annotate("segment",
    x = as.Date("2025-06-20"), xend = as.Date("2025-07-02"),
    y = -1, yend = -1,
    size = 2.5, color = "orange", lineend = "round"
  ) +
  annotate("segment",
    x = as.Date("2025-08-12"), xend = as.Date("2025-08-20"),
    y = -1, yend = -1,
    size = 2.5, color = "orange", lineend = "round"
  ) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  geom_hline(yintercept = 7.5, linetype = "dotted", color = "gray40") +
  # geom_vline(xintercept = lubridate::ymd("2025-06-20"), linetype = "solid", color = "orange") +
  # geom_vline(xintercept = lubridate::ymd("2025-07-02"), linetype = "solid", color = "orange") +
  # geom_vline(xintercept = lubridate::ymd("2025-08-12"), linetype = "solid", color = "orange") +
  # geom_vline(xintercept = lubridate::ymd("2025-08-20"), linetype = "solid", color = "orange") +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.15, color = NA) +
  geom_line(size = 1) +
  facet_grid(
    rows = vars(culture),
    cols = vars(robinia),
    labeller = labeller(robinia = c(
      `without-robinia` = "without robinia",
      `with-robinia` = "with robinia"
    )),
    scales = "fixed"
  ) +
  scale_color_manual(values = c(control = "black", drought = "#D62728"), name = "Precipitation") +
  scale_fill_manual(values = c(control = "black", drought = "#D62728"), name = "Precipitation") +
  scale_linetype_manual(
    values = c(inoc_beech = "dashed", inoc_robinia = "solid"),
    labels = c(inoc_beech = "beech soil", inoc_robinia = "robinia soil"),
    name = "Soil type"
  ) +
  scale_x_date(date_labels = "%m/%y") +
  labs(
    x = "Date", y = "Soil water content (%)",
    title = "Soil water content over time"
  ) +
  coord_cartesian(ylim = c(-2, 30))


# Soil CN Content -----------------------------------------------------------------------------

df <- get_data("box", "cn_isotopes")
df |> glimpse()

# prepare data
df_long <- df %>%
  mutate(
    soiltype = gsub("-", "_", soiltype),
    soiltype = factor(soiltype, levels = c("inoc_beech", "inoc_robinia"))
  ) %>%
  select(soiltype, c_perc, n_perc, d13c_permille, d15n_permille) %>%
  pivot_longer(
    cols = c(c_perc, n_perc, d13c_permille, d15n_permille), # <- not soiltype
    names_to = "variable",
    values_to = "value"
  )

var_labels <- c(
  c_perc = "C (%)",
  n_perc = "N (%)",
  d13c_permille = "d13C (\u2030)",
  d15n_permille = "d15N (\u2030)"
)

ggplot(df_long, aes(x = soiltype, y = value, fill = soiltype)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1.3) +
  stat_compare_means(
    method = "t.test", label = "p.format",
    comparisons = list(c("inoc_beech", "inoc_robinia"))
  ) +
  facet_wrap(~variable,
    scales = "free_y",
    labeller = labeller(variable = var_labels),
    nrow = 2
  ) +
  scale_fill_manual(values = c(inoc_beech = "skyblue3", inoc_robinia = "#E64B35")) +
  labs(x = "Soil type", y = NULL, title = "Isotopes and elemental % by soil type")
