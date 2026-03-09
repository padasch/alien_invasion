# ALINV – Alien Invasion Experiment

Analysis code for the Alien Invasion experiment of the Ecosystem Ecology Group. The experiment investigates how drought, species mixing, *Robinia* presence, and legacy soil type affect tree physiology and growth in a full-factorial mesocosm design with *Fagus sylvatica* and *Quercus petraea*.

## Project Structure

```
functions/       # Reusable function libraries (sourced by notebooks)
scripts/         # Standalone executable scripts (data cleaning, model fitting)
notebooks/       # R Markdown analysis notebooks
data/
  raw/           # Raw data (not tracked – see data/raw/README.md)
  interim/       # Cleaned CSV files (tracked)
  final/         # Final outputs
output/          # Model caches and figures (not tracked)
_archive/        # Superseded files
```

### Notebooks

| File | Description |
|------|-------------|
| `1-treatment-effects.Rmd` | Main analysis: treatment effect time series, biomass GLMM, soil conditions, temporal GLMMs, and piecewise SEMs for all tree traits |
| `2-swc-interpolation.Rmd` | SWC GAM imputation and SEM sensitivity comparison (measured vs. imputed SWC) |

### Scripts

| File | Description |
|------|-------------|
| `1-data-cleaning.R` | Excel → CSV pipeline |
| `2-competition-on-growth.R` | Phenology and SLA exploratory analysis |
| `3-cleaning-sensor-data.R` | Sensor data processing |
| `4-impute-swc-gam.R` | Standalone SWC GAM imputation |

## Reproducing Results

This project uses [`renv`](https://rstudio.github.io/renv/) to manage R package dependencies.

1. **Clone** the repository and open `alien_invasion.Rproj` in RStudio.
2. **Obtain raw data** from the shared Microsoft Teams folder and place the files in `data/raw/` (see `data/raw/README.md`).
3. **Restore packages:**
   ```r
   renv::restore()
   ```
4. **Run data cleaning** to regenerate interim CSVs (if needed):
   ```r
   source("scripts/1-data-cleaning.R")
   ```
5. **Knit notebooks** from the `notebooks/` directory — they resolve the project root automatically.
