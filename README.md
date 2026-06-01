# ALINV – Alien Invasion Experiment

Analysis code for the Alien Invasion experiment of the Ecosystem Ecology Group. The experiment investigates how drought, species mixing, *Robinia* presence, and legacy soil type affect tree physiology and growth in a full-factorial mesocosm design with *Fagus sylvatica* and *Quercus petraea*.

## Project Structure

```
config/          # Shared analysis config (factor baselines, labels, scenarios)
functions/       # Reusable function libraries (sourced by notebooks)
scripts/         # Standalone executable scripts (data cleaning, model fitting)
notebooks/       # R Markdown analysis notebooks
data/
  raw/           # Raw data (not tracked – see data/raw/README.md)
  interim/       # Cleaned CSV files (tracked)
output/          # Model caches and figures (not tracked)
```

### Notebooks

| File | Description |
|------|-------------|
| `1-treatment-effects.Rmd` | Main analysis: summary figures, biomass GLMMs, temporal GLMMs, and direct-effect SEMs |
| `2-swc-interpolation.Rmd` | SWC GAM imputation and SEM sensitivity comparison (measured vs. imputed SWC) |
| `3-sem-aggregation.Rmd` | SEM effect heatmaps and CSV exports |
| `4-data-qc.Rmd` | Data-quality and metadata checks |
| `5-size-trajectories.Rmd` | Growth-only notebook for height, diameter, and volume trajectories plus model panels |
| `6-swp-provenance.Rmd` | Provenance notebook for the SWP/SWC derived product |

### Scripts

| File | Description |
|------|-------------|
| `1-data-cleaning.R` | Excel → CSV pipeline |
| `2-competition-on-growth.R` | Phenology and SLA exploratory analysis |
| `3-cleaning-sensor-data.R` | Sensor data processing |
| `4-impute-swc-gam.R` | Standalone SWC GAM imputation |
| `5-render-analysis-scenarios.R` | Render all notebooks across all soil scenarios |
| `6-plot-sem-heatmaps-from-csv.R` | Rebuild SEM heatmaps from exported CSV files |
| `7-plot-temporal-models-from-csv.R` | Rebuild temporal GLMM plots from exported CSV files |
| `8-plot-sem-models-from-csv.R` | Rebuild SEM path plots from exported CSV files |
| `run_all.R` | Convenience entry point for the default scenario plus the global SWP notebook |

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
5. **Render the default analysis path**:
   ```r
   source("run_all.R")
   ```
6. **Or render all scenarios**:
   ```r
   source("scripts/5-render-analysis-scenarios.R")
   ```
