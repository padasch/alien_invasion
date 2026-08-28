# ALINV - Alien Invasion Experiment

Analysis and publication-figure code for the Alien Invasion experiment of the
Ecosystem Ecology Group. The experiment tests the effects of reduced
precipitation, species mixing, *Robinia* presence, and legacy soil type on
*Fagus sylvatica* and *Quercus ilex*.

## Repository structure

```text
scripts/
  config.R               # Central analysis settings and factor definitions
  data_preparation/      # Executable scripts that rebuild cleaned data
  functions/             # Shared function definitions only
  model_fitting/         # Executable scripts that refit bootstrap models
  main_figures/          # One script per main figure and the main runner
  supplementary_figures/# One script per supplementary figure and its runner
  run_all.R              # Rebuild both publication figure collections
data/
  raw/                   # Raw inputs; see data/raw/README.md
  interim/               # Cleaned analysis-ready CSV files
  final/bootstrap/       # Generated local bootstrap summaries (not tracked)
  final/figure_data/     # Tidy tables containing the values shown in figures
output/
  main_figures/          # Final main-figure PDFs
  supplementary/         # Final supplementary-figure PDFs
docs/                    # Local manuscript files; intentionally not version controlled
_archive/                # Local historical material; intentionally not version controlled
```

Historical dated model caches are stored locally under
`_archive/model-caches/dated-output/`. They are ignored by Git and are retained
only for legacy method comparisons; `output/` contains publication PDFs only.
Local scratch analyses may be kept under `exploration/`, but that directory is
ignored and is not part of the repository. Production bootstrap scripts are
under `scripts/model_fitting/`. The retired GAM-based SWC interpolation is
stored locally under `_archive/gam-swc/`.

## Rebuilding the publication figures

The project uses [`renv`](https://rstudio.github.io/renv/) for package
dependencies.

```r
renv::restore()
```

From the project root, rebuild both collections with:

```sh
Rscript --vanilla scripts/run_all.R
```

Run a collection separately with:

```sh
Rscript --vanilla scripts/main_figures/make_all_figures.R
ALINV_SUPP_BOOT_B=1000 ALINV_SUPP_BOOT_CORES=8 \
  Rscript --vanilla scripts/supplementary_figures/make_all_figures.R
```

The supplementary runner creates all 18 supplementary figures, including
phenology, the exploratory SWC pathways, and the temporal effects on phase-wise
relative volume increments (S18).

Every figure script also refreshes one or more tidy CSV tables under
`data/final/figure_data/`. These contain the scientific values shown in the
figure, without plot styling or ggplot rendering details.

## Rebuilding cleaned data

Cleaned interim files are tracked, so this is normally unnecessary. To rebuild
them from the raw inputs, run:

```sh
Rscript --vanilla scripts/data_preparation/clean_experiment_data.R
Rscript --vanilla scripts/data_preparation/clean_sensor_data.R
```

Shared data access, bootstrap readers, and plotting utilities are under
`scripts/functions/`; analysis settings are in `scripts/config.R`.
