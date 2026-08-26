# ALINV - Alien Invasion Experiment

Analysis and publication-figure code for the Alien Invasion experiment of the
Ecosystem Ecology Group. The experiment tests the effects of reduced
precipitation, species mixing, *Robinia* presence, and legacy soil type on
*Fagus sylvatica* and *Quercus ilex*.

## Repository structure

```text
scripts/
  main_figures/          # One script per main figure and the main runner
  supplementary_figures/# One script per supplementary figure and its runner
  auxiliary/             # Shared config, functions, and data-processing scripts
  run_all.R              # Rebuild both publication figure collections
data/
  raw/                   # Raw inputs; see data/raw/README.md
  interim/               # Cleaned analysis-ready CSV files
  final/bootstrap/       # Promoted bootstrap summaries used by figures
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
under `scripts/auxiliary/bootstrap/`.

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
  Rscript --vanilla scripts/supplementary_figures/make-v1-figures.R
```

The supplementary runner creates the 15 reconstructed supplementary figures.
`scripts/run_all.R` additionally creates the overall phenology-effect and
exploratory SWC-pathway figures.

## Rebuilding cleaned data

Cleaned interim files are tracked, so this is normally unnecessary. To rebuild
them from the raw inputs, run:

```sh
Rscript --vanilla scripts/auxiliary/1-data-cleaning.R
Rscript --vanilla scripts/auxiliary/3-cleaning-sensor-data.R
Rscript --vanilla scripts/auxiliary/4-impute-swc-gam.R
```

Shared data access, bootstrap readers, plotting utilities, factor definitions,
and volume-proxy configuration are under `scripts/auxiliary/`.
