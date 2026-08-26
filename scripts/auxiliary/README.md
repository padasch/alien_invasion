# Auxiliary analysis code

This directory contains only code shared by the active publication pipelines
or needed to rebuild the cleaned interim data.

## Data preparation

- `1-data-cleaning.R`: rebuilds the main cleaned experiment tables.
- `3-cleaning-sensor-data.R`: cleans the soil-sensor data.
- `4-impute-swc-gam.R`: creates the optional GAM-based SWC interpolation.
- `config/analysis-config.R`: central factor levels, response orientation,
  treatment definitions, and volume-proxy settings.

## Shared functions

- `_source.R`: project paths, configuration loading, and core data access.
- `1-summary-figures.R`: descriptive summaries and measurement schedule.
- `2-biomass.R`: biomass data preparation used by the bootstrap and
  supplementary workflows.
- `3-effect-size-factorial.R`: temporal model-data preparation used by the
  bootstrap and supplementary workflows.
- `7-phenology-transition-models.R`: common-shift phenology model.
- `8-size-trajectories.R`: height/diameter/volume trajectory preparation.
- `10-temporal-alignment.R`: past-only SWC-to-response date alignment.
- `11-bootstrap-inference.R`: bootstrap-result readers and interval helpers.

## Bootstrap analyses

`bootstrap/` contains the production model-refitting scripts for the temporal,
biomass, phenology, and SWC-pathway analyses. Compact summaries consumed by the
figures are tracked under `data/final/bootstrap/`; large replicate-level and
model-cache artifacts are regenerated locally and are not publication outputs.

The publication scripts source only the files needed for their collection.
Superseded model and plotting helpers were moved to the ignored `_archive/`
tree on 2026-08-24 instead of remaining on the active source path.
