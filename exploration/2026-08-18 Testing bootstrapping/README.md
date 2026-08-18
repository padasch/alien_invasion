# Testing container-level bootstrap inference

This exploration compares the manuscript's current analytical results with
non-parametric container-cluster bootstrap uncertainty estimates. It does not
replace or overwrite files under `final_figures/`, `supplementary_information/`,
or the dated `output/` analysis directories.

## Shared bootstrap conventions

- Target: 1,000 **successful** bootstrap replicates for every fitted
  species-response model.
- Resampling unit: the experimental container (`boxlabel`). All trees,
  measurement dates, phenological transitions, and aligned SWC observations
  belonging to a selected container are kept together.
- Design preservation: containers are sampled with replacement within
  experimental block. Duplicated containers receive synthetic container and
  tree identifiers before model refitting.
- Model comparison: the fixed- and random-effect formulas used for point
  estimates are retained so that differences isolate the uncertainty method.
  For the older repeated-response SEMs, treatment terms selected on the
  original data are held fixed during bootstrap refitting; model selection is
  not repeated inside each replicate.
- Intervals: 2.5th and 97.5th percentiles of the successful bootstrap
  estimates.
- Two-sided bootstrap probability:

  `2 * min((number <= 0 + 1) / (B + 1), (number >= 0 + 1) / (B + 1))`,

  capped at one.
- Seeds are fixed and recorded by each model-family script.
- Both raw biological signs and the manuscript's oriented display signs are
  retained in source-data exports.

## Model-family directories

- `temporal_lmm/`: date-specific treatment contrasts for repeated responses.
- `biomass_lmm/`: harvest biomass treatment effects.
- `phenology/`: common-shift phenology LMM, sensitivity results, and phenology
  SEM.
- `repeated_response_sem/`: SWC-associated piecewise SEMs for the non-phenology
  responses used in Figure 6.
- `recreated_figures/`: bootstrap versions of the main and supplementary
  figures.
- `comparison/`: machine-readable comparisons and the final validation report.

## Reproduce the exploration

Run the four model-family scripts first, then assemble the figure collection:

```sh
Rscript --no-init-file temporal_lmm/scripts/01_temporal_lmm_bootstrap.R
Rscript --no-init-file biomass_lmm/scripts/01-biomass-container-bootstrap.R
Rscript --no-init-file phenology/run_phenology_bootstrap.R
Rscript --no-init-file repeated_response_sem/scripts/01_repeated_response_sem_bootstrap.R
Rscript --no-init-file recreated_figures/assemble_bootstrap_collection.R
```

The `--no-init-file` flag avoids waiting on a project-library activation lock
when independent bootstrap families are run concurrently. Each script records
its seeds, sample counts, failures, and singular fits in its own output folder.

The integrated scientific interpretation is in
`comparison/validation-report.md`, with its compact machine-readable summary in
`comparison/family-summary.csv`.

The descriptive time-series panels are reproduced from the same data and are
expected to be unchanged because they do not use model-based confidence
intervals.
