# Phenology bootstrap exploration

This directory contains the block-stratified container-bootstrap version of
the overall phenology timing analysis. It does not modify the manuscript's
current figures, supplementary files, scripts, or dated output.

Run from the project root:

```sh
Rscript --vanilla "exploration/2026-08-18 Testing bootstrapping/phenology/run_phenology_bootstrap.R" --bootstrap=1000 --cores=3
```

The script refits, separately for each species:

- the primary common-shift Gaussian LMM for stages 2–4;
- the stage × treatment sensitivity LMM;
- the stage-2-to-stage-4 duration LMM;
- the associative phenology SEM using the existing stage-centred timing index
  and pre-leaf-out SWC definition.

Containers are sampled with replacement within block. Every selected
container retains all trees and transition observations; duplicated copies
receive synthetic container and tree identifiers. Output tables and model
bundles are in `output/`, plots are in `figures/`, and comparisons with the
current analysis are in `comparison/`.
