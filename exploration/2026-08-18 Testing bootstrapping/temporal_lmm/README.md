# Temporal LMM bootstrap exploration

This directory changes only the uncertainty estimator for the manuscript's
date-specific diameter-increment and quantum-yield LMMs. Point estimates use
the same formulas, data filters, factor baselines, and original species-level
response z-scaling as the current analysis.

Run from the repository root:

```sh
RENV_CONFIG_SANDBOX_ENABLED=FALSE Rscript --vanilla \
  "exploration/2026-08-18 Testing bootstrapping/temporal_lmm/scripts/01_temporal_lmm_bootstrap.R"
```

The script requests 1,000 successful bootstrap replicates per
species-response model. It samples containers with replacement within block,
retains all nested observations, assigns synthetic container/tree identifiers
to duplicate selections, and reports percentile 95% intervals and empirical
two-sided bootstrap probabilities.

Outputs are under `output/`; the numerical validation is in `comparison.md`.
