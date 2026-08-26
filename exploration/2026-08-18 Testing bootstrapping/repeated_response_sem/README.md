# Repeated-response SEM bootstrap

This directory tests container-cluster bootstrap inference for the 14
non-phenology species-response slots used in the supplementary SEM analysis.
Before fitting, `00_refresh_response_model_sources.R` reconstructs every
estimable model sample from the current tracked interim data and applies the
prespecified additive mediator and response formulas. Model selection is not
repeated within bootstrap samples.

## Method

- Experimental unit: container (`boxlabel`).
- Resampling: the original number of containers is sampled with replacement
  within each experimental block (`b1`, `b2`, `b3`).
- Dependence: all trees, dates, and aligned SWC observations from a sampled
  container remain together. Repeated draws receive new container and tree IDs.
- Standardization: response, SWC, and centred DOY variables are standardized
  once in the current original-data sample; bootstrap samples are not
  restandardized.
- Paths: direct (`c'`), SWC-associated indirect (`a*b`), and path-summed total
  (`c' + a*b`) effects are recalculated in each replicate.
- Uncertainty: 95% percentile intervals and the shared two-sided empirical
  bootstrap probability are calculated from 1,000 successful replicates.
- Display orientation: both raw coefficient signs and the manuscript-oriented
  signs are exported.

The two Quercus senescence slots (`remaining_green` and `chlavg`) contain no
post-filtering SEM data in the current interim data. They are retained in
the source inventory and model status table as unavailable, matching the
current figure's missing cells.

## Reproduce

From the repository root:

```sh
Rscript --vanilla \
  "exploration/2026-08-18 Testing bootstrapping/repeated_response_sem/scripts/01_repeated_response_sem_bootstrap.R" \
  --bootstrap=1000 --cores=8 --seed=20260818
```

## Outputs

All generated files are under `output/`:

- `repeated-response-sem-source-inventory.csv`: all 14 response slots and
  their refreshed current-data source files.
- `repeated-response-sem-bootstrap-status.csv`: formulas, data sizes,
  successful replicates, attempts, failures, and singular-fit counts.
- `repeated-response-sem-point-delta-effects.csv`: current point estimates and
  delta/Wald uncertainty.
- `repeated-response-sem-bootstrap-replicates.csv`: full replicate-level path
  effects; the equivalent objects are also stored in the results RDS.
- `repeated-response-sem-bootstrap-effects.csv`: percentile intervals and
  empirical probabilities.
- `repeated-response-sem-delta-vs-bootstrap-comparison.csv`: effect-level
  comparison of the old and bootstrap inference.
- `figure6-ready-repeated-response-sem-path-summed-total.csv`: bootstrap total
  effects for integration into the recreated Figure 6.
- `repeated-response-sem-bootstrap-failures.csv`: failed attempts, if any.
- `repeated-response-sem-bootstrap-results.rds`: complete reproducibility
  bundle.

See `comparison.md` for the result-level QA and interpretation.
