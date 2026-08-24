# Exact-date SWC sensitivity analysis

This analysis tests whether the repeated-response SEM decomposition changes
when responses are paired only with measured soil water content (SWC) from the
same container on the exact same date.

## Method

- The current fuzzy-match SEM RDS bundles supply the analysis data, treatment
  coding, response definitions, and selected mediator and response formulas.
- Rows are retained by an inner join to `data/interim/box_soilwater.csv` on
  `boxlabel` and `date`.
- The raw SWC in every retained row is checked against cached `swc_org`; the
  maximum discrepancy was zero.
- Cached standardized `swc`, `y`, `doy_c`, and `doy_c2` values are retained.
  The reduced samples are not re-standardized. Recovered baseline centers,
  scales, and reconstruction errors are exported.
- Original formulas are frozen. A formula is not simplified when exact-date
  restriction removes factor variation or makes a random effect inestimable.
- Containers are sampled with replacement within blocks. All measurements and
  trees from a selected container remain together, and duplicate draws receive
  new container and tree IDs.
- Each estimable SEM uses 1,000 successful container-cluster bootstrap
  replicates. Percentile intervals and empirical two-sided bootstrap
  probabilities use the shared functions in `functions/11-bootstrap-inference.R`.

## Reproduce

From the repository root:

```sh
Rscript --vanilla \
  "exploration/2026-08-18 Testing bootstrapping/swc_matching_sensitivity/exact_date/scripts/01_exact_date_sem_bootstrap.R" \
  --bootstrap=1000 --cores=4 --seed=2026081811
```

## Outputs

- `output/model-status.csv`: estimability, sample sizes, formulas, bootstrap
  completion, failure, and singular-fit counts.
- `output/sample-loss.csv`: baseline and exact-date rows, containers, trees,
  dates, and retained fractions.
- `output/baseline-standardization-constants.csv`: recovered scaling constants
  and reconstruction checks.
- `output/treatment-balance.csv`: factor-level proportions before and after
  exact-date restriction.
- `output/exact-date-point-effects.csv`: exact-date point paths, including the
  SWC-to-response coefficient `b`.
- `output/exact-date-bootstrap-effects.csv`: exact-date bootstrap estimates,
  intervals, and probabilities for `a`, `b`, direct, indirect, and total paths.
- `output/exact-date-vs-fuzzy-comparison.csv`: effect-level comparison with the
  current fuzzy-match bootstrap baseline.
- `output/path-identity-qa.csv`: point and replicate checks of
  `total = direct + indirect`.
- `output/exact-date-vs-fuzzy-sem-paths.pdf`: compact visual comparison for the
  four estimable species-response models.
- `output/exact-date-sem-bootstrap-results.rds`: complete reproducibility
  bundle, including replicate-level estimates.

See `comparison.md` for findings and limitations.
