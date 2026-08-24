# Daily-interpolation SWC sensitivity

This analysis replaces fuzzy response-to-SWC matching in the repeated-response
SEMs with exact response-date SWC. It retains observed SWC when available on the
exact date and otherwise uses the existing treatment-agnostic container-level
daily GAM prediction.

Run from the repository root:

```sh
Rscript "exploration/2026-08-18 Testing bootstrapping/swc_matching_sensitivity/daily_interpolation/scripts/01_daily_interpolation_sem_bootstrap.R" --bootstrap=1000 --cores=4
```

The script audits the daily interpolation, performs leave-measurement-date-out
validation, preserves the original model formulas and scaling, runs the same
block-stratified container-cluster bootstrap used for the fuzzy baseline, and
exports effect-level comparisons and diagnostics under `output/`.
