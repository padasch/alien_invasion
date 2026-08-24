# Antecedent SWC-window sensitivity

Run from the repository root:

```sh
Rscript "exploration/2026-08-18 Testing bootstrapping/swc_matching_sensitivity/antecedent_windows/scripts/01_antecedent_window_sem_bootstrap.R" --bootstrap=1000 --cores=4
```

The script compares modeled trailing 7-day and 14-day mean SWC mediators with the fuzzy-matched baseline. Both windows end on the response date and never use future information. Outputs include coverage, status, point and bootstrap effects, all replicates, pairwise comparisons, path-identity QA, a compact PDF, and `comparison.md`.
