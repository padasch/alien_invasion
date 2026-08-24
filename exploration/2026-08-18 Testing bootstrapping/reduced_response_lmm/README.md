# Reduced-form response models

This analysis supplies the standardized treatment-effect heatmaps in Figure 6.
For each available species-response combination, it retains the covariates,
random effects, response scaling, and analysis sample of the corresponding
repeated-response model but removes SWC from the response formula. Treatment
effects are therefore estimated directly rather than reconstructed as SEM path
sums.

Uncertainty is estimated with a block-stratified container-cluster bootstrap.
Complete containers, including all trees and repeated observations, are
resampled within experimental block. The publication workflow currently
requires 1,000 successful refits per estimable model.

Run from the project root:

```sh
Rscript --vanilla \
  "exploration/2026-08-18 Testing bootstrapping/reduced_response_lmm/scripts/01_reduced_response_bootstrap.R" \
  --bootstrap=1000 --cores=4
```

The script also compares the directly fitted estimates with the former
path-summed SEM totals. That comparison is diagnostic only; Figure 6 reads the
reduced-form estimates.
