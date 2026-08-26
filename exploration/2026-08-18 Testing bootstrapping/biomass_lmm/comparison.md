# Harvest-biomass bootstrap comparison

## Result

All 18 treatment coefficients retained their sign; 2 changed their 95%-interval significance classification. The current Wald analysis classified 10/18 effects as significant, compared with 12/18 using percentile container-cluster bootstrap intervals.

Bootstrap-to-Wald CI-width ratios ranged from 0.85 to 1.01 (median 0.96).

## Method

For each species and biomass metric, whole containers were sampled with replacement within experimental block. All trees from each selected container were retained and duplicated clusters received synthetic container and tree identifiers. The original Figure 5 model was then refitted: standardized response ~ precipitation + culture + Robinia + (1 | container). Response standardization was fixed to the original model sample. Percentile 95% intervals and empirical two-sided probabilities were calculated from 1,000 successful replicates.

## Fit diagnostics

The six models produced 6000 successful bootstrap refits (1000 per model) with 0 failed attempts. Five of six point models were singular; 4479 successful bootstrap refits estimated a zero or near-zero container variance.

The analysed harvest workbook contained 32 containers from two represented blocks (b2 and b3) per species; block b1 had no rows in this biomass dataset. Model samples contained 128 Fagus trees and 127 Quercus trees.

Singularity is scientifically informative here: for most biomass responses, between-container residual variance was estimated at zero after the additive treatment effects. Cluster resampling nevertheless retains the container as the treatment-assignment and dependence unit.

## Classification changes

- fagus, root_biomass, precipitationdrought: Wald not significant; bootstrap significant (bootstrap p = 0.052).
- fagus, root_biomass, culturemixed: Wald not significant; bootstrap significant (bootstrap p = 0.052).

## Interpretation

Bootstrap and Wald estimates use the same fitted model and therefore answer the same additive treatment-effect question. Differences reflect uncertainty estimation, not a change in the estimand. Because treatments were assigned to containers, the block-stratified container bootstrap is the more design-aligned sensitivity analysis. Results remain descriptive of the fitted factorial main effects; treatment interactions were not tested.

Point-model singularity count: 5/6.
