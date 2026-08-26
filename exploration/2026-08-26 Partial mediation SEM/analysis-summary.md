# Partial-mediation SEM exploration

## Question

For each imposed treatment, separate the overall standardized response into an SWC-adjusted component and an SWC-associated indirect component, while retaining their sum as the total treatment effect.

## Model

For treatment k, the two linked equations are:

`SWC_z = time terms + sum(a_k X_k) + random effects + error`

`Y_z = b SWC_z + time terms + sum(c'_k X_k) + random effects + error`

The SWC-adjusted component is `c'_k`, the SWC-associated indirect component is `a_k b`, and the path-summed total is `c'_k + a_k b`.

Repeated responses use same-day or the latest preceding measured SWC within seven days. Leaf-out uses mean pre-leaf-out SWC at container level (4 March to 2 April 2025). Negative effects denote deterioration or delayed leaf-out; positive effects denote improvement or earlier leaf-out.

## Bootstrap uncertainty

Yes, bootstrapping is part of this analysis. Containers were sampled with replacement within experimental block, retaining all associated trees and observations. Each estimable model contributed 1,000 successful refits. Direct, indirect, and total effects were recalculated within every replicate; the confidence interval for the total therefore comes from the empirical distribution of `c' + ab`.

The present exploration reuses the completed replicate-level fits and independently recalculates their percentile intervals and two-sided bootstrap probabilities. It does not perform a second set of identical model refits.

## Numerical QA

All bootstrap path-sum checks passed: the maximum absolute difference between the stored total and `c' + ab` was 4.440892e-16.

## Main findings

- Fagus: 21 of 29 SWC-adjusted effects, 10 of 29 indirect effects, and 19 of 29 total effects excluded zero.
- Quercus: 11 of 23 SWC-adjusted effects, 9 of 23 indirect effects, and 9 of 23 total effects excluded zero.
- Direct and indirect estimates opposed one another in 18 of 29 Fagus combinations and 11 of 23 Quercus combinations.
- The indirect estimate was larger in absolute magnitude than the SWC-adjusted estimate in 2 Fagus and 2 Quercus combinations.

## Interpretation

This is a clean additive decomposition for presentation: the total panel answers what each treatment did overall, while the other panels show how that estimate splits into an SWC-adjusted remainder and an SWC-associated pathway. However, it does not support a blanket claim that Fagus responses were mainly mediated by measured SWC. Many Fagus indirect paths oppose the corresponding SWC-adjusted paths, particularly for chlorophyll, vitality, and quantum yield. Those patterns are better described as statistical suppression or competing associations than as straightforward positive water mediation.

Because the response-SWC matching and the SWC-to-response path are observational within the experiment, `ab` should be called an SWC-associated indirect component rather than a proven causal mediation effect. A percentage mediated should only be reported when direct and indirect paths point in the same direction; the diagnostics table therefore leaves that quantity missing for opposing paths.

## Files

- `output/sem-total-and-path-decomposition.pdf`: shared-scale heatmap of total, SWC-adjusted, and indirect effects.
- `output/sem-effects-total-direct-indirect.csv`: complete effect table with bootstrap intervals and probabilities.
- `output/sem-decomposition-diagnostics.csv`: component dominance and path-direction diagnostics.
- `output/bootstrap-path-sum-qa.csv`: replicate counts and numerical identity checks.
