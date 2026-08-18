# Phenology bootstrap comparison

All models used **1000 successful non-parametric bootstrap replicates per species**. Containers were sampled with replacement within block; all trees and transition records in each selected container were retained under synthetic IDs.

## Primary common-shift LMM

- Direction changes: 0 of 6.
- 95% interval-exclusion changes: 0 of 6.
- P < 0.05 classification changes: 0 of 6.
- Point estimates are unchanged by construction; only uncertainty estimation differs.

## Stage-specific and duration sensitivities

- Stage-specific contrasts: 0 direction, 1 interval, and 1 P-value classification changes among 18 contrasts.
- Stage 2–4 duration: 0 direction, 0 interval, and 0 P-value classification changes among 6 contrasts.

- The only stage-specific threshold change is *Fagus*, mixed culture, stage 4: Kenward–Roger 95% CI [-2.308, 0.096], P = 0.0706; bootstrap 95% CI [-2.184, -0.014], P = 0.0500. This is a borderline sensitivity result, not the primary common-shift result.

## Phenology SEM

- Direction changes: 0 of 24 component effects.
- 95% interval-exclusion changes: 1 of 24.
- P < 0.05 classification changes: 1 of 24.
- The only threshold change is the *Quercus* reduced-precipitation direct path: unstratified 95% CI [-0.788, 0.024], P = 0.0659; block-stratified 95% CI [-0.772, -0.005], P = 0.0480.
- The Figure 6 target—the six path-summed total effects—has 0 threshold changes among 6; all six intervals still include zero.
- Point estimates are unchanged; differences arise from preserving block composition in every bootstrap sample.

## Fit accounting

- LMM bootstrap failures before reaching the target: 0.
- SEM bootstrap failures before reaching the target: 0.
- Singular primary-LMM fits: 416 of 2000; singular SEM outcome fits: 410 of 2000.
- Singular fits were retained because the random-intercept variance can reach zero after cluster resampling while the requested fixed treatment contrasts remain estimable. This high frequency, especially for *Fagus*, is an important bootstrap caveat rather than a convergence failure.

## Interpretation

Direction is based on the manuscript-oriented sign: negative values indicate delayed/slower phenology and positive values earlier/faster phenology. Bootstrap intervals are percentile intervals and P values follow the pre-specified sign-count formula in the exploration README.

Overall, the bootstrap analysis agrees with the current phenology conclusions: treatment-effect directions are unchanged, all six primary common-shift effects remain non-significant, all duration effects retain their classification, and all six Figure 6 path-summed totals remain non-significant. The two isolated P ≈ 0.05 changes should be reported as sensitivity-level borderline findings.
