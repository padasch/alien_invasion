# Repeated-response SEM bootstrap check

## Current run

The source samples were rebuilt from the current tracked interim data before
bootstrapping. The growth models now contain the combined 1 September growth
round: 2,520 Fagus and 2,664 Quercus observations. The two Quercus senescence
responses remain unavailable because no usable response rows exist.

All 12 estimable response models reached 1,000 successful block-stratified
container-bootstrap refits without failed attempts. The run used eight workers
and took 255.9 seconds. The SWC models had no singular refits; 3,374 of 12,000
response-model refits were singular and are reported rather than discarded.

The updated event window provided no event-period observations for the two
Fagus senescence responses. Their event terms were therefore omitted instead
of fitting an unidentifiable contrast.

## Uncertainty comparison

Across 184 estimable treatment-path combinations, bootstrap inference changed
the probability-based significance classification for 5 direct, 5 indirect,
6 total, and 5 treatment-to-SWC effects relative to the point-model delta/Wald
approximation. The six changed path-summed totals were:

| Species | Response | Treatment | Estimate | Delta p | Bootstrap p | Bootstrap 95% interval |
|---|---|---|---:|---:|---:|---:|
| Fagus | Vitality | Extreme event | -0.063 | 0.015 | 0.266 | [-0.167, 0.057] |
| Quercus | Chlorophyll | Extreme event | -0.125 | 0.059 | 0.004 | [-0.229, -0.039] |
| Quercus | Vitality | Culture | 0.135 | 0.286 | 0.044 | [0.002, 0.358] |
| Quercus | Quantum yield | Precipitation | 0.185 | 0.084 | 0.028 | [0.028, 0.372] |
| Quercus | Volume (total) | Extreme event | -0.038 | 0.146 | 0.002 | [-0.060, -0.013] |
| Quercus | Volume (increment) | Extreme event | -0.136 | 0.003 | 0.318 | [-0.376, 0.128] |

## Interpretation

The bootstrap is design-aligned because it resamples whole containers within
block while retaining their trees and repeated observations. It estimates
uncertainty across experimental containers, not across alternative event
timings. The high singular-refit frequency for some response models remains a
limitation and should accompany any inferential interpretation of the
exploratory path decomposition.
