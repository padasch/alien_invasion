# Exact-date versus fuzzy SWC matching

## Feasibility and sample loss

Only four of the 14 repeated-response SEM slots were estimable with exact-date
SWC while preserving the original formulas:

| Model | Baseline rows | Exact-date rows | Retained | Containers | Trees |
|---|---:|---:|---:|---:|---:|
| Fagus chlorophyll | 1,680 | 1,120 | 66.7% | 47 | 140 |
| Fagus quantum yield | 1,524 | 839 | 55.1% | 47 | 140 |
| Quercus chlorophyll | 1,776 | 1,184 | 66.7% | 49 | 148 |
| Quercus quantum yield | 1,628 | 888 | 54.5% | 49 | 148 |

Vitality had no response observations on measured-SWC dates. Fagus
senescence retained 15.4--22.0% of rows, but the exact-date samples contained
only one extreme-event level, so the frozen formulas were not estimable.
Growth retained 12.5% of rows, likewise lost extreme-event variation, and had
only one retained observation per tree, which also prevents estimation of the
tree random intercept. The two Quercus senescence slots were already
unavailable in the fuzzy-match baseline. These models were recorded as
unavailable rather than simplified.

## Comparison of estimable models

Among Fagus and Quercus chlorophyll and quantum yield:

- All 16 direct-effect directions agreed with the fuzzy baseline.
- All 16 indirect-effect directions agreed. None changed bootstrap
  significance, although magnitudes sometimes changed substantially.
- Fifteen of 16 total-effect directions agreed.
- Two direct effects and two total effects changed significance at 0.05.
- The exception in total-effect direction was the Fagus quantum-yield extreme-
  event estimate: it changed from -0.265 (95% CI -0.340 to -0.185) under fuzzy
  matching to +0.091 (0.017 to 0.165) under exact-date matching. Its indirect
  component increased from +0.286 to +0.701, outweighing the negative direct
  component in the exact-date subset.
- Quercus quantum-yield precipitation remained positive in direction, but its
  direct and total intervals included zero after exact-date restriction.
- The Quercus chlorophyll extreme-event total and Quercus quantum-yield
  extreme-event direct effect also lost significance.

This provides qualified support for the original qualitative decomposition in
the response types that can be tested: direct and indirect directions are
stable, but several magnitudes, net totals, and uncertainty conclusions are
sensitive to temporal matching. It does not establish robustness for vitality,
growth, or senescence because exact-date data cannot fit their frozen models.

## Treatment balance and time selection

Exact-date restriction barely changed the randomized treatment composition in
the four fitted models: the largest absolute proportion changes were 0.003 for
precipitation, 0.005 for Robinia, and 0.002 for culture. It changed the
extreme-event proportion by as much as 0.142 because exact measurement dates
are not uniformly distributed through time. The extreme-event sensitivity is
therefore plausibly driven partly by date selection rather than treatment
imbalance.

## Bootstrap and QA

All four models reached 1,000 successful replicates without failed attempts.
The response submodel was singular in 11 Fagus quantum-yield, 50 Quercus
chlorophyll, and 377 Quercus quantum-yield replicates; this boundary behavior
should temper inference for those models, especially Quercus quantum yield.
No SWC submodel replicate was singular.

For every point estimate and bootstrap replicate,
`total - direct - indirect` was checked. The maximum absolute numerical error
was `1.11e-16`. Baseline standardization was preserved exactly rather than
recomputed on the smaller samples; reconstruction errors were below `4e-14`.

## Interpretation

Exact-date matching avoids using SWC measured after a response and removes
temporal mismatch, but it is not a generally usable primary analysis here
because it discards too much information and selectively samples dates. Its
main value is as a sensitivity check. For chlorophyll and quantum yield, the
SWC pathway is qualitatively similar but not fully quantitatively stable. A
preceding-only match and climate-assisted daily interpolation are needed to
assess whether that instability reflects sample restriction or exposure
measurement error.
