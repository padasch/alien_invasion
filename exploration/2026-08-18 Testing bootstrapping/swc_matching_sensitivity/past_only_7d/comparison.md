# Past-only 7-day SWC sensitivity analysis

## Analysis

For each repeated-response SEM, a response observation was retained only when measured SWC was available on the same day or during the preceding seven days. Future SWC observations and the original unbounded past fallback were excluded. The original fitted bundle supplied the frozen model formulas, factor coding, and existing response and SWC standardisation; the restricted sample was not re-standardised.

Uncertainty was estimated using 1,000 successful block-stratified container-cluster bootstrap replicates per estimable species-response model. Containers were sampled with replacement within experimental blocks. The direct, indirect (`a * b`), total (`direct + indirect`), treatment-to-SWC (`a`), and SWC-to-response (`b`) paths were recalculated in every replicate.

## Sample and fit diagnostics

- All 12 baseline-estimable repeated-response SEMs remained estimable. The two Quercus senescence models remained unavailable because their baseline bundles contained no post-filtering data.
- Sample loss ranged from 18.2% to 54.5%. Quantum-yield models lost about 18%; chlorophyll, volume, and volume-increment models about 25%; Fagus senescence 38.5%; and vitality 54.5%.
- Retained matches spanned lags of 0–7 days. Across reused model records, same-day observations were the most common, but the exact proportions should not be interpreted as independent observations because schedules recur across models.
- Precipitation, Robinia, and culture balance was generally preserved. Exclusion varied much more between extreme-event and non-event dates, because matching availability is intrinsically date-dependent.
- Every estimable model reached 1,000 successful replicates with no fitting failures. SWC submodels were never singular; response submodels were singular in 3,387 of 12,000 successful fits, concentrated in models already prone to boundary random-effect estimates.
- All 48,000 replicate-level checks satisfied `total = direct + indirect` to numerical precision (maximum absolute discrepancy < 1.2e-16).

## Comparison with the current fuzzy-match baseline

Across the 48 treatment-response combinations:

| Component | Direction changes | Significance changes | Median absolute estimate change | Non-overlapping 95% intervals |
|---|---:|---:|---:|---:|
| Direct | 3 | 6 | 0.031 SD | 2 |
| Indirect via SWC | 8 | 10 | 0.023 SD | 6 |
| Total | 3 | 7 | 0.020 SD | 5 |

The qualitative total-effect conclusions for the two main imposed treatments were comparatively stable:

- Robinia: all 12 total-effect directions were retained, and none changed bootstrap significance status.
- Precipitation: all 12 total-effect directions were retained; one changed significance status.
- Culture: changes were small and mostly concerned estimates near zero; two total-effect significance classifications changed and one near-zero total crossed zero.
- Extreme event: this was the least stable term, with four total-effect significance changes and two direction changes. This sensitivity is consistent with the strong date-dependent loss of event observations.

The pathway decomposition was less stable than the net total effects. For example, the Quercus quantum-yield extreme-event effect remained clearly negative in total, but the negative component shifted from the direct path in the fuzzy analysis to the SWC-mediated path in the past-only analysis. Fagus senescence under the extreme-event term was more consequential: its total changed from significantly negative to small positive/borderline non-significant.

## Interpretation

The main Robinia and precipitation total-effect story is broadly robust to prohibiting future SWC matches. However, attribution between direct and SWC-associated indirect pathways is materially more sensitive, and extreme-event pathway estimates are especially dependent on the matching rule. The fuzzy baseline should therefore not be used alone to make strong claims about mediation. The past-only result is temporally defensible but loses substantial data—particularly for vitality—so it should be considered alongside exact-date, interpolated-daily, and antecedent-window analyses.

## Outputs

- `output/past-only-sem-bootstrap-status.csv`: model, sample-loss, bootstrap, failure, and singularity status.
- `output/past-only-lag-diagnostics.csv`: retained lag distribution.
- `output/past-only-sample-balance-diagnostics.csv`: treatment/date balance before and after restriction.
- `output/past-only-sem-bootstrap-effects.csv`: point estimates, bootstrap intervals, and p-values for `a`, `b`, direct, indirect, and total effects.
- `output/past-only-vs-fuzzy-bootstrap-comparison.csv`: effect-level comparison against the current fuzzy baseline.
- `output/past-only-qa-total-identity.csv`: numerical path-sum QA.
- `output/past-only-vs-fuzzy-effect-deltas.pdf`: diagnostic heatmap of effect changes.
