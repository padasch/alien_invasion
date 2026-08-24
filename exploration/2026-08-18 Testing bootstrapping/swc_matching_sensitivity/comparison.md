# SWC alignment sensitivity of the repeated-response SEMs

## Question and common design

This exploration tests whether the original response-to-SWC date matching makes the SEM underestimate SWC-associated indirect effects. All branches retain the original response-specific fixed-effect formulas, factor coding, response scaling, and SWC scaling. Uncertainty is estimated with 1,000 successful block-stratified container-cluster bootstrap replicates for every estimable model. The two Quercus senescence models remain unavailable because no post-filtering response data exist in the baseline analysis.

The original measured-SWC matching first selects the latest SWC observation in the preceding seven days, then the earliest observation in the following seven days, and finally the latest earlier observation without a strict maximum lag. Across model-by-container-date records, 29.6% were same-day, 39.7% used SWC measured 1–7 days earlier, 29.8% used SWC measured 1–7 days later, and 0.9% used SWC measured 8–12 days earlier.

## Sensitivity definitions

| Method | SWC definition | Main strength | Main limitation |
|---|---|---|---|
| Fuzzy baseline | Original measured-SWC rule | Retains the response sample | Includes future matches and an unbounded past fallback |
| Past-only measured | Same-day or latest preceding measured SWC within 7 days | Preserves temporal ordering and measured SWC | Loses 18.2–54.5% of rows |
| Exact-date measured | Same container and date only | No temporal mismatch or interpolation | Only four frozen-formula models remain estimable |
| Exact-date daily GAM | Observed same-day SWC, otherwise daily GAM estimate | Full exact-date coverage | Conditional on a moderately predictive, treatment-agnostic GAM |
| Antecedent 7-day GAM mean | Modeled response day through day -6 mean | Represents recent integrated exposure without future information | Uses modeled rather than measured daily SWC |
| Antecedent 14-day GAM mean | Modeled response day through day -13 mean | Tests a longer exposure window | Same interpolation limitation and less plausible for fast responses |

## Climate-informed interpolation

The existing daily GAM already uses measured container SWC, site-level soil water potential, air temperature, VPD, ambient precipitation, radiation, a date smooth, and a container random intercept. It is climate-informed but not treatment-informed. Ambient climate is shared by all containers and cannot reconstruct the experimental precipitation manipulation between SWC observations.

Leave-measurement-date-out validation gives RMSE 5.42 SWC percentage points, MAE 4.15, negligible mean bias, correlation 0.712, and predictive R-squared 0.40. This supports its use as a sensitivity analysis, but not as a definitive reconstruction of the mediator. Bootstrap intervals are conditional on the predicted SWC series and do not propagate GAM prediction uncertainty.

Adding further ambient climate variables is not recommended. A defensible treatment-aware interpolation would require dated container-level irrigation, exclusion, or water-input records.

## Results

### Total treatment effects

- Past-only measured SWC retained all Robinia total-effect directions and significance classifications. All precipitation directions were retained, with one significance change. Culture and especially the time-defined extreme-event term were less stable because exclusion was date-dependent.
- Exact-date measured SWC was estimable only for Fagus and Quercus chlorophyll and quantum yield. Direct and indirect directions agreed in all 16 comparisons; 15 of 16 total-effect directions agreed. The exception was the Fagus quantum-yield extreme-event total.
- Exact-date daily GAM SWC retained all 48 total-effect directions, 47 of 48 significance classifications, and overlapping intervals in all 48 comparisons.
- Both antecedent-window definitions retained all 48 total-effect directions. The 7-day window retained 47 of 48 significance classifications and the 14-day window retained all 48.

Stable path-summed totals in analyses using the same rows are not independent evidence of robustness. In a linear path model with the same covariates, `direct + a*b` is nearly the corresponding reduced-form treatment coefficient by regression algebra. Across the antecedent point fits, the median absolute difference from the reduced-form coefficient was 0.0001 standardized units and the maximum was 0.0109. The past-only and exact-date restrictions are therefore more informative total-effect checks because they change the analyzed observations.

### Direct versus SWC-associated indirect effects

| Method | Direct larger than indirect | Indirect share of summed absolute component magnitude | Opposing direct/indirect signs |
|---|---:|---:|---:|
| Fuzzy baseline | 87.5% | 16.5% | 70.8% |
| Past-only measured, 0–7 days | 91.7% | 21.7% | 68.8% |
| Exact-date measured | 87.5% | 30.6% | 68.8% |
| Exact-date daily GAM | 91.7% | 17.9% | 72.9% |
| Antecedent GAM mean, 7 days | 95.8% | 19.8% | 60.4% |
| Antecedent GAM mean, 14 days | 95.8% | 10.1% | 54.2% |

The broad direct-dominance result is therefore not created solely by fuzzy matching. It persists under every tested SWC definition. However, individual indirect pathways are definition-sensitive. Relative to the fuzzy baseline, indirect-effect direction agreement was 83.3% for past-only measured SWC, 87.5% for daily exact-date GAM SWC, 66.7% for the 7-day antecedent mean, and 41.7% for the 14-day antecedent mean. Exact-date measured comparisons had 100% direction agreement but covered only four models.

## Interpretation and recommendation

The sensitivity analyses support two conclusions at different strengths:

1. The randomized-treatment conclusions and the broad pattern that most estimated component magnitude remains in the SWC-adjusted treatment coefficient are reasonably robust.
2. The exact size, sign, and significance of individual SWC-associated indirect effects are not robust enough for strong causal mediation claims.

For a concise physiological paper, the cleanest approach is to treat the SEM as an associative decomposition rather than proof of mechanism. Describe `c'` as the **SWC-adjusted remaining treatment association** and `a*b` as the **SWC-associated indirect component**. If one measured-SWC specification must be preferred, use same-day or preceding SWC within seven days because it respects temporal ordering, and present the daily and antecedent GAM definitions as sensitivity analyses. Avoid interpreting extreme-event mediation strongly because the event variable and SWC availability are both date-dependent.

## Reproducible outputs

- `exact_date/`: exact-date measured-SWC script, model status, bootstrap effects, comparison, and PDF.
- `past_only_7d/`: past-only measured-SWC script, lag and sample diagnostics, bootstrap effects, comparison, and PDF.
- `daily_interpolation/`: interpolation audit, leave-date-out validation, bootstrap effects, comparison, and PDFs.
- `antecedent_windows/`: 7- and 14-day trailing-window script, coverage and reduced-form QA, bootstrap effects, comparison, and PDF.
- `scripts/01_compare_swc_matching_methods.R`: cross-method synthesis script.
- `output/swc-matching-sensitivity-summary.pdf`: cross-method coefficient-change heatmap.
- `output/swc-matching-comparison-summary.csv`: agreement and effect-change metrics by component.
- `output/swc-matching-pathway-dominance.csv`: direct/indirect magnitude comparison by method.
