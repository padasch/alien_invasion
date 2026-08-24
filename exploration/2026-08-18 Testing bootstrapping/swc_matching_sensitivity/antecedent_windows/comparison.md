# Antecedent modeled-SWC sensitivity

## Design

Repeated-response SEMs were refitted using modeled trailing SWC means ending on each response date. The 7-day definition uses response day through day -6; the 14-day definition uses response day through day -13. The daily latent input is the existing climate-informed, treatment-agnostic container-level GAM estimate (`swc_hat`). Neither definition uses future information. These mediators are modeled antecedent means, not measured means.

Frozen original formulas and model decisions were retained. Each model retained its fuzzy-baseline response observations, response scaling, and original SWC mean/SD, so comparisons isolate the mediator definition. Uncertainty used 1,000 successful block-stratified container-cluster bootstrap replicates per estimable species-response model and window.

## Coverage and fitting

All estimable baseline response records had complete 7- and 14-day modeled windows (minimum model-level coverage 100.0%). The two Quercus senescence models remained unavailable because their original SEM bundles contained no post-filtering response data. All 24 estimable model-window tasks reached 1,000 successful replicates in 1,000 attempts; there were no failed refits.

The SWC submodels were never singular. Response-model singularity was frequent for several endpoints: 69.1% of Quercus total-volume replicates under the 7-day definition, 64.7% under 14 days, 60.6% for Fagus incremental volume under 14 days, and approximately half of Fagus vitality replicates. These are retained, as in the baseline bootstrap, but they indicate weak or boundary random-effect variance and require cautious interval interpretation.

## Pathway comparison

**Path-summed totals were effectively unchanged, but this is largely algebraic.** Relative to fuzzy matching, total-effect direction and interval overlap were both 100% for both antecedent windows. Median absolute changes were only 2e-04 (7 days) and 2e-04 (14 days), with maxima of 0.0089 and 0.0095 standardized units.

**The direct/indirect allocation was sensitive.** Direct-effect direction agreed with the fuzzy baseline for 93.8% of paths under 7 days and 91.7% under 14 days; corresponding interval overlap was 95.8% and 93.8%. Indirect-effect direction agreement was only 66.7% and 41.7%, and interval overlap was 58.3% and 68.8%.

The shared SWC-to-response path (`b`) was especially unstable: its direction agreed with the fuzzy baseline for 66.7% of the 12 response models under 7 days and only 41.7% under 14 days. Significant `b` paths numbered 8 for fuzzy matching, 9 for 7 days, and 5 for 14 days. Significant indirect paths numbered 23, 29, and 17, respectively.

Examples show why a single mechanistic conclusion is unsafe: the oriented Fagus chlorophyll `b` path changed from -0.198 under fuzzy matching to -0.074 under 7 days and +0.055 under 14 days; Quercus quantum yield changed from +0.076 to -0.083 and -0.054. These sign changes are conditional associations after DOY and treatment adjustment, not necessarily biological reversals.

The broad magnitude statement that most modeled treatment effects were not allocated to the SWC-mediated path was nevertheless retained: median absolute indirect shares were 14.1% (fuzzy), 19.6% (7 days), and 7.3% (14 days). Indirect magnitude exceeded direct magnitude in 6, 2, and 2 of 48 paths. This supports only a broad descriptive statement, not stable attribution for individual treatment-response combinations.

## Regression-algebra check

Stable path-summed totals are not independent evidence of robustness. In linear path models using the same observations and covariates, `c' + a*b` can be nearly invariant to mediator redefinition by regression algebra. The point path sum differed from the corresponding reduced-form treatment coefficient (response model with SWC omitted) by a median of 1e-04 standardized units; the 95th percentile was 0.0035 and the maximum was 0.0109. Mixed-model random-effects estimation prevents exact equality, but the near-identity confirms that decomposition sensitivity, not total-effect stability, is the meaningful test.

## Scientific recommendation

Do not select a window because it produces more significant indirect paths. For a simple, pre-specified sensitivity analysis, the prior 7-day mean is the more defensible common window for rapidly responsive physiological traits; retain the 14-day mean as a longer-exposure check. Growth and senescence could plausibly integrate water availability over longer or stage-specific periods, but adopting endpoint-specific windows post hoc would add researcher degrees of freedom and should require an a priori biological rationale.

More importantly, neither modeled window should be treated as definitive causal mediation. The daily GAM is climate-informed but treatment-agnostic, and its prediction uncertainty was not propagated. Ambient climate describes common temporal forcing but cannot reconstruct container-specific experimental watering or exclusion without dated treatment-operation records. Until a treatment-aware daily SWC reconstruction is available, report the total treatment patterns separately and describe direct-versus-SWC-mediated allocation as exploratory and sensitive to mediator timing.

## QA

All estimable point and bootstrap SEMs retained exact `total = direct + indirect`; the maximum absolute numerical error was 1.110223e-16.
