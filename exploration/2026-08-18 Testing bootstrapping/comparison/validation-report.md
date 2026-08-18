# Validation of block-stratified container bootstrap inference

## Bottom line

The bootstrap results **agree with the current analysis on effect direction and the main biological story**, but they do not reproduce every threshold-based significance call. No compared point estimate reversed sign. The strongest agreement occurs for harvest biomass and the primary phenology model. Differences concentrate in isolated date-specific temporal contrasts and in the older repeated-response SEM, where uncertainty for a sum of direct and indirect paths is not well represented by the previous delta/Wald approximation.

This is a 1,000-replicate sensitivity analysis, not yet a recommendation to replace every manuscript interval. The container bootstrap is better aligned with treatment assignment and repeated measurements, but the repeated-response SEMs and several mixed models frequently become singular in resampled data. The final manuscript analysis should therefore use a single predeclared uncertainty framework, increase the final bootstrap run to at least 5,000 successful replicates, and test whether the SEM changes persist under a scientifically defensible simplified random-effects structure.

## What was held constant

Across families, the original analysis data, response scaling, estimands, and fixed/random-effect formulas were retained. Whole experimental containers were sampled with replacement within block, and all nested trees, dates, transitions, and aligned soil-water observations were carried together. Duplicated containers and trees received synthetic identifiers. Existing repeated-response SEM formulas were frozen; model selection was not repeated inside bootstrap samples.

Intervals are the 2.5th and 97.5th percentiles of 1,000 successful estimates. The reported empirical two-sided probability is based on the fraction of estimates at or beyond zero, with a finite-sample correction. Consequently, differences from the current analysis reflect the uncertainty method rather than a change in the fitted biological question.

## Agreement by model family

| Model family | Successful resamples | Compared effects | Direction changes | 95% interval-classification changes | Assessment |
|---|---:|---:|---:|---:|---|
| Temporal LMMs | 4,000 | 120 | 0 | 10 | Same trajectories; five intervals newly exclude zero and five cease to exclude zero |
| Harvest-biomass LMMs | 6,000 | 18 | 0 | 0 | Complete agreement |
| Primary phenology common-shift LMM | 2,000 | 6 | 0 | 0 | Complete agreement; all six effects remain non-significant |
| Phenology duration sensitivity | 2,000 | 6 | 0 | 0 | Complete agreement |
| Phenology stage-specific sensitivity | 2,000 | 18 | 0 | 1 | One borderline Fagus stage-4 contrast changes classification |
| Phenology SEM components | 2,000 | 24 | 0 | 1 | One borderline direct path changes; all six Figure 6 totals are unchanged |
| Repeated-response SEM totals | 12,000 | 48 | 0 | 7 | Magnitudes/directions stable; uncertainty-sensitive significance calls |

The successful-resample counts for the three phenology rows refer to the same 2,000 species-level bootstrap samples, each of which refitted the relevant primary, interaction, and duration models. They should not be added as independent resamples.

## Scientifically relevant differences

### Temporal responses

Ten of 120 pointwise treatment contrasts changed whether their 95% interval excluded zero. These changes were balanced: five newly excluded zero and five no longer did. None changed direction, and the dominant seasonal trajectories remain intact. Because the 120 intervals are pointwise and unadjusted for repeated testing across dates and treatments, these isolated threshold changes should be described as sensitivity results rather than new temporal events. The early Fagus diameter contrast is especially fragile because baseline diameter increment is constructed near zero.

### Phenology

The primary model of overall phenological timing is fully stable: all six treatment directions and all six non-significant classifications are unchanged. The stage 2–4 duration analysis is also unchanged. Two secondary results sit at the 0.05 boundary:

- Fagus mixed culture at stage 4 changes from Kenward–Roger P = 0.0706 to bootstrap P = 0.0500.
- The Quercus reduced-precipitation direct path in the phenology SEM changes from the earlier unstratified P = 0.0659 to block-stratified bootstrap P = 0.0480.

Neither change alters the paper's intended overall-timing conclusion, and all six path-summed phenology effects used in Figure 6 remain non-significant. The common-shift LMM and phenology SEM also agree on the displayed direction of all six treatment effects (six of six), despite answering different questions: the LMM estimates the total treatment association with timing, whereas the SEM decomposes that association into SWC-associated and remaining paths.

### Harvest biomass

All 18 coefficients retain both their direction and their interval classification. Ten effects exclude zero under both approaches. This is the clearest evidence that the biomass conclusions are not artifacts of the original Wald intervals.

### Repeated-response SEM

All 48 path-summed Figure 6 point estimates match the current sources exactly, and all bootstrap distributions retain the original direction. Seven total-effect classifications change:

- Fagus vitality under the extreme event ceases to exclude zero.
- Fagus total volume under the extreme event newly excludes zero.
- Quercus chlorophyll under the extreme event newly excludes zero.
- Quercus vitality under mixed culture newly excludes zero.
- Quercus quantum yield under reduced precipitation newly excludes zero.
- Quercus total volume under the extreme event newly excludes zero.
- Quercus relative volume increment under reduced precipitation newly excludes zero.

These are not contradictory point estimates. They arise because bootstrap uncertainty includes the covariance between each direct path and its indirect product path. Strong negative covariance can make the uncertainty of their sum smaller even when uncertainty in the indirect product alone is larger. That is a substantive advantage of joint resampling for the SEM totals.

## Diagnostics and limitations

- Temporal models had 16 rejected nonconvergent attempts before reaching 4,000 successful samples. Among accepted fits, 345 were singular, mainly Quercus diameter increment.
- Five of six biomass point models were singular, and 4,459 of 6,000 resampled fits estimated the container variance at or near zero. This indicates little remaining between-container variance after the fixed effects, not failure to calculate the treatment contrasts.
- Phenology produced no failed attempts, but 416 of 2,000 primary-LMM fits and 410 of 2,000 SEM outcome fits were singular, mainly for Fagus.
- Repeated-response SEM mediator models were never singular, whereas 2,504 of 12,000 outcome models were singular. The concentration in Quercus total volume and several Fagus responses makes the seven changed total effects provisional.
- Two Quercus senescence response slots contain no analysable SEM rows and remain explicitly missing, not null effects.
- A 1,000-replicate run has coarse tail resolution near P = 0.05 and a minimum adjusted two-sided probability of about 0.002.

## Recommendation

The design-aligned container bootstrap is a strong candidate for a common manuscript uncertainty framework, especially for SEM path sums. Before replacing the current intervals, however:

1. rerun the candidate analysis with at least 5,000 successful replicates per estimable model;
2. predefine whether inference is based on percentile intervals, empirical probabilities, or both, because these are not exact finite-sample inverses;
3. evaluate a justified simpler random-effects structure for model families with frequent singular bootstrap fits;
4. keep temporal inference explicitly pointwise or add a multiplicity strategy if individual dates are emphasized; and
5. report effect sizes and intervals as primary evidence, treating P = 0.05 crossings as sensitivity rather than biological reversals.

On the current 1,000-replicate evidence, the manuscript's broad treatment-response story is robust. The elements that need further work are a small set of borderline temporal contrasts and seven repeated-response SEM totals—not the estimated direction or magnitude of the main effects.

## Reproducibility map

- Shared conventions: `../README.md`
- Temporal models and comparison: `../temporal_lmm/`
- Harvest biomass and comparison: `../biomass_lmm/`
- Phenology models and comparison: `../phenology/`
- Repeated-response SEMs and comparison: `../repeated_response_sem/`
- Recreated Figure 2–6 and phenology supplementary figures: `../recreated_figures/`
- Machine-readable family summary: `family-summary.csv`
