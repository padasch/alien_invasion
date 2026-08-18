# Repeated-response SEM: delta/Wald versus container bootstrap

## Overall assessment

**Share with caveats.** The bootstrap preserves every original Figure 6 point
estimate and all bootstrap-mean effect directions. It changes uncertainty
enough to alter 7 of the 48 path-summed treatment conclusions. These changes
are scientifically relevant, but the high frequency of singular outcome-model
fits in some bootstrap samples and the 1,000-replicate pilot resolution should
be addressed before treating the new significance calls as final.

## Scope and implementation checks

- The inventory contains all 14 non-phenology Figure 6 species-response slots.
- Twelve slots were estimable. Quercus `remaining_green` and `chlavg` had zero
  rows after the current SEM filtering and remain unavailable, as in Figure 6.
- Each estimable slot reached 1,000 successful replicates in exactly 1,000
  attempts: 12,000 complete container resamples and 24,000 fitted LMM
  submodels. No fitting attempt failed.
- Runtime was 375.7 seconds using four parallel workers.
- Each bootstrap reused the original standardized data, resampled whole
  containers within block, assigned duplicate container/tree IDs, and refitted
  the frozen original mediator and outcome formulas.
- The 48 exported path-summed point effects agree exactly with the current
  Figure 6 matrix sources (maximum absolute difference = 0).
- All 192 effect-level bootstrap distributions contain exactly 1,000 estimates.
- Bootstrap-mean and original point-estimate directions agree for all 192
  treatment/path combinations. The median absolute difference between the
  original estimate and bootstrap mean was 0.0015 (maximum 0.0140).

## Inference comparison

Among the 144 direct, indirect, and total effects, 11 significance calls
changed under the empirical bootstrap probability. One additional
treatment-to-SWC path changed.

| Component | Became significant | Became non-significant | Median bootstrap/old interval width |
|---|---:|---:|---:|
| Direct | 1 | 0 | 0.95 |
| Indirect (`a*b`) | 0 | 3 | 1.11 |
| Path-summed total | 6 | 1 | 0.88 |
| Treatment to SWC | 1 | 0 | not summarized |

The bootstrap therefore tended to widen product-path uncertainty while
narrowing total-effect uncertainty. This is plausible because the old delta
calculation did not include covariance between the direct and indirect paths.
For example, bootstrap direct and indirect effects were strongly negatively
correlated for the Fagus volume extreme-event total (`r = -0.986`) and the
Quercus volume extreme-event total (`r = -0.962`), substantially reducing the
uncertainty of their sums.

### Changed path-summed totals

| Species | Response | Treatment | Estimate | Old p | Bootstrap p | Bootstrap 95% interval | Change |
|---|---|---|---:|---:|---:|---:|---|
| Fagus | Vitality | Extreme event | -0.073 | 0.022 | 0.058 | [-0.143, 0.001] | lost |
| Fagus | Volume (total) | Extreme event | 0.078 | 0.063 | 0.002 | [0.053, 0.103] | gained |
| Quercus | Chlorophyll | Extreme event | -0.125 | 0.060 | 0.004 | [-0.229, -0.040] | gained |
| Quercus | Vitality | Culture | 0.135 | 0.287 | 0.044 | [0.002, 0.358] | gained |
| Quercus | Quantum yield | Precipitation | 0.185 | 0.084 | 0.028 | [0.028, 0.372] | gained |
| Quercus | Volume (total) | Extreme event | -0.061 | 0.244 | 0.002 | [-0.080, -0.043] | gained |
| Quercus | Volume (increment) | Precipitation | -0.168 | 0.053 | 0.032 | [-0.308, -0.016] | gained |

The three indirect effects that became non-significant were precipitation,
Robinia, and extreme-event paths for Fagus relative phase volume increment;
each had bootstrap `p = 0.058` and a percentile interval that narrowly crossed
zero. The one changed direct effect was the Quercus vitality culture path,
which became significant (`p = 0.044`). The Quercus vitality culture effect on
SWC also became significant (`p = 0.030`).

## QA caveats

1. **Singular bootstrap outcome fits.** The SWC mediator models had no singular
   fits, but 2,504 of 12,000 outcome-model refits were singular. Counts were
   concentrated in Quercus volume (705), Fagus relative phase volume increment
   (449), Fagus vitality (429), and Fagus volume (402). All original-data fits
   were non-singular, and all retained bootstrap fits produced finite path
   estimates. Boundary estimates are possible under cluster resampling, so they
   were retained and explicitly counted rather than silently discarded. A
   sensitivity analysis should assess whether the changed conclusions persist
   with a simpler justified random-effects structure or another cluster-aware
   bootstrap implementation.
2. **Pilot precision.** With 1,000 replicates, the smallest possible adjusted
   two-sided probability is approximately 0.002. Use at least 5,000 successful
   replicates for final manuscript inference, particularly for effects near
   0.05.
3. **One interval/p-value boundary mismatch.** The Quercus quantum-yield culture
   direct effect had a percentile interval barely above zero
   `[0.00006, 0.326]` but empirical `p = 0.052`. This can occur because the
   percentile interval and the specified finite-sample tail probability are
   not exact inverses. Figure integration should consistently use the agreed
   inferential rule and show intervals transparently.
4. **Extreme-event scope.** Container resampling conditions on the observed
   measurement dates and extreme-event schedule. It estimates uncertainty
   across experimental containers, not across hypothetical alternative extreme
   events or event timings.
5. **Missing Quercus senescence cells.** No bootstrap inference exists for the
   two unavailable species-response combinations; they must remain visibly
   missing rather than be interpreted as null effects.

## Recommendation

Use this pilot to recreate Figure 6 and compare the visual story, but do not yet
replace the manuscript inference solely on the new p-values. First rerun the
final candidate with at least 5,000 replicates and resolve the sensitivity of
the seven changed totals to singular outcome fits. The most robust current
conclusion is that effect directions and magnitudes are stable, whereas several
significance classifications depend on the uncertainty method.
