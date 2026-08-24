# Exact-date daily-SWC sensitivity

## Design

Repeated-response SEMs were rebuilt on exact response dates using observed SWC when measured on that date and the existing container-level daily GAM prediction otherwise. Frozen formulas, response scaling, and each baseline model's original SWC mean/SD were retained. Uncertainty used 1,000 successful block-stratified container-cluster bootstrap replicates per estimable model.

## Interpolation audit

The existing GAM is climate-informed but treatment-agnostic: ambient precipitation, air temperature, VPD, radiation and site-level soil-water potential vary through time, while a container random intercept represents persistent container differences. It has no precipitation-treatment input or treatment-specific temporal smooth. Leave-measurement-date-out validation, with predictions constrained to the physical SWC range as in the production export, gave RMSE 5.42, MAE 4.15, bias 0.01, and predictive R-squared 0.40.

## SEM comparison

Across 144 comparable direct, indirect and total paths, effect direction agreed for 135 (93.8%), bootstrap significance agreed for 129 (89.6%), and intervals overlapped for 135 (93.8%).

The net total effects were exceptionally stable: all 48 retained their direction, 47/48 retained the same bootstrap significance classification, all intervals overlapped, and the mean absolute estimate change was 0.001 SD. The decomposition was less stable. Direct-effect direction agreed for 45/48 paths and indirect-effect direction for 42/48; indirect-effect significance agreed for 37/48. Under exact-date daily SWC, the direct component remained larger in absolute magnitude than the indirect component for 44/48 treatment-response combinations. Thus, the net treatment conclusions are robust to this SWC definition, whereas the quantitative attribution between direct and SWC-associated pathways is definition-sensitive.

## Climate-covariate recommendation

Do not add ambient climate as a substitute for manipulated precipitation: the current daily GAM already uses ambient climate, which is common to all containers and cannot reconstruct container-specific water exclusion/addition. Without dated container-level irrigation or throughfall-exclusion inputs, a more complex climate model would add assumptions rather than treatment information. Treat this exact-date GAM analysis as a sensitivity check. If pathway decomposition changes materially, prefer same-day/prior measured SWC or a transparent treatment-aware interpolation supported by actual watering records. Bootstrap intervals here are conditional on the interpolated series and do not propagate GAM prediction uncertainty.

## QA

All estimable SEMs retained exact total = direct + indirect identity (maximum absolute numerical error 1.110223e-16).
