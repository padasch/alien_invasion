# Methods

## Data analysis

Analyses were conducted separately for *Fagus sylvatica* and *Quercus ilex*. The primary analyses included observations from both soil-inoculation treatments but did not include soil type as a predictor. Precipitation treatment (control or reduced precipitation), *Robinia* treatment (without or with *R. pseudoacacia*), and culture type (monoculture or mixed culture) were coded as categorical predictors, with control precipitation, absence of *Robinia*, and monoculture as the respective reference levels.

### Derived response variables

Height and basal diameter were combined into a species-specific allometric growth proxy following Annighöfer et al. (2016):

$$
G = \beta_1\left(D^2H\right)^{\beta_2},
$$

where $D$ is basal diameter (mm), $H$ is height (cm), and $G$ is the estimated allometric proxy (g). We used the published coefficients for *F. sylvatica* ($\beta_1 = 0.62342$, $\beta_2 = 0.87409$). Because an equation for *Q. ilex* was unavailable in the selected source, the *Q. robur* coefficients ($\beta_1 = 0.67311$, $\beta_2 = 0.85202$) were used as a surrogate. The analysis code and current figures refer to $G$ as a volume proxy for continuity, although it should be interpreted as an allometric biomass or size proxy rather than geometric stem volume.

For height, diameter, and the allometric proxy, we calculated absolute and relative increments from the first measurement. We additionally calculated phase-specific increments for spring (up to and including June), summer (July-August), and autumn (September onward). The spring baseline was the first measurement, whereas summer and autumn increments were referenced to the final measurement of the preceding phase. Phase-specific trajectories therefore restart at zero at each phase boundary.

Leaf-senescence analyses used remaining green canopy, calculated as $100 -$ visually estimated senesced canopy area (%), so that negative coefficients denote less remaining green canopy and hence earlier or faster senescence. Vitality was analysed in the orientation used by the analysis data, in which higher values denote better condition; this orientation should be kept consistent with the measurement-scale description elsewhere in the manuscript.

For spring phenology, non-monotonic stage sequences were corrected before determining the first date on which each tree reached stages 2–4. Transition dates were expressed as day of year (DOY). On this natural scale, negative treatment coefficients indicate earlier and positive coefficients later leaf development. For presentation alongside the other response variables, the signs of phenology effects were reversed so that negative displayed values indicate delayed or slower development; this reorientation did not affect model fitting or inference.

### Temporal mixed-effects models

Seasonal treatment effects on repeated responses were estimated with species-specific Gaussian linear mixed-effects models. Each response was standardized to a $z$ score within the analysed species-trait subset. Sampling date was treated as a categorical fixed effect and was allowed to modify the effect of each treatment:

$$
y_z = \mathrm{date} + \mathrm{date}:\mathrm{precipitation}
+ \mathrm{date}:\mathrm{Robinia}
+ \mathrm{date}:\mathrm{culture}
+ (1\mid\mathrm{container}) + (1\mid\mathrm{tree}).
$$

Container and tree identity were included as random intercepts to account for the experimental-unit and repeated-measure structure, respectively. Models were fitted by restricted maximum likelihood. For each date, estimated marginal contrasts were calculated as treatment minus reference using `emmeans`: reduced minus control precipitation, with minus without *Robinia*, and mixed culture minus monoculture. Contrasts are reported with pointwise 95% asymptotic Wald confidence intervals. No multiplicity adjustment was applied to the date-specific contrasts. These models were additive across treatments and did not test precipitation × *Robinia*, culture × *Robinia*, or other treatment-by-treatment interactions.

### Phenological transition models

For each species, the DOY at which trees reached stages 2–4 was analysed jointly using a Gaussian linear mixed-effects model:

$$
\mathrm{DOY} = \mathrm{stage} + \mathrm{precipitation} + \mathrm{Robinia}
+ \mathrm{culture} + \mathrm{block}
+ (1\mid\mathrm{container}) + (1\mid\mathrm{tree}).
$$

Phenological stage, precipitation treatment, *Robinia* treatment, culture type, and experimental block were included as fixed effects. Container and tree identity were included as random intercepts to account for the experimental-unit structure and the repeated transition dates from each tree, respectively. The primary model assumes that each treatment produces a common shift in the timing of stages 2–4; the stage effect accounts for their different mean dates. Models were fitted by restricted maximum likelihood, and treatment contrasts were estimated in days with 95% confidence intervals using Kenward–Roger degrees of freedom.

As a sensitivity analysis, interactions between phenological stage and each treatment were added to determine whether treatment effects varied among stages. Stage-specific contrasts were estimated from this model, and interaction terms were evaluated by likelihood-ratio comparisons of models fitted by maximum likelihood. As a secondary measure of developmental rate, the elapsed time from stage 2 to stage 4 was calculated for trees with both transition dates and analysed with an additive Gaussian mixed-effects model containing precipitation, *Robinia*, culture, and block as fixed effects and container as a random intercept. Positive duration coefficients denote a longer, and hence slower, progression from stage 2 to stage 4.

### Harvest biomass

Shoot biomass, root biomass, and root:shoot ratio were analysed separately for each species. Each response was standardized within species and fitted with a Gaussian linear mixed-effects model containing the additive main effects of precipitation, *Robinia*, and culture and a random intercept for container:

$$
y_z = \mathrm{precipitation} + \mathrm{Robinia} + \mathrm{culture}
+ (1\mid\mathrm{container}).
$$

Models were fitted by restricted maximum likelihood. Fixed-effect estimates are presented with 95% Wald confidence intervals calculated as the estimate $\pm 1.96$ standard errors. The biomass models did not include treatment interactions.

### Piecewise structural equation models

We used species- and trait-specific piecewise structural equation models (SEMs) to distinguish treatment-response associations represented as direct paths from associations statistically transmitted through soil water content (SWC). For the repeatedly measured physiological responses, each tree observation was matched to the most recent SWC measurement from the same container within the preceding seven days; when none was available, the earliest following observation within seven days was used. All observations retained in these SEMs had an SWC match within this $\pm7$-day window. Spring phenology was analysed using the separate procedure described below.

The extreme-event indicator identified measurements made during either summer drought event or its 14-day post-event window. In the implemented analysis, these windows were 20 June-17 July and 12 August-3 September 2025. The indicator is consequently a time-window contrast rather than an experimentally randomized treatment.

The response, SWC, and centred DOY were standardized before model fitting. Linear and quadratic DOY terms were included to represent nonlinear seasonal trajectories. For each response, the SEM consisted of an SWC submodel and a response submodel:

$$
\mathrm{SWC}_z = \mathrm{DOY}_z + \mathrm{DOY}_z^2
+ \mathrm{precipitation} + \mathrm{Robinia} + \mathrm{culture}
+ \mathrm{event} + (1\mid\mathrm{container}),
$$

$$
y_z = \mathrm{SWC}_z + \mathrm{DOY}_z + \mathrm{DOY}_z^2
+ \mathrm{precipitation} + \mathrm{Robinia} + \mathrm{culture}
+ \mathrm{event} + (1\mid\mathrm{container}) + (1\mid\mathrm{tree}).
$$

The repeatedly measured response SEMs were additive and did not include treatment interactions. Candidate treatment terms were subjected separately in the two submodels to backward AIC selection fitted by maximum likelihood. A term was removed only when doing so reduced AIC by at least two units; linear and quadratic DOY were retained in both submodels, and SWC was retained in the response submodel. Selected models were then refitted by restricted maximum likelihood.

For the repeatedly measured responses and treatment $x$, the direct association was the coefficient $c'$ from the response submodel, the SWC-mediated component was the product $ab$, where $a$ is the $x\rightarrow$ SWC coefficient and $b$ is the SWC $\rightarrow y$ coefficient, and the total association was $c' + ab$. Approximate standard errors were obtained using first-order variance propagation, $\mathrm{Var}(ab)=b^2\mathrm{Var}(a)+a^2\mathrm{Var}(b)$ and $\mathrm{Var}(c'+ab)=\mathrm{Var}(c')+\mathrm{Var}(ab)$. Two-sided Wald tests used the standard normal distribution. Because SWC and the responses were observationally aligned in time, these path decompositions describe statistical associations and do not alone establish causal mediation.

The phenology SEM used a stage-centred overall timing index. Within each species, transition DOY was centred on the mean date of its respective stage, and each tree's mean deviation was calculated from the available transitions; trees were retained when at least two of stages 2–4 had been observed. The index was standardized within species and sign-reversed for presentation, such that negative effects denote delayed overall phenological timing.

To ensure that SWC preceded the phenological transitions, pre-leaf-out SWC was calculated as the mean measured SWC in each container during a fixed window from 4 March to 2 April 2025, which ended before the earliest focal transition. The SWC submodel was fitted to one observation per container and included precipitation, *Robinia*, culture, and block as additive fixed effects. The phenology submodel was fitted at tree level and included standardized SWC and the same treatment and block terms as fixed effects, with a random intercept for container. Direct, indirect ($ab$), and path-summed total ($c' + ab$) associations were calculated without model selection. Confidence intervals and two-sided probabilities were obtained from 1,000 bootstrap samples in which containers, together with their constituent trees, were resampled as clusters. These paths describe conditional associations and should not be interpreted as evidence of causal mediation.

### Statistical interpretation and software

Unless stated otherwise, treatment estimates are expressed in within-species standard-deviation units. Coefficients standardized separately by species should therefore not be interpreted as direct tests of between-species differences. Statistical evidence was evaluated from effect estimates, 95% confidence intervals, and two-sided $P$ values, with $P<0.05$ used for the significance labels in the figures. The primary analyses did not adjust for the number of traits, dates, or SEM paths and should be interpreted together with effect magnitude, consistency, and sensitivity analyses.

Analyses were conducted in R. Mixed-effects models were fitted with `lme4`, marginal contrasts were obtained with `emmeans`, piecewise SEM objects were assembled with `piecewiseSEM`, and data processing and figures used packages from the `tidyverse`, including `dplyr`, `tidyr`, and `ggplot2`.
