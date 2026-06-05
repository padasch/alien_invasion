# Methods

## Data analysis

### Data preparation

To quantify treatment effects on tree growth, we used height and diameter measurements to calculate a simple cylindrical stem-volume proxy,

$$
V = \pi \left(\frac{d}{2}\right)^2 h,
$$

where $d$ is stem diameter (mm) and $h$ is height (mm). For all three metrics, we further derived the absolute increment relative to the first measurement, the relative increment relative to the first measurement, the absolute increment within phase, and the relative increment within phase. Phase-specific increments were referenced to the first measurement in phase 1 and to the last measurement of the previous phase in later phases.

Leaf senescence was retained in its measured form as senesced canopy area (%) for shown timeseries. For analyses, however, we modelled remaining green canopy as "100 - senescence (%)", so that negative treatment effects represent earlier or faster senescence progression.


TODO:
To analyse phenology, we used two complementary approaches because treatment effects can be expressed either as developmental status at a given date or as shifts in the timing of specific stage transitions. The first data set retained the repeated stage observations for each tree and was used in the same temporal mixed-model framework described above, thereby testing whether trees were more or less developmentally advanced at a given date.

The second data set contained the cleaned transition dates at which stages 2-4 were reached. These transitions were modelled directly as day of year, first in separate stage-specific mixed-effects models and then in a pooled model across stages with stage as a fixed effect and tree identity retained as a random intercept. In these transition-timing models, negative coefficients indicate earlier attainment of a phenological stage and positive coefficients indicate delayed attainment.

- Any other wrangling steps needed to describe here?
    -   For example, minor cleaning steps, for example where a tree was noted to "have gone back" in phenology from (eg: 1 - 2 - 3 - 2 - 3 was cleaned into 1 - 2 - 3 - 3 - 3).

### Temporal mixed-effects models

To test how treatment effects changed through the season, we fitted species-specific mixed-effects models to repeated measurements of each response variable. The data were structured as repeated observations per tree, linked to box identity and sampling date. Before fitting, each response was standardized to a $z$-score within the analysed subset so that effect sizes were comparable across traits.

The temporal model for each trait $y$ was

$$
y = \mathrm{date} + \mathrm{date}:\mathrm{precipitation} + \mathrm{date}:\mathrm{Robinia} + \mathrm{date}:\mathrm{culture} + (1 | \mathrm{box}) + (1 | \mathrm{tree}).
$$

Date-specific treatment effects were extracted as treatment minus baseline contrasts, such that the plotted coefficients describe the change from control to drought, from without *Robinia* to with *Robinia*, from monoculture to mixture, and, where shown separately (Figures TODO). Negative coefficients therefore indicate reduced performance relative to the corresponding reference level.

Endpoint biomass traits were analysed separately from the repeated trajectories because they were measured only once at harvest. Shoot biomass, root biomass, and root:shoot ratio were fitted as species-specific mixed-effects models with precipitation, *Robinia*, and culture treatments as fixed effects and box identity as a random intercept.

### Structural equation modelling

To separate direct treatment effects from effects mediated through soil moisture, we fitted species- and trait-specific piecewise structural equation models (SEMs). The data were structured at the tree-observation level, with each response matched to the nearest measured SWC value within a $\pm 7$ d window, giving preference to the most recent preceding SWC observation. We also interpolated SWC measurements to daily conditions using climate and soil water potential sensor data, however results did not change substantially and we therefore kept using actual measured SWC data (see Supplementary Information, TOOD). All continuous variables used in the SEMs, including SWC, the focal response, and centred day-of-year terms, were standardized before fitting. We included both linear and quadratic day-of-year terms to account for non-linear seasonal trajectories in soil moisture and tree responses over the course of the growing season. Each SEM consisted of two mixed-effects submodels:

$$
\mathrm{SWC} = \mathrm{DOY} + \mathrm{DOY}^2 + \mathrm{precipitation} + \mathrm{Robinia} + \mathrm{soil} + \mathrm{culture} + \mathrm{extreme\ event} + (1 | \mathrm{box}),
$$

and

$$
y = \mathrm{SWC} + \mathrm{DOY} + \mathrm{DOY}^2 + \mathrm{precipitation} + \mathrm{Robinia} + \mathrm{soil} + \mathrm{culture} + \mathrm{extreme\ event} + (1 | \mathrm{box}) + (1 | \mathrm{tree}).
$$

For each treatment, we decomposed the standardized effect into a direct path ($c'$), an indirect path through SWC ($a \times b$), and a total effect ($c' + a \times b$). Direct-path uncertainty was taken from the fitted mixed-model coefficients, whereas indirect and total-effect uncertainty was propagated from the component path variances. To summarize these results across traits, we converted the SEM outputs into heatmaps of direct, indirect, and total effects, together with complementary matrices of "treatment $\rightarrow$ SWC" and "SWC $\rightarrow$ response" (Figures TODO).

- Probably more for the figure caption: "All heatmaps use a shared diverging colour scale so that effect magnitudes were directly comparable across responses. Cells were blanked when the corresponding effect was not significant at $P \geq 0.05$."

- We can show either total effect of treatments (which includes the SWC-pathway), the direct effect (not including SWC-pathway), or both.

### Sensitivity analyses

We carried out two sensitivity analyses to test whether the main conclusions depended on presentation choices. First, because SWC was not measured daily, we compared the measured SWC series against a treatment-agnostic daily interpolation fitted with a generalized additive model using day of year, soil water potential, air temperature, vapour pressure deficit, precipitation, radiation, and box identity. Re-fitting the SEMs with interpolated SWC did not materially alter the qualitative pattern of direct versus SWC-mediated treatment effects.

Second, the main figures pool both soil types and omit soil as a treatment factor to simplify interpretation. To test whether this masked important legacy effects, we repeated the analyses with soil type separated or modelled explicitly. These analyses showed that one soil type was consistently somewhat drier and amplified the negative effects of drought and *Robinia*, but they did not alter the overall direction of the main treatment responses. These soil-specific results are therefore best presented as supplementary sensitivity figures (Supplementary Fig. Sx-Sy).

Note that data was pooled across both soil types due increase sample size for analysis. In the Supporting Information, we present equivalent results as described below, using soil type as another treatment factor alongside the others. Overall, soil type had a comparatively small influence, primarily through one being sandier and less water-retaining than the other, which intensified drought-stress induced by the precipitation treatment and robinia presence (TODO, not sure if we have the data to proof the soil texture difference).


### Software

All analyses were performed in R using packages from the `tidyverse`; mixed-effects models were fitted with `lme4`, temporal contrasts were extracted with `emmeans`, piecewise SEMs were assembled with `piecewiseSEM`, and SWC interpolation was fitted with `mgcv`.

## Notes for manuscript placement

- Some of the variable-derivation details in the `Data preparation` section, especially the growth-derived variables and the remaining-green transformation, could be shortened further by moving them into the corresponding measurement subsections.
- Supplementary figure placeholders (`Fig. Sx-Sy`) should be replaced once the soil-separated and SWC-sensitivity panels are finalized.
