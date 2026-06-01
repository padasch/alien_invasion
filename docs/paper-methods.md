# Methods

## Data preparation

Raw spreadsheet- and sensor-based measurements were cleaned into structured interim tables for experimental metadata, growth, chlorophyll, quantum yield, canopy condition, senescence, phenology, soil water content, soil respiration, and isotope chemistry. Treatment factors were harmonized before analysis and coded with the following reference levels: control precipitation, monoculture, absence of *Robinia*, and the legacy *Robinia* soil as the wetter soil reference. In analyses where soil type was retained as an explicit factor, the beech soil represented the drier comparison level.

Repeated growth measurements were organized as longitudinal records of height and stem diameter for each tree. In addition to absolute height and diameter, we derived a cylindrical stem-volume proxy as

$$
V = \pi \left(\frac{d}{20}\right)^2 h,
$$

where $d$ is stem diameter in mm and $h$ is height in cm, yielding an approximate stem volume in cm$^3$. For height, diameter, and volume, we calculated five response forms used in the analyses: the absolute trait value, the absolute increment relative to the first measurement, the relative increment relative to the first measurement, the absolute increment within phase, and the relative increment within phase. Phase-specific increments were referenced to the first measurement in phase 1 and to the last measurement of the previous phase in subsequent phases. Artificial reset-to-zero points were added only for visualization of within-phase trajectories and were never used for model fitting.

Leaf-condition scores were transformed so that larger values represented better canopy condition. Senescence records were retained in their measured form as percentage senesced canopy area for descriptive plots, but inferential models used a derived variable, remaining green canopy, calculated as $100 -$ senescence (%), so that negative effects corresponded to earlier or faster senescence progression. For senescence-associated chlorophyll, two leaf-level chlorophyll measurements were averaged for each observation date.

Phenology was cleaned from the recorded dates of stages 1-4. Transition dates were rebuilt directly from the raw stage dates, day of year was recomputed from the cleaned calendar dates, and conflicting stage sequences were resolved so that stage timing remained monotonic within each tree. We retained two phenology views: a descriptive transition table giving the date and day of year at which each stage was reached, and a longitudinal progression table describing the stage observed at each measurement date.

For all tree-level analyses, observations were linked to the experimental metadata by tree identity and box identity. Measurements were also assigned to two predefined extreme-event periods, with a short post-event buffer retained to capture delayed responses in later measurements. The main manuscript presentation pooled both legacy soil types and omitted soil as an explicit treatment factor for simplicity, but soil-specific analyses were retained as a sensitivity check (see Supplementary Fig. Sx-Sy).

## Temporal mixed-effects models

Time-varying treatment effects were estimated separately for each species and response variable using linear mixed-effects models fitted across all measurement dates. Response variables were standardized to $z$-scores within the analysed subset so that effect sizes were directly comparable among traits. Sampling date was treated as a categorical fixed effect, and treatment effects were estimated through date-by-treatment interactions. In the full specification, the fitted model was

$$
y_{ijkt} = \mu + \mathrm{Date}_t + \mathrm{Date}_t \times P_j + \mathrm{Date}_t \times R_j + \mathrm{Date}_t \times C_j + \mathrm{Date}_t \times S_j + b_{\mathrm{box}} + b_{\mathrm{tree}} + \varepsilon_{ijkt},
$$

where $y_{ijkt}$ is the standardized response, $P$ denotes precipitation treatment, $R$ denotes *Robinia* treatment, $C$ denotes culture treatment, $S$ denotes soil type, $b_{\mathrm{box}}$ is a random intercept for mesocosm box, and $b_{\mathrm{tree}}$ is a random intercept for tree identity. In the pooled-soil presentation used in the main text, both soil types were retained in the data set but soil type was omitted from the fixed treatment contrasts. Accordingly, the plotted temporal effects represent the contrasts control $\rightarrow$ drought, without *Robinia* $\rightarrow$ with *Robinia*, and monoculture $\rightarrow$ mixture; where soil was analysed explicitly, the contrast was wetter *Robinia* soil $\rightarrow$ drier beech soil.

Date-specific treatment effects were extracted as treatment minus baseline contrasts at each sampling date. The plotted values therefore represent the standardized change associated with the named treatment relative to its stated reference level at that date. For growth, chlorophyll, quantum yield, condition, and remaining green canopy, negative coefficients indicate reduced performance relative to the reference level and positive coefficients indicate improved performance. Because senescence was modelled as remaining green canopy rather than percentage senesced area, negative estimates in that analysis correspond to earlier or faster senescence progression.

Endpoint biomass traits were analysed separately from the temporal trajectories. Shoot biomass, root biomass, and root:shoot ratio were fitted as species-specific linear mixed-effects models on standardized responses, with precipitation, *Robinia*, and culture as fixed effects and box identity as a random intercept; soil type was included only in soil-separated sensitivity analyses. These endpoint models were used to summarize harvest-based allocation responses and were not included in the temporal effect trajectories.

## Structural equation models

To partition direct and soil-water-mediated treatment effects, we fitted piecewise structural equation models separately for each species and response variable. The mediator was soil water content (SWC) matched to each tree observation date from the nearest measured SWC record within a $\pm 7$ d window, with preference given to the most recent preceding observation. In sensitivity analyses, this measured SWC series was replaced by a treatment-agnostic daily interpolation (see below).

Each piecewise SEM consisted of two mixed-effects submodels. The SWC submodel was

$$
\mathrm{SWC}_{it} = \beta_0 + \beta_1 \mathrm{DOY}_{it} + \beta_2 \mathrm{DOY}_{it}^2 + \beta_3 P_i + \beta_4 R_i + \beta_5 S_i + \beta_6 C_i + \beta_7 E_t + b_{\mathrm{box}} + \varepsilon_{it},
$$

and the response submodel was

$$
Y_{it} = \gamma_0 + \gamma_1 \mathrm{SWC}_{it} + \gamma_2 \mathrm{DOY}_{it} + \gamma_3 \mathrm{DOY}_{it}^2 + \gamma_4 P_i + \gamma_5 R_i + \gamma_6 S_i + \gamma_7 C_i + \gamma_8 E_t + b_{\mathrm{box}} + b_{\mathrm{tree}} + \varepsilon_{it},
$$

where $E_t$ denotes the extreme-event indicator and day of year (DOY) was centred before inclusion as linear and quadratic terms. Continuous variables used in the SEMs, including SWC, the focal response, and centred day of year, were standardized before fitting. Both submodels were fitted with `lme4`, and the piecewise model was assembled with `piecewiseSEM`.

The main manuscript workflow used additive, direct-effect SEMs only; interaction SEMs were explored during analysis development but were not retained for the primary presentation. Treatment effects were decomposed into a direct path ($c'$), an indirect path through SWC ($a \times b$), and a total effect ($c' + a \times b$). Direct-path uncertainty was taken from the fitted mixed-model coefficients, whereas uncertainty for indirect and total effects was propagated from the variance of the component coefficients in the SWC and response submodels. Model performance was summarized from the two component mixed models using AIC, BIC, and marginal and conditional $R^2$, together with piecewise SEM fit summaries where available.

## Heatmap synthesis of SEM effects

To summarize mediation results across traits, we transformed the species-specific SEM outputs into standardized effect matrices. Separate heatmaps were generated for the direct effect of treatment on the response, the indirect effect mediated by SWC, and the total effect, and we additionally visualized treatment $\rightarrow$ SWC and SWC $\rightarrow$ response paths. All heatmaps within a scenario shared the same diverging colour scale so that effect magnitudes were directly comparable across panels and response variables. Warmer colours denote negative effects and cooler colours denote positive effects.

Heatmap cells were based on the standardized SEM coefficients described above. Cells were blanked when the corresponding path or decomposed effect was not significant at $P \geq 0.05$, so the displayed matrices show only statistically supported direct, indirect, and total effects. The heatmaps therefore provide a compact view of how drought, *Robinia*, culture, and, where analysed, soil legacy affected each response either directly or via changes in soil water availability.

## Phenology analyses

We used two complementary phenology analyses. First, the descriptive transition summaries quantified the mean day of year at which each treatment group reached phenological stages 1-4. Second, to estimate treatment effects on developmental timing more directly, we modelled the day of year at which stages 2-4 were attained. These transition-timing models were fitted separately for each species and stage as

$$
\mathrm{DOY}_{ij} = \alpha_0 + \alpha_1 P_i + \alpha_2 R_i + \alpha_3 C_i + \alpha_4 S_i + b_{\mathrm{box}} + \varepsilon_{ij},
$$

with soil type included only when soil-separated analyses were requested. We additionally fitted a pooled transition model across stages 2-4 with stage as a fixed effect and both box and tree identity as random intercepts, allowing average treatment effects on transition timing to be estimated while accounting for repeated stage observations within trees. In these transition-timing models, negative coefficients indicate earlier stage attainment and positive coefficients indicate delayed attainment.

In parallel, the main temporal GLMM and SEM framework retained a progression-model view of phenology, in which the response variable was the current phenological stage recorded at a given measurement date. This progression analysis answers a different biological question than the transition models: it tests whether a treatment group is more or less developmentally advanced at a given point in time, whereas the transition-timing models estimate whether stage attainment itself was advanced or delayed.

## Sensitivity analyses

We performed two sensitivity analyses to assess whether the main conclusions depended on modelling choices made for presentation. First, we compared SEM results obtained with measured SWC against SEMs re-fitted with a daily, treatment-agnostic SWC interpolation. The interpolation was generated with a generalized additive model (`mgcv`) driven by site-level hydrometeorological predictors, including day of year, soil water potential, air temperature, vapour pressure deficit, precipitation, radiation, and a random intercept for box identity. Because no treatment variables entered this interpolation model, the imputed SWC series did not encode treatment structure. Re-fitting the SEMs with interpolated SWC did not materially alter the qualitative pattern of direct versus indirect effects.

Second, although the main figures present both soil types together and omit soil type as a treatment factor for simplicity, we repeated the analyses with soil legacy separated or modelled explicitly. These complementary analyses showed that one soil type was consistently somewhat drier and that this drier legacy amplified the negative effects of drought and *Robinia* on several responses. Because these soil-type differences do not alter the direction of the main treatment conclusions but add complexity to the presentation, they are best shown as supplementary sensitivity figures (Supplementary Fig. Sx-Sy).

## Software

All data processing, modelling, and figure generation were conducted in R using packages from the `tidyverse`, with mixed-effects models fitted in `lme4`, temporal contrasts obtained with `emmeans`, piecewise structural equation models assembled with `piecewiseSEM`, and SWC interpolation fitted with `mgcv`.
