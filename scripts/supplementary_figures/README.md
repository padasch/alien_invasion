# Supplementary information

This directory contains reproducible scripts for the supplementary figures,
including phenology, temporal responses, biomass, and SEM decompositions. Run
the scripts from the repository root:

```sh
Rscript --vanilla scripts/supplementary_figures/fig_s16.R
Rscript --vanilla scripts/supplementary_figures/fig_s17.R
Rscript --vanilla scripts/supplementary_figures/fig_s18.R
```

The scripts locate the repository root automatically and use the project
`renv` library when it is available. Publication outputs are vector PDFs in
`output/supplementary/`; the output directory contains no figure-specific
subdirectories, PNGs, CSVs, reports, or model objects. Small tidy tables
containing the values shown in each figure are stored separately under
`data/final/figure_data/`.

## Reconstruction of `Supplementary_v1.pdf`

Each figure has a dedicated numbered entry script.
Shared data preparation, bootstrap, modelling, plotting, and export code is in
`scripts/supplementary_figures/figure_helpers.R`; the figure scripts only specify the relevant
dataset, response, species, and labels. Run the complete collection with:

```sh
ALINV_SUPP_BOOT_B=1000 ALINV_SUPP_BOOT_CORES=8 \
  Rscript --vanilla scripts/supplementary_figures/make_all_figures.R
```

All 18 outputs are written as `fig_s1.pdf` through `fig_s18.pdf` directly to
`output/supplementary/`. No figure-specific subdirectories,
PNGs, CSVs, reports, or model objects are written there. During a run, fitted
bootstrap objects may be cached only in the R session's temporary directory;
`ALINV_SUPP_MODEL_CACHE` can optionally point to a cache outside the
supplementary output directory.

| PDF | Entry script | Analysis shown |
|---|---|---|
| `fig_s1.pdf` | `fig_s1.R` | Soil C, N, and C:N contrasts |
| `fig_s2.pdf` | `fig_s2.R` | Measurement schedule |
| `fig_s3.pdf` | `fig_s3.R` | Diameter-increment time series |
| `fig_s4.pdf` | `fig_s4.R` | Height-increment time series |
| `fig_s5.pdf` | `fig_s5.R` | Quantum-yield time series |
| `fig_s6.pdf` | `fig_s6.R` | Chlorophyll time series |
| `fig_s7.pdf` | `fig_s7.R` | Vitality time series |
| `fig_s8.pdf` | `fig_s8.R` | Spring phenology progression |
| `fig_s9.pdf` | `fig_s9.R` | Autumn-senescence time series |
| `fig_s10.pdf` | `fig_s10.R` | Height temporal-LMM effects |
| `fig_s11.pdf` | `fig_s11.R` | Vitality temporal-LMM effects |
| `fig_s12.pdf` | `fig_s12.R` | Remaining-green temporal-LMM effects |
| `fig_s13.pdf` | `fig_s13.R` | Chlorophyll temporal-LMM effects |
| `fig_s14.pdf` | `fig_s14.R` | *Fagus* biomass distributions |
| `fig_s15.pdf` | `fig_s15.R` | *Quercus* biomass distributions |
| `fig_s16.pdf` | `fig_s16.R` | Overall phenology timing effects |
| `fig_s17.pdf` | `fig_s17.R` | Exploratory SWC-path decomposition |
| `fig_s18.pdf` | `fig_s18.R` | Phase-wise relative volume-increment temporal-LMM effects |

### Bootstrap used for the reconstructed figures

For the temporal inferential figures, response values are z-standardized once
using the original species-specific analysis sample. At each measurement date,
the model estimates drought-minus-control, with-minus-without-Robinia, and
mixed-minus-monoculture effects:

```text
y_z ~ date + date:precipitation + date:Robinia + date:culture
      + (1 | container) + (1 | tree)
```

Uncertainty is estimated with 1,000 successful nonparametric bootstrap refits.
Whole containers are sampled with replacement within experimental block; all
trees and dates from a sampled container remain together and receive synthetic
container and tree identifiers. Nonconverged refits are rejected and retried;
converged singular fits are retained and counted during the run. The plotted
intervals are 2.5th and 97.5th percentiles of the successful bootstrap
estimates.

The descriptive time-series figures do not fit the LMM. Their uncertainty bands
use the same design-aligned, block-stratified container resampling principle on
container-level group means. Figure S2.6 retains mean +/- SE to match the
phenology-progression description, and Figures S2 and S14-S15 are purely
descriptive.

## Figures

- `fig_s8.pdf`: observed treatment-group mean transition
  day of year (DOY) with standard errors. Stages 1–4 are shown to describe the
  complete spring progression; stage 1 is not included in the inferential
  model.
- `fig_s16.pdf`: treatment contrasts from the
  primary common-shift phenology model, with 95% block-stratified
  container-bootstrap intervals.
- `fig_s17.pdf`: exploratory small-multiple heatmaps of
  SWC-adjusted treatment-response associations, SWC-associated indirect paths,
  and their uncertainty. Repeated responses use same-day or latest preceding
  measured SWC within seven days; future SWC is never assigned to an earlier
  response. Both components use one shared diverging scale. The figure is an
  associative decomposition and is not evidence of causal mediation. Cell
  intervals and probabilities use a 1,000-replicate block-stratified
  container-cluster bootstrap. Figure 6 reports their path-summed total
  (`c' + ab`).
- `fig_s18.pdf`: date-specific treatment effects on phase-wise relative volume
  increments, using the same response as the Figure 6 time series and the
  shared temporal LMM above. The figure data are saved in
  `data/final/figure_data/fig_s18_effects.csv`.

### Phase-wise relative volume increments (S18)

The response is `volume / phase_volume_baseline - 1`. The initial measurement
(7 March) is the spring baseline; the last measurements before the next phases
(11 June and 1 September) are the summer and autumn baselines. The 1 September
growth round is still assigned to summer. The initial all-zero increment is
excluded before standardization. The artificial zero-reset rows drawn in
Figure 6 are not used for modelling. Genuine observed negative increments are
retained.

One model per species estimates the three treatment contrasts at each of eight
post-baseline measurement dates, with container and tree random intercepts.
The response SD is fixed across dates and bootstrap samples within each species.
S18 uses 1,000 successful container-bootstrap refits per species and eight
workers by default. Block is a resampling stratum, not a model predictor.
Lines and confidence bands stop at phase boundaries because the denominator
changes; the figure has no added zero-reset observations.

Suggested caption: **Figure S18.** Temporal treatment effects on phase-wise
relative volume increments in Fagus and Quercus. Points show standardized
treatment-minus-reference contrasts from species-specific linear mixed-effects
models; bands show 95% container-bootstrap confidence intervals (1,000 refits).
Negative values indicate smaller increments. Vertical lines mark phase
boundaries and orange bars indicate extreme-event periods. Lines are separated
where the phase baseline changes.

These are accumulated increments within each phase, not instantaneous growth
rates. A summer-specific reduction can leave a persistent size deficit even
when later phase-wise contrasts are near zero. However, the autumn denominator
is already treatment-affected, so an autumn contrast near zero does not by
itself demonstrate recovery or catch-up growth. Evidence of a stronger summer
effect requires comparing coefficients across dates, not just their separate
significance labels.

## Statistical model

One Gaussian linear mixed-effects model is fitted separately for *Fagus
sylvatica* and *Quercus ilex* (shown as Fagus and Quercus in the figures):

```text
transition DOY ~ phenological stage + precipitation + Robinia + culture
               + (1 | container) + (1 | tree)
```

Only transitions into stages 2–4 are used for inference. Phenological stage,
precipitation, *Robinia* presence, and culture are fixed effects; container and
tree identity are random intercepts. The additive
treatment coefficients estimate a common shift in the stage-2–4 transition
profile. Treatment contrasts are treatment minus reference (drought minus
control, with minus without *Robinia*, and mixed culture minus monoculture).

Raw DOY contrasts are retained in the bootstrap audit output: negative values
mean earlier and positive values mean later attainment. Displayed phenology
effects are multiplied by -1 so that negative values mean later and positive
values earlier leaf-out in Figures 6, S16, and the leaf-out cells of S17.
Because DOY is continuous, this is an LMM rather than a GLMM.

Production model-fitting scripts are under `scripts/model_fitting/`, and the
compact summaries used by the figures are under `data/final/bootstrap/`.
Replicate-level audit artifacts remain local and are not stored with the
publication PDFs.
