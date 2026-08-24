# Supplementary information

This directory contains reproducible scripts for the supplementary figures,
including phenology, temporal responses, biomass, and SEM decompositions. Run
the scripts from the repository root:

```sh
Rscript --vanilla supplementary_information/scripts/01-phenology-progression.R
Rscript --vanilla supplementary_information/scripts/02-phenology-overall-timing-lmm.R
Rscript --vanilla supplementary_information/scripts/03-sem-path-decomposition-heatmaps.R
```

The scripts locate the repository root automatically and use the project
`renv` library when it is available. Publication outputs are vector PDFs in
`supplementary_information/output/v1/`; the output directory contains no
figure-specific subdirectories, PNGs, CSVs, reports, or model objects.

## Reconstruction of `Supplementary_v1.pdf`

The 15 figures in `Supplementary_v1.pdf` each have a dedicated entry script.
Shared data preparation, bootstrap, modelling, plotting, and export code is in
`scripts/_v1-figure-helpers.R`; the figure scripts only specify the relevant
dataset, response, species, and labels. Run the complete collection with:

```sh
ALINV_SUPP_BOOT_B=1000 ALINV_SUPP_BOOT_CORES=4 \
  Rscript --vanilla supplementary_information/scripts/make-v1-figures.R
```

All 15 outputs are written as PDFs directly to
`supplementary_information/output/v1/`. No figure-specific subdirectories,
PNGs, CSVs, reports, or model objects are written there. During a run, fitted
bootstrap objects may be cached only in the R session's temporary directory;
`ALINV_SUPP_MODEL_CACHE` can optionally point to a cache outside the
supplementary output directory.

| Figure | Entry script | Analysis shown |
|---|---|---|
| S1.1 | `fig-s1-1-soil-elemental.R` | Soil C, N, and C:N; block-stratified container bootstrap contrasts |
| S1.2 | `fig-s1-2-measurement-schedule.R` | Measurement schedule (descriptive) |
| S2.1 | `fig-s2-1-diameter-timeseries.R` | Diameter-increment group means and container-bootstrap intervals |
| S2.2 | `fig-s2-2-height-timeseries.R` | Height-increment group means and container-bootstrap intervals |
| S2.3 | `fig-s2-3-quantum-yield-timeseries.R` | Fv/Fm group means and container-bootstrap intervals |
| S2.4 | `fig-s2-4-chlorophyll-timeseries.R` | Chlorophyll group means and container-bootstrap intervals |
| S2.5 | `fig-s2-5-vitality-timeseries.R` | Vitality group means and container-bootstrap intervals |
| S2.6 | `fig-s2-6-spring-phenology.R` | Spring phenology progression (mean +/- SE, descriptive) |
| S2.7 | `fig-s2-7-autumn-senescence-timeseries.R` | Senesced-crown group means and container-bootstrap intervals |
| S3.1 | `fig-s3-1-height-effects-bootstrap.R` | Height temporal LMM effects with bootstrap intervals |
| S3.2 | `fig-s3-2-vitality-effects-bootstrap.R` | Vitality temporal LMM effects with bootstrap intervals |
| S3.3 | `fig-s3-3-senescence-effects-bootstrap.R` | Remaining-green temporal LMM effects with bootstrap intervals |
| S4.1 | `fig-s4-1-chlorophyll-effects-bootstrap.R` | Chlorophyll temporal LMM effects with bootstrap intervals |
| S5.1 | `fig-s5-1-biomass-fagus.R` | Fagus biomass distributions (descriptive) |
| S5.2 | `fig-s5-2-biomass-quercus.R` | Quercus biomass distributions (descriptive) |

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
phenology-progression description, and Figures S1.2 and S5.1-S5.2 are purely
descriptive.

## Figures

- `figure-s1-phenology-progression`: observed treatment-group mean transition
  day of year (DOY) with standard errors. Stages 1–4 are shown to describe the
  complete spring progression; stage 1 is not included in the inferential
  model.
- `figure-s2-phenology-overall-timing-effects`: treatment contrasts from the
  primary common-shift phenology model, with 95% block-stratified
  container-bootstrap intervals.
- `figure-s3-sem-path-decomposition`: exploratory small-multiple heatmaps of
  SWC-adjusted treatment-response associations, SWC-associated indirect paths,
  and treatment-to-SWC paths. Repeated responses use same-day or latest
  preceding measured SWC within seven days; future SWC is never assigned to an
  earlier response. All components use one shared diverging scale. The figure
  is an associative decomposition and is not evidence of causal mediation.
  Cell intervals and probabilities use a 1,000-replicate block-stratified
  container-cluster bootstrap. Figure 6 instead reports directly fitted
  reduced-form treatment effects without SWC.

## Statistical model

One Gaussian linear mixed-effects model is fitted separately for *Fagus
sylvatica* and *Quercus ilex* (shown as Fagus and Quercus in the figures):

```text
transition DOY ~ phenological stage + precipitation + Robinia + culture + block
               + (1 | container) + (1 | tree)
```

Only transitions into stages 2–4 are used for inference. Phenological stage,
precipitation, *Robinia* presence, culture, and experimental block are fixed
effects; container and tree identity are random intercepts. The additive
treatment coefficients estimate a common shift in the stage-2–4 transition
profile. Treatment contrasts are treatment minus reference (drought minus
control, with minus without *Robinia*, and mixed culture minus monoculture).

Raw DOY contrasts are retained in the bootstrap audit output: negative values
mean earlier and positive values mean later attainment. The plotted contrast is multiplied
by -1 solely to align its interpretation with the other stress responses:
negative plotted values mean delayed phenology and positive values mean earlier
phenology. Because DOY is continuous, this is an LMM rather than a GLMM.

Bootstrap model artifacts and audit tables are retained in
`exploration/2026-08-18 Testing bootstrapping/`; they are intentionally kept
separate from the publication-only supplementary PDF directory.
