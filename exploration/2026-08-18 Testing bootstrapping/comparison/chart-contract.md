# Chart contract

## Recreated manuscript figures

- **Question:** Do the published treatment-effect patterns remain visually and
  inferentially similar when uncertainty is estimated by container-cluster
  bootstrap?
- **Takeaway sought:** Direction and magnitude should remain unchanged because
  point-estimate models are retained; interval width and zero-exclusion may
  change.
- **Forms:** Original descriptive time-series layouts; faceted dot-and-interval
  plots for temporal, biomass, and phenology effects; diverging heatmap for
  Figure 6 path-summed SEM totals.
- **Data sufficiency:** All original dates and treatment contrasts are retained;
  every inferential interval targets 1,000 successful replicates.
- **Palette:** Preserve the manuscript treatment and signed-effect palettes.
  Line type, point shape, facets, and zero lines provide non-colour distinction.
- **Outputs:** Vector PDF and publication-resolution PNG under
  `../recreated_figures/`, with source CSVs retained in each model-family
  directory.
- **QA:** Inspect exported PNGs for clipping, sign orientation, consistent
  scales, legible legends, and correspondence with their CSV source data.

## Comparison summary

- **Question:** Which treatment conclusions change under bootstrap inference?
- **Form:** Compact model-family table plus machine-readable effect-level CSV;
  use a comparison plot only if the number of changed conclusions is large
  enough to benefit from one.
- **Primary comparisons:** point-estimate sign; old and bootstrap CI
  zero-exclusion; old and bootstrap two-sided significance; interval width;
  number of successful bootstrap replicates.
- **Caveat:** Temporal date-specific contrasts and whole-season SEM totals have
  different estimands. Cross-model agreement is evaluated by direction and
  qualitative pattern, not equality of coefficients.
