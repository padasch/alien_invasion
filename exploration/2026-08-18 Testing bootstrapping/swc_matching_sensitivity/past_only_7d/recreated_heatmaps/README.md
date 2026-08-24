# Past-only SWC heatmap comparison

These PDFs replace only the repeated-response SEM input with the past-only measured-SWC sensitivity (same-day or latest preceding SWC within seven days). All available models use 1,000 successful block-stratified container-cluster bootstrap replicates.

- `fig6-volume-sem-past-only-swc.pdf` should be compared with `output/main_figures/fig6_volume_sem.pdf`.
- `figure-s3-sem-path-decomposition-past-only-swc.pdf` should be compared with the current version at `output/supplementary/fig_s17.pdf`.

The phenology/leaf-out SEM is unchanged because it does not use the repeated-response fuzzy-matching rule. Both past-only variants retain the baseline heatmap color limits so color intensity is directly comparable. Coefficients are printed only where the relevant 1,000-replicate bootstrap p-value is below 0.05; grey cells are unavailable.

Older supplementary SEM exports were moved to `_archive`; they predate the current repeated-response bootstrap source and should not be used for the uncertainty-consistent comparison.
