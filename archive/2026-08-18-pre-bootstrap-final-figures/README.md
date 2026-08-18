# Pre-bootstrap final-figure workflow

This directory is a read-only snapshot of the final-figure workflow immediately
before the publication figures were migrated to design-aligned bootstrap
uncertainty on 18 August 2026.

The snapshot preserves:

- the complete `final_figures/main` R pipeline and its rendered PDFs;
- the project functions called directly by that pipeline;
- the analysis configuration; and
- the contemporaneous statistical Methods text.

The archived workflow used several uncertainty estimators: standard errors for
descriptive trajectories, asymptotic/Wald intervals for temporal and biomass
mixed models, variance-propagated intervals for most repeated-response SEMs,
and a container bootstrap for the newer phenology SEM. It is retained for
provenance and result comparison, not as the active publication workflow.

The copied scripts still expect the repository's original `data/` and `output/`
trees. Run the active scripts under `final_figures/main/` for the maintained
bootstrap implementation.
