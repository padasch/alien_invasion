# Partial-mediation SEM exploration

This exploration evaluates the proposed figure structure using one coherent SEM:

- main total effect: `c' + a*b`;
- SWC-adjusted component: `c'`;
- SWC-associated indirect component: `a*b`.

The repeated-response analysis uses measured SWC from the same day or the latest preceding measurement within seven days. The phenology analysis uses the pre-leaf-out container mean SWC already defined by the phenology bootstrap analysis.

Run from the repository root:

```sh
Rscript "exploration/2026-08-26 Partial mediation SEM/scripts/01_sem_total_and_decomposition.R"
```

The script reuses the completed 1,000-replicate block-stratified container-bootstrap draws, recalculates intervals and probabilities from those draws, validates `total = direct + indirect` within every replicate, and writes the figure, tables, and interpretation to this directory.
