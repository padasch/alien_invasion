# Harvest-biomass LMM bootstrap

This exploration retains the six models currently used for Figure 5—three
harvest metrics fitted separately for *Fagus* and *Quercus*—and replaces the
Wald uncertainty calculation with a block-stratified container-cluster
bootstrap.

Run from the repository root:

```sh
Rscript --vanilla "exploration/2026-08-18 Testing bootstrapping/biomass_lmm/scripts/01-biomass-container-bootstrap.R"
```

The script targets 1,000 successful refits per species–metric model, writes
point and bootstrap source data, status diagnostics, an RDS model bundle, a
Wald-versus-bootstrap comparison, and bootstrap PDF/PNG versions of Figure 5
under `output/`.
