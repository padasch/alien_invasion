# Production bootstrap models

These scripts implement the container-cluster bootstrap analyses used by the
publication figures:

- `01_temporal_lmm_bootstrap.R`: temporal treatment effects;
- `02_biomass_lmm_bootstrap.R`: harvest biomass effects;
- `03_phenology_bootstrap.R`: phenology timing and its SWC decomposition;
- `04_refresh_response_model_sources.R`: current-data SEM source models;
- `05_repeated_response_sem_bootstrap.R`: repeated-response SEM paths;
- `06_past_only_7d_sem_bootstrap.R`: SEM paths using same-day or preceding SWC
  within seven days.

The plotting pipeline normally reads validated 1,000-replicate summaries from
`data/final/bootstrap/`. If a required summary is missing or incomplete, the
shared bootstrap reader reruns the corresponding production script. Large
replicate-level draws, fitted objects, diagnostic graphics, and refreshed SEM
source models are local analysis artifacts rather than publication outputs.
