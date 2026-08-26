# Model-fitting scripts

These scripts implement the container-cluster bootstrap analyses used by the
publication figures:

- `temporal_models.R`: temporal treatment effects;
- `biomass_models.R`: harvest biomass effects;
- `phenology_model.R`: phenology timing and its SWC decomposition;
- `response_sem_sources.R`: current-data SEM source models;
- `response_sem.R`: repeated-response SEM paths;
- `past_swc_sem.R`: SEM paths using same-day or preceding SWC
  within seven days.

The plotting pipeline normally reads validated 1,000-replicate summaries from
`data/final/bootstrap/`. If a required summary is missing or incomplete, the
shared bootstrap reader reruns the corresponding production script. Large
replicate-level draws, fitted objects, diagnostic graphics, and refreshed SEM
source models are local analysis artifacts rather than publication outputs.
