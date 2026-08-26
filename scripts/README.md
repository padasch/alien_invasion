# Analysis scripts

The active code is separated by purpose. Every directory below is flat; files
do not depend on a numerical execution prefix.

- `data_preparation/`: executable scripts that rebuild tracked interim data
  from local raw inputs.
- `functions/`: reusable function definitions sourced by analyses and figures.
- `model_fitting/`: executable bootstrap and SEM model-refitting scripts.
- `main_figures/`: one entry script per main figure plus a collection runner.
- `supplementary_figures/`: one entry script per supplementary figure plus a
  collection runner.
- `config.R`: shared factor levels, treatment definitions, response
  orientations, and volume-proxy settings.
- `run_all.R`: rebuilds both publication figure collections.

The active pipeline uses measured SWC only. The retired GAM-based interpolation
script and its generated daily dataset are stored locally in
`_archive/gam-swc/`.

Compact model summaries consumed by the figures are tracked under
`data/final/bootstrap/`. Large replicate-level draws and model caches are
regenerated locally and are not publication outputs.
