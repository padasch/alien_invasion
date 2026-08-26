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

Bootstrap summaries, replicate-level draws, and model caches are generated
locally under `data/final/bootstrap/` and are not version controlled.

Every main and supplementary figure also writes one or more small, tidy CSV
tables under `data/final/figure_data/`. Only the scientific values shown in the
figure are retained; plot styling and other ggplot internals are excluded.
