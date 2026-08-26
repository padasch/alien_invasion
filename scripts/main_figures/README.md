# Main figures

Run the complete main-figure collection from the repository root:

```sh
Rscript --vanilla scripts/main_figures/make_all_figures.R
```

The scripts share setup and plotting helpers through `setup.R`. Final PDFs
are written directly to `output/main_figures/`; no model objects or tabular
intermediates are written to that publication-output directory.
