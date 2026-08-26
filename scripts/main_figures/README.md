# Main figures

Run the complete main-figure collection from the repository root:

```sh
Rscript --vanilla scripts/main_figures/make_all_figures.R
```

The scripts share setup and plotting helpers through `setup.R`. Final PDFs
are written directly to `output/main_figures/`. Each script also writes one or
more tidy tables containing the values shown in the figure under
`data/final/figure_data/`; the publication-output directory remains PDF-only.
