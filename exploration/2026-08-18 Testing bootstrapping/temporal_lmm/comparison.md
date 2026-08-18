# Temporal LMM bootstrap comparison

The comparison uses 1000 successful non-parametric container-cluster bootstrap replicates per species-response model. Containers were resampled within block, with all nested trees and dates retained. The response z-scaling and the current temporal LMM formula were held fixed.

## Validation summary

- Successful models: 4/4.
- Successful replicates: 4000; failed attempts: 16; singular successful fits: 345.
- Maximum absolute difference between refitted and saved current point estimates: 0.000000000040631942 SD.
- Point-estimate sign disagreements: 0.
- CI-based significance classifications changed for 10 of 120 date-specific contrasts.
- Among those changes, 5 bootstrap intervals newly excluded zero and 5 no longer excluded zero.

## By response and species

| Response | Species | Contrasts | Median width ratio | Width-ratio range | Significance flips |
|---|---:|---:|---:|---:|---:|
| diameter_inc_t0 | fagus | 27 | 1.20 | 0.107–1.45 | 1 |
| diameter_inc_t0 | quercus | 27 | 0.875 | 0.0247–1.58 | 3 |
| qy | fagus | 33 | 1.17 | 0.424–1.49 | 1 |
| qy | quercus | 33 | 1.00 | 0.611–1.25 | 5 |

## Changed CI classifications

| Response | Species | Date | Effect | Current excludes zero | Bootstrap excludes zero | Current CI | Bootstrap CI |
|---|---|---|---|---:|---:|---:|---:|
| diameter_inc_t0 | fagus | 2025-03-07 | robinia | FALSE | TRUE | [-0.209, 0.113] | [-0.0785, -0.0101] |
| diameter_inc_t0 | quercus | 2025-05-08 | robinia | FALSE | TRUE | [-0.221, 0.0283] | [-0.158, -0.0451] |
| diameter_inc_t0 | quercus | 2025-09-29 | culture | TRUE | FALSE | [-0.326, -0.0635] | [-0.398, 0.00855] |
| diameter_inc_t0 | quercus | 2025-12-02 | culture | TRUE | FALSE | [-0.305, -0.0431] | [-0.360, 0.00619] |
| qy | fagus | 2025-08-14 | culture | TRUE | FALSE | [-0.545, -0.0308] | [-0.613, 0.0264] |
| qy | quercus | 2025-07-31 | precipitation | FALSE | TRUE | [-0.543, 0.0553] | [-0.471, -0.00287] |
| qy | quercus | 2025-08-29 | robinia | TRUE | FALSE | [-0.600, -0.00164] | [-0.617, 0.00266] |
| qy | quercus | 2025-09-10 | culture | FALSE | TRUE | [-0.0841, 0.547] | [0.0262, 0.446] |
| qy | quercus | 2025-10-20 | culture | FALSE | TRUE | [-0.0121, 0.619] | [0.0432, 0.569] |
| qy | quercus | 2025-11-12 | precipitation | TRUE | FALSE | [0.0169, 0.615] | [-0.00455, 0.709] |

## Interpretation

The two approaches estimate the same treatment contrasts: point estimates are unchanged apart from numerical fitting tolerance. Differences concern uncertainty only. The dominant temporal trajectories are therefore unchanged, and rows above are borderline inferential changes rather than reversals of the estimated biological direction. Date-specific intervals remain pointwise; no multiplicity correction across dates and treatments is applied.

345 of 4000 accepted bootstrap fits were singular, primarily the Quercus diameter model (260/1,000). These fits converged and returned complete contrasts, but at least one random-intercept variance was estimated at the boundary. This signals limited support for both random effects in some resamples and is a reason to interpret isolated date-specific CI flips cautiously.

## Files

- `output/data/temporal-lmm-bootstrap-effects.csv`: plotted bootstrap summaries.
- `output/data/temporal-lmm-bootstrap-draws.csv`: all successful replicate estimates.
- `output/data/temporal-lmm-bootstrap-vs-current.csv`: row-level comparison.
- `output/status/temporal-lmm-bootstrap-status.csv`: convergence and sample-size audit.
- `output/models/`: complete per-model and combined RDS objects.
- `output/figures/`: recreated Figures 2–4 in PDF and PNG.
