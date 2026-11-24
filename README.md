# Holiday- and seasonality-aware renewal forecasts of respiratory viruses

Code accompanying the manuscript:

> Tobin RJ, Eales O, Plank MJ, Morris DJ, McCaw JM, Shearer FM. **Improving respiratory virus forecasts by capturing the effects of school holidays, public holidays and seasonal changes in transmission.** *medRxiv pre-print* (2026). [doi: 10.64898/2026.11.15.26367964](https://doi.org/10.64898/2026.11.15.26367964)

The model is a Bayesian renewal model, with the log effective reproduction number decomposed into additive components:

```
log R(t) = ω_S(t) + ω_P(t) + f_seasonal(t) + f_noise(t) + c
```

- `ω_S(t)` and `ω_P(t)` are school and public holiday effects.
- `f_seasonal(t)` is an annually periodic Gaussian process.
- `f_noise(t)` is a Matérn 5/2 Gaussian process that captures remaining variation.
- `c` is a baseline.

Infections are convolved with an estimated infection-to-notification delay to give expected case counts. Expected counts are further adjusted by school holiday, public holiday and day-of-week reporting effects. The model is fitted in Stan through `cmdstanr`. The study produced 28-day-ahead forecasts of RSV, influenza and SARS-CoV-2 daily case incidence in five Australian jurisdictions (ACT, NSW, NT, QLD, VIC) across the 2025 season.

## Data availability

The case data come from line lists in the National Notifiable Disease Surveillance System (NNDSS) and, for Victoria, the Public Health Event Surveillance System (PHESS). Their use was approved by the University of Melbourne Health Sciences Human Ethics Sub-Committee (ID 2024-26949-50575-3). **These data cannot be shared publicly.** Requests for access should go to the data custodians: the Interim Australian Centre for Disease Control for NNDSS, and the Victorian Department of Health for PHESS.

This repository instead includes:

- **Anonymised case counts.** These are the final case time series used in the study, with noise added. Wherever a figure shows observed case counts, it shows these anonymised counts, so those points will differ slightly from the manuscript.
- **Summarised outputs of the retrospective forecasts.** These are forecast scores, forecast quantiles and posterior draws, produced from model fits to the true data. They contain no case counts. The scores were computed against the true counts, so they are released only as weekly means (see [summarise_performance.R](R/paper/summarise_performance.R)). Per-day scores are withheld because they can be inverted to recover the counts.

Together these regenerate every figure and table in the manuscript. Re-fitting the model requires the restricted data (see [Running the full pipeline](#running-the-full-pipeline)).

### Included files

| File | Contents | Written by |
|------|----------|------------|
| `data/combined-count-smoothed_final_anonymised.csv` | Anonymised daily case counts, January 2023 to December 2025 | (provided) |
| `saved_results/perf_summary.parquet` | Mean CRPS (log scale) per model version, pathogen, jurisdiction, round and horizon week, for all nine model versions and the random-walk baseline | [summarise_performance.R](R/paper/summarise_performance.R) |
| `saved_results/baseline_quants.parquet` | Forecast quantiles (50%, 80% and 95% intervals) for every forecast | [collect_trajectories.R](R/paper/collect_trajectories.R) |
| `data/paper/holiday_draws.parquet` | Posterior draws of the school holiday effects (final round) | [collect_effect_draws.R](R/paper/collect_effect_draws.R) |
| `data/paper/public_holiday_draws.parquet` | Posterior draws of the public holiday effects (final round) | [collect_effect_draws.R](R/paper/collect_effect_draws.R) |
| `data/paper/gp_draws.parquet` | Posterior draws of `log R(t)`, `f_seasonal` and `f_noise` across 2025 (final round) | [collect_effect_draws.R](R/paper/collect_effect_draws.R) |
| `data/paper/delay_draws.parquet` | Posterior draws of the infection-to-notification delay (final round) | [collect_effect_draws.R](R/paper/collect_effect_draws.R) |
| `data/paper/plot_data/*.parquet` | Posterior draws for the individual forecasts shown in Figs 1 and 3 and Supp. Fig 6 | `load_fit_trajectories()` in [common.R](R/paper/common.R) |

The holiday dates used in the study are hard-coded in [R/holidays.R](R/holidays.R), with a source URL and verification date for each.

## Requirements

- R ≥ 4.5.0
- The R packages below. `magick` needs the ImageMagick system libraries (`libmagick++-dev` on Debian/Ubuntu); it is only used to composite overlays onto Figs 1–3 after they are saved.
- To fit the model: [CmdStan](https://mc-stan.org/cmdstanr/), installed through `cmdstanr`. The study was run on the University of Melbourne Spartan HPC cluster (Slurm), and a single fit takes roughly an hour or more with 4 cores.

```r
install.packages(c(
  "tidyverse", "tidybayes", "bayesplot", "patchwork", "arrow", "fs", "glue",
  "scales", "posterior", "matrixStats", "scoringRules", "distributional",
  "fable", "fabletools", "tsibble", "ggdist", "ggtext", "ggh4x", "ggokabeito",
  "RColorBrewer", "scico", "pammtools", "magick", "tinytable", "processx"
))

# Used for plotting credible-interval ribbons
remotes::install_github("idd-unit/fastribbon")

# Only needed to fit the model
install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
cmdstanr::install_cmdstan()
```

Run all scripts from the repository root (for example, by opening `periodicgpcast.Rproj`).

## Reproducing the figures and tables

Each script in [R/paper/](R/paper/) writes one figure or table to `paper_figures/` or `paper_figures/supp/`, using only the included files:

```bash
Rscript R/paper/f_2_performance.R
```

| Output | Script | Inputs |
|--------|--------|--------|
| Figure 1 | [f_1_example_forecast.R](R/paper/f_1_example_forecast.R) | plot data, case counts |
| Figure 2 | [f_2_performance.R](R/paper/f_2_performance.R) | score summary (`baseline`, random-walk baseline), case counts |
| Figure 3 | [f_3_examples.R](R/paper/f_3_examples.R) | plot data, `gp_draws`, case counts |
| Figure 4 | [f_4_model_comparison.R](R/paper/f_4_model_comparison.R) | score summary (`baseline`, `no_holidays`, `no_annual`, `no_holidays_no_annual`) |
| Figure 5 | [f_5_estimated_effects.R](R/paper/f_5_estimated_effects.R) | `holiday_draws`, `gp_draws` |
| Supp. Table 1 | [t_supp_estimated_effects.R](R/paper/t_supp_estimated_effects.R) | `holiday_draws` |
| Supp. Fig 1 (holiday dates) | [f_supp_public_holidays.R](R/paper/f_supp_public_holidays.R) | case counts |
| Supp. Figs 2–4 (case incidence) | [f_supp_case_incidence.R](R/paper/f_supp_case_incidence.R) | case counts |
| Supp. Fig 5 (CRPS by horizon) | [f_supp_perf_heatmap.R](R/paper/f_supp_perf_heatmap.R) | score summary (`baseline`), case counts |
| Supp. Fig 6 (forecasts with and without components) | [f_supp_forecast_components.R](R/paper/f_supp_forecast_components.R) | plot data, case counts |
| Supp. Fig 7 (skill heatmaps) | [f_supp_skill_heatmap.R](R/paper/f_supp_skill_heatmap.R) | score summary (`baseline`, `no_holidays_no_annual`, random-walk baseline), case counts |
| Supp. Fig 8 (component ablations) | [f_supp_model_comparison.R](R/paper/f_supp_model_comparison.R) | score summary (all nine model versions) |
| Supp. Figs 9–10 (holiday effects) | [f_supp_estimated_effects.R](R/paper/f_supp_estimated_effects.R) | `holiday_draws`, `public_holiday_draws` |
| Supp. Fig 11 (delay distributions) | [f_supp_delay.R](R/paper/f_supp_delay.R) | `delay_draws` |
| Supp. Fig 12 (HSGP validation) | [f_supp_hsgp.R](R/paper/f_supp_hsgp.R) | none |
| Supp. Figs 13–17 (all forecasts) | [f_supp_all_forecasts.R](R/paper/f_supp_all_forecasts.R) | `baseline_quants`, case counts |

## Running the full pipeline

These steps regenerate the included outputs from scratch. They need the restricted case data, organised as one folder per forecast round:

```
data/retro_rounds/<round_id>/
├── combined-count-<round_id>.csv            # defines the jobs (location × pathogen × target) for the round
└── combined-count-smoothed-<round_id>.csv   # counts used for fitting
```

`.gitignore` excludes these folders and the raw fits in `saved_results/<model_version>/`, so they can't be committed by accident.

### 1. Fit the retrospective forecasts

[R/forecast_renewal_retro.R](R/forecast_renewal_retro.R) takes a single integer job index. The job list is every round in `data/retro_rounds/` × location × pathogen × target, crossed with the model versions listed in the `expand_grid(model_version = ...)` call near the top of the script. Uncomment the versions you want to run:

| `model_version`         | Description                                                  |
|-------------------------|--------------------------------------------------------------|
| `baseline`              | Full model                                                   |
| `no_holidays_no_annual` | No holiday effects and no seasonal GP (the null model)       |
| `no_holidays`           | No school or public holiday effects (on R(t) or reporting)   |
| `no_annual`             | No seasonal GP                                               |
| `no_school_holidays`    | No school holiday effects                                    |
| `no_public_holidays`    | No public holiday effects                                    |
| `no_asc_holidays`       | No holiday effects on reporting                              |
| `no_Rt_holidays`        | No holiday effects on R(t)                                   |
| `no_delay`              | No infection-to-notification delay                           |

Components are removed by setting the corresponding flags in the Stan data (see `assign_retro_varietals()` in [R/functions_fitting.R](R/functions_fitting.R)), which fixes those components at zero. All versions use [src/gp_renewal_varietal.stan](src/gp_renewal_varietal.stan).

To run one job locally:

```bash
Rscript R/forecast_renewal_retro.R 1
```

Each job writes the following to `saved_results/<model_version>/<round_id>/`: the `cmdstanr` fit object (`.rds`), the training data, forecast trajectories, per-day scores (`*_performance.parquet`) and diagnostic plots. Scores are computed against the counts file named in `calculate_performance_metrics()` in [R/functions.R](R/functions.R). It points at the anonymised counts in this repository, but the manuscript's scores were computed against the final restricted counts, so point it at those to reproduce them exactly.

Rounds are fitted with 2,000 warm-up and 2,000 sampling iterations, except the `baseline` fit for the final round (`2025-11-27`), which uses 5,000 warm-up and 20,000 sampling iterations, thinned by 10.

**On a cluster:** [run_retro.slurm](run_retro.slurm) submits the jobs as a Slurm array. The working directory, `R_LIBS_USER`, the partitions and the `my-job-stats` call are specific to Spartan, so edit them before using the script elsewhere. Set `--array` to cover the number of jobs printed by the script ("Generated N jobs"). Indices above N exit immediately.

### 2. Refit jobs that did not converge

[R/forecast_renewal_retro_refit.R](R/forecast_renewal_retro_refit.R) takes the same job indices and model-version list. It loads each existing fit and checks R̂. If R̂ ≥ 1.05, it refits with 10,000 warm-up and 40,000 sampling iterations, thinned by 20. Submit it with [run_retro_refit.slurm](run_retro_refit.slurm), repeating until every job converges. [R/paper/check_convergence.R](R/paper/check_convergence.R) summarises the convergence of a run.

### 3. Baseline model and posterior summaries

```bash
Rscript R/paper/produce_snaive_forecasts.R   # random-walk baseline scores → data/snaive_perf_extended.parquet (not released)
Rscript R/paper/collect_effect_draws.R       # final-round posterior draws → data/paper/*.parquet
Rscript R/paper/collect_trajectories.R       # forecast quantiles → saved_results/baseline_quants.parquet
```

### 4. Summarise the scores

```bash
Rscript R/paper/summarise_performance.R      # → saved_results/perf_summary.parquet
```

This reads the per-day scores for every model version from `saved_results/<model_version>/` and the random-walk baseline scores from `data/snaive_perf_extended.parquet`. It writes their weekly means, which are all the figures use. It stops with an error if any model is missing a forecast date that the others have.


### 5. Plot data for individual forecasts

Figs 1 and 3 and Supp. Fig 6 read their posterior draws through `load_fit_trajectories()` in [R/paper/common.R](R/paper/common.R). If the matching file in `data/paper/plot_data/` exists, it is used. Otherwise, the draws are extracted from the raw fit in `saved_results/` and saved there. To regenerate these files, delete them and re-run the three figure scripts.

## Repository map

| Path | Contents |
|------|----------|
| [src/](src/) | Stan models. `gp_renewal_varietal.stan` is the manuscript model. `gp_functions.stan` contains the HSGP functions adapted from [EpiNow2](https://github.com/epiforecasts/EpiNow2). |
| [R/functions*.R](R/) | Data processing, Stan data construction, fitting, forecast output, scoring and plotting helpers. |
| [R/holidays.R](R/holidays.R) | School and public holiday dates for Australian states and territories and New Zealand. |
| [R/forecast_renewal_retro*.R](R/) | Retrospective fitting entry points (see above). |
| [R/paper/](R/paper/) | Manuscript figures, tables and summaries. |
| `paper_figures/` | Rendered manuscript figures. |
| `workflow.R`, `reporting.R`, `R/forecast_renewal.R`, `run_forecasts.slurm` | Operational real-time pipeline that submits weekly forecasts to the Forecasting Hub. Not needed for the manuscript. |
| `R/mediaflux_integration.R` | Upload of results to University of Melbourne Mediaflux storage. |
| `R/scratch/`, `R/plots/`, `R/presentation/`, `R/forecast_renewal_scratch.R`, `src/gp_renewal_sus.stan` | Exploratory model variants (hierarchical, immunity-based, locally periodic) and presentation plots. Not used in the manuscript, and many read data files that are not included. |

## Citation

If you use this code, please cite the manuscript:

```
Tobin RJ, Eales O, Plank MJ, Morris DJ, McCaw JM, Shearer FM (2026).
Improving respiratory virus forecasts by capturing the effects of school holidays,
public holidays and seasonal changes in transmission. [doi: 10.64898/2026.11.15.26367964](https://doi.org/10.64898/2026.11.15.26367964)
```
