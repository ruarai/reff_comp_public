# Holiday- and seasonality-aware renewal forecasts of respiratory viruses

Code accompanying the manuscript:

> Tobin RJ, Eales O, Plank MJ, Morris DJ, McCaw JM, Shearer FM. **Improving respiratory virus forecasts by capturing the effects of school holidays, public holidays and seasonal changes in transmission.** *Preprint* (2026). doi: TODO

The model is a Bayesian renewal model. The log effective reproduction number is decomposed into additive components:

```
log R(t) = ω_S(t) + ω_P(t) + f_seasonal(t) + f_noise(t) + c
```

- `ω_S(t)` and `ω_P(t)` are school and public holiday effects.
- `f_seasonal(t)` is an annually periodic Gaussian process.
- `f_noise(t)` is a Matérn 5/2 Gaussian process that captures remaining variation.
- `c` is a baseline.

Infections are convolved with an estimated infection-to-notification delay to give expected case counts. Expected counts are further adjusted by school holiday, public holiday and day-of-week reporting effects, and observed counts are negative binomial. Both Gaussian processes use Hilbert space (HSGP) approximations. The model is fitted in Stan through `cmdstanr`. The study produced 28-day-ahead forecasts of RSV, influenza and SARS-CoV-2 daily case incidence in five Australian jurisdictions (ACT, NSW, NT, QLD, VIC) across the 2025 season.

The model also produces real-time forecasts for the Australia-Aotearoa Forecasting Hub under the model name `uom-reff_comp`.

## Data availability

The case data come from line lists in the National Notifiable Disease Surveillance System (NNDSS) and, for Victoria, the Public Health Event Surveillance System (PHESS). Their use was approved by the University of Melbourne Health Sciences Human Ethics Sub-Committee (ID 2024-26949-50575-3). These data cannot be shared publicly; requests for access should go to the data custodians.

So that the code can still be run end to end, `data/` contains **anonymised case time series**. These match the structure of the true data but have noise added. Results produced from them will be close to, but not identical to, those in the manuscript.

### Expected layout of `data/`

```
data/
├── retro_rounds/<round_id>/
│   ├── combined-count-<round_id>.csv            # defines the jobs (location × pathogen × target) for a round
│   └── combined-count-smoothed-<round_id>.csv   # counts used for fitting
├── combined-count-smoothed_final_anonymised.csv # final (anonymised) case counts: used for scoring and all figures
└── paper/                                       # posterior summaries written by R/paper/collect_effect_draws.R
```

The holiday dates used in the study are hard-coded in [R/holidays.R](R/holidays.R), with a source URL and verification date for each.

## Requirements

- R ≥ 4.5.0
- [CmdStan](https://mc-stan.org/cmdstanr/) installed through `cmdstanr`
- The study was run on the University of Melbourne Spartan HPC cluster (Slurm). Fitting a single job takes roughly an hour or more with 4 cores.

```r
install.packages(c(
  "tidyverse", "tidybayes", "bayesplot", "patchwork", "arrow", "fs", "glue",
  "scales", "posterior", "matrixStats", "scoringRules", "distributional",
  "fable", "fabletools", "tsibble", "ggdist", "ggtext", "ggh4x", "ggokabeito",
  "RColorBrewer", "scico", "pammtools", "magick", "tinytable", "processx"
))

install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
cmdstanr::install_cmdstan()

# Used for plotting credible-interval ribbons
remotes::install_github("idd-unit/fastribbon")
```

Run all scripts from the repository root (for example, by opening `periodicgpcast.Rproj`).

## Reproducing the manuscript

### 1. Fit the retrospective forecasts

[R/forecast_renewal_retro.R](R/forecast_renewal_retro.R) takes a single integer job index. The job list is every round in `data/retro_rounds/` × location × pathogen × target, crossed with the model versions listed in the `expand_grid(model_version = ...)` call near the top of the script. Uncomment the versions you want to run:

| `model_version`         | Used in                       | Description                                                  |
|-------------------------|-------------------------------|--------------------------------------------------------------|
| `baseline`              | Figs 1–5, Supp. Table 1       | Full model                                                   |
| `no_holidays_no_annual` | Fig 4                         | No holiday effects and no seasonal GP (the null model)       |
| `no_holidays`           | Fig 4, Supp. Fig 8            | No school or public holiday effects (on R(t) or reporting)   |
| `no_annual`             | Fig 4, Supp. Fig 8            | No seasonal GP                                               |
| `no_school_holidays`    | Supp. Fig 8                   | No school holiday effects                                    |
| `no_public_holidays`    | Supp. Fig 8                   | No public holiday effects                                    |
| `no_asc_holidays`       | Supp. Fig 8                   | No holiday effects on reporting                              |
| `no_Rt_holidays`        | Supp. Fig 8                   | No holiday effects on R(t)                                   |
| `no_delay`              | Supp. Fig 8                   | No infection-to-notification delay                           |

Components are removed by setting the corresponding flags in the Stan data (see `assign_retro_varietals()` in [R/functions_fitting.R](R/functions_fitting.R)), which fixes those components at zero. All versions use [src/gp_renewal_varietal.stan](src/gp_renewal_varietal.stan).

To run one job locally:

```bash
Rscript R/forecast_renewal_retro.R 1
```

Each job writes the following to `saved_results/<model_version>/<round_id>/`:

- the `cmdstanr` fit object (`.rds`)
- the training data
- forecast trajectories
- CRPS performance metrics
- diagnostic plots

Rounds are fitted with 2,000 warm-up and 2,000 sampling iterations, except the `baseline` fit for the final round (`2025-11-27`), which uses 5,000 warm-up and 20,000 sampling iterations.

**On a cluster:** [run_retro.slurm](run_retro.slurm) submits the jobs as a Slurm array. The working directory, `R_LIBS_USER`, the partitions and the `my-job-stats` call are specific to Spartan, so edit them before using the script elsewhere. Set `--array` to cover the number of jobs printed by the script ("Generated N jobs"). Indices above N exit immediately.

### 2. Refit jobs that did not converge

[R/forecast_renewal_retro_refit.R](R/forecast_renewal_retro_refit.R) takes the same job indices and model-version list. It loads each existing fit and checks R̂. If R̂ ≥ 1.05, it refits with longer chains. Submit it with [run_retro_refit.slurm](run_retro_refit.slurm), repeating until every job converges. [R/paper/check_convergence.R](R/paper/check_convergence.R) summarises the convergence of a run.

### 3. Baseline forecasts and summaries

```bash
Rscript R/paper/produce_snaive_forecasts.R   # random-walk (seasonal naive) baseline → data/snaive_perf_extended.parquet
Rscript R/paper/collect_effect_draws.R       # posterior draws from the final round → data/paper/*.parquet
Rscript R/paper/collect_trajectories.R       # forecast quantiles → saved_results/baseline_quants.parquet
```

### 4. Figures and tables

Each script in `R/paper/` writes a single figure to `paper_figures/` or `paper_figures/supp/`:

| Output                                   | Script                                                                 |
|------------------------------------------|------------------------------------------------------------------------|
| Figure 1                                 | [f_1_example_forecast.R](R/paper/f_1_example_forecast.R)               |
| Figure 2                                 | [f_2_performance.R](R/paper/f_2_performance.R)                         |
| Figure 3                                 | [f_3_examples.R](R/paper/f_3_examples.R)                               |
| Figure 4                                 | [f_4_model_comparison.R](R/paper/f_4_model_comparison.R)               |
| Figure 5                                 | [f_5_estimated_effects.R](R/paper/f_5_estimated_effects.R)             |
| Supp. Table 1                            | [t_supp_estimated_effects.R](R/paper/t_supp_estimated_effects.R)       |
| Supp. Fig 1 (holiday dates)              | [f_supp_public_holidays.R](R/paper/f_supp_public_holidays.R)           |
| Supp. Figs 2–4 (case incidence)          | [f_supp_case_incidence.R](R/paper/f_supp_case_incidence.R)             |
| Supp. Fig 5 (CRPS by horizon)            | [f_supp_perf_heatmap.R](R/paper/f_supp_perf_heatmap.R)                 |
| Supp. Fig 6 (forecasts with and without components) | [f_supp_forecast_components.R](R/paper/f_supp_forecast_components.R) |
| Supp. Fig 7 (skill heatmaps)             | [f_supp_skill_heatmap.R](R/paper/f_supp_skill_heatmap.R)               |
| Supp. Fig 8 (component ablations)        | [f_supp_model_comparison.R](R/paper/f_supp_model_comparison.R)         |
| Supp. Figs 9–10 (holiday effects)        | [f_supp_estimated_effects.R](R/paper/f_supp_estimated_effects.R)       |
| Supp. Fig 11 (delay distributions)       | [f_supp_delay.R](R/paper/f_supp_delay.R)                               |
| Supp. Fig 12 (HSGP validation)           | [f_supp_hsgp.R](R/paper/f_supp_hsgp.R)                                 |
| Supp. Figs 13–17 (all forecasts)         | [f_supp_all_forecasts.R](R/paper/f_supp_all_forecasts.R)               |

## Repository map

| Path | Contents |
|------|----------|
| [src/](src/) | Stan models. `gp_renewal_varietal.stan` is the manuscript model. `gp_functions.stan` contains the HSGP functions adapted from [EpiNow2](https://github.com/epiforecasts/EpiNow2). |
| [R/functions*.R](R/) | Data processing, Stan data construction, fitting, forecast output, scoring and plotting helpers. |
| [R/holidays.R](R/holidays.R) | School and public holiday dates for Australian states and territories and New Zealand. |
| [R/forecast_renewal_retro*.R](R/) | Retrospective fitting entry points (see above). |
| [R/paper/](R/paper/) | Manuscript figures, tables and summaries. |
| `paper_figures/` | Rendered manuscript figures. |
| `R/paper/pnas/`, `pnas_figures/` | Alternative figure layouts for a different journal format. Not needed for the manuscript. |
| `workflow.R`, `reporting.R`, `R/forecast_renewal.R`, `run_forecasts.slurm` | Operational real-time pipeline that submits weekly forecasts to the Forecasting Hub. Not needed for the manuscript. |
| `R/mediaflux_integration.R` | Upload of results to University of Melbourne Mediaflux storage. |
| `R/scratch/`, `R/plots/`, `R/presentation/`, `R/forecast_renewal_scratch.R`, `src/gp_renewal_sus.stan` | Exploratory model variants (hierarchical, immunity-based, locally periodic) and presentation plots. Not used in the manuscript. |

## Citation

If you use this code, please cite the manuscript:

```
Tobin RJ, Eales O, Plank MJ, Morris DJ, McCaw JM, Shearer FM (2026).
Improving respiratory virus forecasts by capturing the effects of school holidays,
public holidays and seasonal changes in transmission. doi: TODO
```

## Licence

MIT. See [LICENSE](LICENSE).
