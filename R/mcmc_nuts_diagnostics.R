


source("R/functions.R")

fit_files <- fs::dir_ls("saved_results/results_baseline/2025-10-09/", recurse = TRUE, glob = "*.rds")

fit_diagnostics <- fit_files |> 
  map(function(fit_file) {
    fit <- read_rds(fit_file)
    fit$diagnostic_summary() |> 
      as_tibble() |> 
      mutate(chain = row_number()) |> 
      mutate(round_id = ymd(str_extract(fit_file, "\\d{4}-\\d{2}-\\d{2}")),
             location_pathogen = fs::path_file(fit_file) |> fs::path_ext_remove())
  }, .progress = TRUE) |> 
  bind_rows()




fit_diagnostics |> 
  ggplot() +
  geom_point(aes(y = location_pathogen, x = num_divergent))


fit_diagnostics |> 
  arrange(desc(num_divergent))

fit <- read_rds("saved_results/results_baseline/2025-10-09/RSV_QLD.rds")


fit_posterior <- as.array(fit$draws())
fit_np <- nuts_params(fit)


pars <- c(
  "scaling_annual", "scaling_weekly_asc", "scaling_noise", "ls_noise", "ls_annual", "phi",
  "delay_mean", "delay_sd"
)


fit_posterior[is.na(fit_posterior)] <- -1

mcmc_pairs(fit_posterior, np = fit_np, pars = pars,
           off_diag_args = list(size = 0.5))


mcmc_nuts_energy(fit_np)
