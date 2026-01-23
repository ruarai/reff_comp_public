
fit_posterior <- as.array(fit$draws())
fit_np <- nuts_params(fit)

mcmc_parcoord(fit_posterior, np = fit_np, pars = c("delay_mean", "delay_sd", "scaling_annual"))


pars <- c(
  "scaling_annual", "scaling_weekly_asc", "scaling_noise", "ls_noise", "ls_annual", "phi",
  "delay_mean", "delay_sd"
)


mcmc_pairs(fit_posterior, np = fit_np, pars = pars,
           off_diag_args = list(size = 0.5))


mcmc_nuts_energy(fit_np)