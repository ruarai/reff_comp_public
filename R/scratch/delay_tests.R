
n_draws <- 50

tibble(
  delay_mean = rnorm(n_draws, 9, 2.0),
  delay_sd = rlnorm(n_draws, 0, 0.5),
  draw = 1:n_draws
) |> 
  filter(delay_mean > 3,
         delay_sd > 1) |> 
  mutate(
    sigma = sqrt(log(1 + (delay_sd ^ 2) / (delay_mean ^ 2))),
    mu = log(delay_mean) - 0.5 * sigma ^ 2
  ) |> 
  
  expand_grid(
    t = 0:21
  ) |> 
  
  mutate(p = plnorm(t, mu, sigma)) |> 
  
  ggplot() +
  
  geom_line(aes(x = t, y = p, group = draw),
            alpha = 0.3) +
  
  theme_bw()
  


x_location <- "QLD"
x_pathogen <- "RSV"
x_target <- "case_incidence"
x_version <- "baseline"
round_id <- "2025-11-27"



draws <- fit$draws()

pars <- c(
  "scaling_annual", "scaling_noise", "ls_noise", "ls_annual", "phi",
  "delay_mean", "delay_sd", "log_R_0"
) |> 
  intersect(fit$metadata()$model_params)

p_hist <- mcmc_hist(posterior::subset_draws(draws, pars), pars = pars) +
  theme_bw()

p_trace <- mcmc_trace(posterior::subset_draws(draws, pars), pars = pars) +
  theme_bw()

p_hist | p_trace


mcmc_hist(posterior::subset_draws(draws, "delay_mean")) +
  theme_bw() +
  
  scale_x_continuous(breaks = 1:21)
