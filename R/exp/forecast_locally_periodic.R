
source("R/functions.R")

round_id <- "2026-04-29"

x_location <- "VIC"
x_pathogen <- "RSV"
x_target <- "case_incidence"

year <- as.character(year(round_id))
combined_counts <- read_csv(str_c("../forecast-hub/processed-data/", year, "/", round_id, "/combined-count-smoothed-", round_id, ".csv"))

origin_date <- get_origin_date(combined_counts, x_pathogen, x_location, x_target)

eval_counts <- read_csv("../forecast-hub/processed-data/2026/2026-05-20/combined-count-smoothed-2026-05-20.csv") |> 
  filter(notification_date > origin_date, location == x_location, pathogen == x_pathogen)

n_days_ahead <- 365

training_counts <- process_counts(combined_counts, x_location, x_pathogen, x_target, n_days_ahead, origin_date)

data_list <- create_stan_data(training_counts, n_days_ahead, x_pathogen, x_location)

data_list$M_annual <- 600
data_list$c_factor_annual <- 5

library(cmdstanr)
model <- cmdstan_model("src/gp_renewal_lp.stan")

fit <- model$optimize(
  data = data_list, seed = 3, iter = 2e4, init = 0.1,
  
  algorithm = "lbfgs",
  tol_obj     = 1e-15,
  tol_rel_obj = 1,    
  tol_param   = 1e-12,
  tol_grad    = 1e-12,
  tol_rel_grad = 1e4,
  jacobian = TRUE
)

fit <- model$pathfinder(
  data = data_list, seed = 3, init = 0.1
)
fit <- model$variational(
  data = data_list, seed = 6, iter = 1e4
)



fit <- model$sample(
  data = data_list,
  seed = 1,
  chains = 4,
  parallel_chains = 4,
  refresh = 50,
  iter_warmup = 250, iter_sampling = 250
)


pars <- c(
  "scaling_annual",
  "ls_annual", "ls_annual_decay_years"
)

draws <- fit$draws() |> 
  posterior::subset_draws(pars)

draws <- draws[,1:3,]

mcmc_hist(draws) +
  theme_bw()
mcmc_trace(draws) +
  theme_bw()


trajectory_data <- fit$draws()[,1:3,] |> 
  spread_draws(log_R[ix_t], inf[ix_t], inf_delayed[ix_t],
               cases[ix_t], f_annual[ix_t], f_noise[ix_t],
               latent_cases_log[ix_t]) |> 
  
  mutate(date = training_counts$date[ix_t],
         latent_cases = exp(latent_cases_log),
         
         cases = replace_na(cases, 0))


trajectory_data |> 
  filter((.draw-1) %% 100 == 0) |>
  ggplot() +
  geom_line(aes(x = date, y = f_annual, group = .draw),
            alpha = 0.2)

trajectory_data |> 
  filter((.draw-1) %% 100 == 0) |>
  ggplot() +
  geom_line(aes(x = date, y = latent_cases, group = .draw),
            alpha = 0.2) +
  
  geom_line(aes(x = notification_date, y = value_trend),
            eval_counts, colour = "orangered") +
  
  scale_fill_brewer() +
  
  theme_bw() +
  
  coord_cartesian(ylim = c(0, 500))


public_holidays <- get_public_holidays() |> 
  filter(location == x_location)
school_holidays <- get_school_holidays() |> 
  filter(location == x_location)

make_forecast_plot(
  trajectory_data, training_counts, origin_date,
  school_holidays, public_holidays, eval_counts,
  NULL
)



p_reff <- plot_reff_components(trajectory_data, origin_date, school_holidays, public_holidays) +
  xlim(c(ymd(min(trajectory_data$date)), ymd(origin_date + days(365))))

y_upper_latent <- trajectory_data |> 
  filter(date <= origin_date + days(26 * 7)) |> 
  group_by(date) |> summarise(q95 = quantile(latent_cases, 0.95)) |> summarise(q95 = max(q95)) |> pull(q95)

p_latent_inf <- plot_latent_inf(trajectory_data, training_counts, origin_date, school_holidays, public_holidays) +
  
  xlim(c(ymd(min(trajectory_data$date)), ymd(origin_date + days(365)))) +
  coord_cartesian(ylim = c(0, 1000)) +
  
  geom_line(aes(x = notification_date, y = value_trend), eval_counts, colour = "orangered")

p_reff / p_latent_inf


ggplot(trajectory_data) +
  
  fastribbon::stat_fastribbon(aes(x = yday(date), y = f_annual)) +
  
  facet_wrap(~year(date), ncol = 1) +
  
  theme_bw() +
  scale_fill_brewer()
  
  
ggplot(trajectory_data |> filter((.draw - 1) %% 100 == 0)) +
    
  geom_line(aes(x = yday(date), y = f_annual, colour = factor(.draw))) +
  
  facet_wrap(~year(date), ncol = 1) +

  theme_bw() +
  scale_fill_brewer()




