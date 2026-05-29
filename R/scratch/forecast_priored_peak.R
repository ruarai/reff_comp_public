


source("R/functions.R")
source("R/peak_priors.R")


round_id <- "2026-05-27"

x_location <- "VIC"
x_pathogen <- "RSV"
x_target <- "case_incidence"

year <- as.character(year(round_id))
combined_counts <- read_csv(str_c("../forecast-hub/processed-data/", year, "/", round_id, "/combined-count-smoothed-", round_id, ".csv"))

origin_date <- get_origin_date(combined_counts, x_pathogen, x_location, x_target)

eval_counts <- read_csv("../forecast-hub/processed-data/2026/2026-04-22/combined-count-smoothed-2026-04-22.csv") |> 
  filter(notification_date > origin_date, location == x_location, pathogen == x_pathogen)

n_days_ahead <- 365

training_counts <- process_counts(combined_counts, x_location, x_pathogen, x_target, n_days_ahead, origin_date)

data_list <- create_stan_data(training_counts, n_days_ahead, x_pathogen, x_location)

peak_prior <- get_peak_priors() |> 
  filter(location == x_location, pathogen == x_pathogen)

data_list$log_peak_size_mean <- peak_prior$log_peak_size_mean
data_list$log_peak_size_sd <- peak_prior$log_peak_size_sd
data_list$peak_day_mean <- peak_prior$peak_yday_mean
data_list$peak_day_sd <- peak_prior$peak_yday_sd

data_list$peak_season_start <- which(training_counts$date == ymd("2026-01-01"))
data_list$peak_season_end <- which(training_counts$date == ymd("2027-01-01")) - 1


n_samples <- 1000
peak_prior_samples <- peak_prior |> 
  rowwise() |> 
  mutate(
    day_of_year = list(rnorm(n_samples, peak_yday_mean, peak_yday_sd)),
    peak_size = list(exp(rnorm(n_samples, log_peak_size_mean, log_peak_size_sd)))
  ) |> 
  unnest(c(day_of_year, peak_size))

ggplot(training_counts) +
  geom_line(aes(x = date, y = y_trend)) +
  
  geom_density2d(aes(x = ymd("2026-01-01") + day_of_year - 1, y = peak_size),
                 peak_prior_samples) +
  
  theme_bw()

library(cmdstanr)
model <- cmdstan_model("src/gp_renewal_peak.stan")

fit <- model$optimize(
  data = data_list, seed = 4, iter = 2e4, init = 0.1,
  
  algorithm = "lbfgs",
  tol_obj     = 1e-15,
  tol_rel_obj = 1,    
  tol_param   = 1e-12,
  tol_grad    = 1e-12,
  tol_rel_grad = 1e4,
  jacobian = TRUE
)

fit <- model$laplace(
  data = data_list, seed = 3, init = 0.1,
  mode = fit, jacobian = TRUE
)


fit <- model$pathfinder(
  data = data_list, seed = 4, init = 0.5
)
fit <- model$variational(
  data = data_list, seed = 4, iter = 1e4,
  algorithm = "fullrank"
)



fit <- model$sample(
  data = data_list,
  seed = 1,
  chains = 4,
  parallel_chains = 4,
  refresh = 50,
  iter_warmup = 250, iter_sampling = 250
)


trajectory_data <- fit |> 
  spread_draws(log_R[ix_t], inf[ix_t], inf_delayed[ix_t],
               cases[ix_t], f_annual[ix_t], f_noise[ix_t],
               latent_cases_log[ix_t]) |> 
  
  mutate(date = training_counts$date[ix_t],
         latent_cases = exp(latent_cases_log),
         
         cases = replace_na(cases, 0))

peak_posterior_samples <- trajectory_data |> 
  filter(year(date) == year(round_id)) |> 
  group_by(.draw) |> 
  summarise(peak_case_incidence = max(latent_cases),
            peak_date = first(date[latent_cases == max(latent_cases)]))

library(pammtools)

trajectory_data |> 
  filter(.draw %% 100 == 0) |>
  ggplot() +
  
  geom_density2d(aes(x = ymd("2026-01-01") + day_of_year - 1, y = peak_size),
                 colour = "lightblue",
                 peak_prior_samples) +
  
  geom_density2d(aes(x = peak_date, y = peak_case_incidence),
                 peak_posterior_samples) +
  
  geom_line(aes(x = notification_date, y = value_trend),
            eval_counts, colour = "orangered") +
  
  geom_line(aes(x = date, y = latent_cases, group = .draw),
            alpha = 0.2) +
  
  scale_fill_brewer() +
  
  theme_bw() +
  
  coord_cartesian(ylim = c(0, 400))


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
  coord_cartesian(ylim = c(0, y_upper_latent * 1.25)) +
  
  geom_line(aes(x = notification_date, y = value_trend), eval_counts, colour = "orangered")

p_reff / p_latent_inf




