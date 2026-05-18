


source("R/functions.R")


round_id <- "2025-05-08"

x_location <- "TAS"
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

get_normal_params <- function(q, prob = 0.99) {
  p1 <- (1 - prob) / 2
  p2 <- 1 - (1 - prob) / 2
  
  z1 <- qnorm(p1)
  z2 <- qnorm(p2)
  
  sd_val <- (q[2] - q[1]) / (z2 - z1)
  mean_val <- q[1] - (sd_val * z1)
  
  return(c(mean = mean_val, sd = sd_val))
}

peak_log_size_range <- log(c(7, 20))
peak_timing_range <- as.numeric(ymd(c("2026-03-01"), ymd("2026-10-01")) - ymd("2026-01-01")) + 1

normal_log_peak_size <- get_normal_params(peak_log_size_range)
normal_peak_timing <- get_normal_params(peak_timing_range)

peak_prior_samples <- tibble(
  day_of_year = rnorm(1000, normal_peak_timing[1], normal_peak_timing[2]),
  peak_size = exp(rnorm(1000, normal_log_peak_size[1], normal_log_peak_size[2]))
)

ggplot(training_counts) +
  geom_line(aes(x = date, y = y_trend)) +
  
  geom_density2d(aes(x = ymd("2025-01-01") + day_of_year - 1, y = peak_size),
                 peak_prior_samples) +
  
  theme_bw()

data_list$log_peak_size_mean <- normal_log_peak_size[1]
data_list$log_peak_size_sd <- normal_log_peak_size[2]
data_list$peak_day_mean <- normal_peak_timing[1]
data_list$peak_day_sd <- normal_peak_timing[2]

data_list$peak_season_start <- which(training_counts$date == ymd("2025-01-01"))
data_list$peak_season_end <- which(training_counts$date == ymd("2026-01-01"))

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
  data = data_list, seed = 5, iter = 1e4
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
  geom_line(aes(x = date, y = latent_cases, group = .draw),
            alpha = 0.2) +
  
  geom_density2d(aes(x = ymd("2025-01-01") + day_of_year - 1, y = peak_size),
                 colour = "lightblue",
                 peak_prior_samples) +
  
  geom_density2d(aes(x = peak_date, y = peak_case_incidence),
                 peak_posterior_samples) +
  
  geom_line(aes(x = notification_date, y = value_trend),
            eval_counts, colour = "orangered") +
  
  scale_fill_brewer() +
  
  theme_bw() +
  
  coord_cartesian(ylim = c(0, 50))


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
  
  geom_line(aes(x = notification_date, y = value_trend), eval_counts, colour = "orangered") +
  
  labs(caption = job_name)

p_reff / p_latent_inf




