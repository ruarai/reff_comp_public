


source("R/functions.R")


round_id <- "2026-04-29"

x_location <- "VIC"
x_pathogen <- "RSV"
x_target <- "case_incidence"

year <- as.character(year(round_id))
combined_counts <- read_csv(str_c("../forecast-hub/processed-data/", year, "/", round_id, "/combined-count-smoothed-", round_id, ".csv"))

origin_date <- get_origin_date(combined_counts, x_pathogen, x_location, x_target)

eval_counts <- read_csv("../forecast-hub/processed-data/2026/2026-08-19/combined-count-smoothed-2026-08-19.csv") |> 
  filter(notification_date > origin_date, location == x_location, pathogen == x_pathogen)

n_days_ahead <- 365

weather_data <- read_csv("data/weather_data.csv")

training_counts <- process_counts(combined_counts, x_location, x_pathogen, x_target, n_days_ahead, origin_date) |> 
  left_join(weather_data)

data_list <- create_stan_data(training_counts, n_days_ahead, x_pathogen, x_location)

max_window_size <- 360

data_list$N_max_window <- max_window_size

data_list$past_infections <- combined_counts |> 
  filter(location == x_location, pathogen == x_pathogen, target == x_target) |> 
  complete(notification_date = seq(ymd("2020-01-01"), max(notification_date), "days"), fill = list(value_trend = 0)) |> 
  arrange(notification_date) |> 
  filter(notification_date >= min(training_counts$date) - days(max_window_size),
         notification_date < min(training_counts$date)) |> 
  pull(value_trend)

library(cmdstanr)
model <- cmdstan_model("src/gp_renewal_sus.stan")


fit <- model$sample(
  data = data_list,
  seed = 1,
  chains = 4,
  parallel_chains = 4,
  refresh = 50,
  iter_warmup = 500, iter_sampling = 500,
  init = 0.1
)

fit <- read_rds("saved_results/gp_renewal_sus/2026-04-29/RSV_VIC_case_incidence.rds")


draws <- fit$draws(c(
  "scaling_noise", "ls_noise", "phi", "ls_annual", "scaling_annual",
  "cumulative_infections_effect", "sus_window_size"
))

mcmc_hist(draws) +
  theme_bw()
mcmc_trace(draws) +
  theme_bw()

trajectory_data <- fit |> 
  spread_draws(log_R[ix_t], inf[ix_t], inf_delayed[ix_t],
               cases[ix_t], f_noise[ix_t],
               latent_cases_log[ix_t],
               cumulative_effect[ix_t],
  ) |> 
  mutate(date = training_counts$date[ix_t],
         latent_cases = exp(latent_cases_log),
         cases = replace_na(cases, 0)) |> 
  
  filter(.chain != 3)


wrap_plots(
  ggplot(trajectory_data) + 
    fastribbon::stat_fastribbon(aes(x = date, y = exp(log_R))) + 
    scale_fill_brewer() + 
    geom_hline(yintercept = 1) +
    
    coord_cartesian(ylim = c(0.5, 1.8)) +
    theme_bw(),
  ggplot(trajectory_data) + 
    fastribbon::stat_fastribbon(aes(x = date, y = exp(-cumulative_effect))) + 
    scale_fill_brewer() + 
    geom_hline(yintercept = 1) +
    
    coord_cartesian(ylim = c(0.5, 1.8)) +
    theme_bw(),
  ggplot(trajectory_data) + 
    fastribbon::stat_fastribbon(aes(x = date, y = exp(log_R - cumulative_effect))) + 
    scale_fill_brewer() + 
    geom_hline(yintercept = 1) +
    
    coord_cartesian(ylim = c(0.5, 1.8)) +
    theme_bw(),
  ncol = 1
)

ggplot(trajectory_data) + 
  fastribbon::stat_fastribbon(aes(x = date, y = latent_cases),
                              linewidth = 0) +
  geom_line(aes(x = date, y = y_trend), training_counts,
            colour = "black") +
  
  geom_line(aes(x = notification_date, y = value_trend), eval_counts,
            colour = "orangered") +
  
  coord_cartesian(ylim = c(0, 400)) +
  
  scale_fill_brewer() +
  
  theme_bw()


library(pammtools)
ggplot(trajectory_data) + 
  fastribbon::stat_fastribbon(aes(x = date, y = cases),
                              geom = "stepribbon") +
  
  geom_line(aes(x = date, y = cases, group = .draw),
            trajectory_data |> filter((.draw - 1) %% 2000 == 0),
            alpha = 0.1) +
  geom_line(aes(x = date, y = y), training_counts,
            colour = "black") +
  geom_line(aes(x = notification_date, y = value), eval_counts,
            colour = "orange") +
  
  coord_cartesian(ylim = c(0, 600), xlim = c(origin_date - days(28), origin_date + days(7 * 8))) +
  
  scale_fill_brewer() +
  
  theme_bw()


trajectory_data <- fit |> 
  spread_draws(log_R[ix_t], inf[ix_t], inf_delayed[ix_t],
               cases[ix_t], f_annual[ix_t], f_noise[ix_t],
               latent_cases_log[ix_t]) |> 
  
  mutate(date = training_counts$date[ix_t],
         latent_cases = exp(latent_cases_log),
         
         cases = replace_na(cases, 0)) |> 
  
  filter(.chain != 3)


school_holidays <- get_school_holidays() |> 
  filter(location == x_location, end_date >= ymd("2023-01-01"), start_date <= origin_date + 365)

public_holidays <- get_public_holidays() |> 
  filter(location == x_location, date >= ymd("2023-01-01"), date <= origin_date + 365)


p_reff <- plot_reff_components(trajectory_data, origin_date, school_holidays, public_holidays) +
  xlim(c(ymd(min(trajectory_data$date)), ymd(origin_date + days(365))))

y_upper_latent <- trajectory_data |> 
  filter(date <= origin_date + days(26 * 7)) |> 
  group_by(date) |> summarise(q95 = quantile(latent_cases, 0.95)) |> summarise(q95 = max(q95)) |> pull(q95)

p_latent_inf <- plot_latent_inf(trajectory_data, training_counts, origin_date, school_holidays, public_holidays) +
  
  xlim(c(ymd(min(trajectory_data$date)), ymd(origin_date + days(365)))) +
  coord_cartesian(ylim = c(0, y_upper_latent * 1.25)) +
  
  geom_line(aes(x = notification_date, y = value_trend_and_remainder), eval_counts, colour = "orangered")

p_reff / p_latent_inf
