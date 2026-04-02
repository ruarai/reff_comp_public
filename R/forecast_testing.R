


source("R/functions.R")


round_id <- "2025-05-08"

x_location <- "VIC"
x_pathogen <- "RSV"

year <- as.character(year(round_id))
combined_counts <- read_csv(str_c("../forecast-hub/processed-data/", year, "/", round_id, "/combined-count-", round_id, ".csv"))

eval_counts <- read_csv("data/combined-count-smoothed-2026-01-15.csv") |> 
  filter(notification_date > origin_date, location == x_location, pathogen == x_pathogen)


origin_date <- get_origin_date(combined_counts, x_pathogen, x_location)

n_days_ahead <- 365

training_counts <- process_counts(combined_counts, x_location, x_pathogen, n_days_ahead) |> 
  # Save some extra data in case we need it later
  mutate(origin_date = origin_date)

data_list <- create_stan_data(training_counts, n_days_ahead, x_pathogen, x_location)


library(cmdstanr)
model <- cmdstan_model("src/gp_renewal.stan")

data_list$k <- 8
data_list$init_factor <- 0.5e-2


data_list$rel_hum <- mean_weather$rel_hum[,1]
data_list$abs_hum <- mean_weather$abs_hum[,1]

fit <- model$optimize(
  data = data_list, seed = 5, iter = 1e4, init = 0.1
)
fit <- model$pathfinder(
  data = data_list, seed = 8, init = 0.01
)
fit <- model$variational(
  data = data_list, seed = 10, iter = 3e4,
  algorithm="fullrank"
)

fit <- model$sample(
  data = data_list,
  seed = 1,
  chains = 4,
  parallel_chains = 4,
  refresh = 50,
  iter_warmup = 1000, iter_sampling = 1000
)


draws <- fit |>
  tidybayes::tidy_draws() |> 
  select(-starts_with("cases"))

pars <- c(
  "scaling_annual", "scaling_noise", "ls_noise", "ls_annual", "phi",
  "delay_mean", "delay_sd", "log_R_0"
)

p_hist <- mcmc_hist(draws |> mutate(chain = .chain), pars = pars) +
  theme_bw()

p_trace <- mcmc_trace(draws |> mutate(chain = .chain), pars = pars) +
  theme_bw()

p_hist | p_trace

mcmc_trace(draws |> mutate(chain = .chain),
           pars = c("alpha", "rel_hum_factor", "abs_hum_factor", "immunity_half_life")) +
  theme_bw()

mcmc_pairs(draws |> mutate(chain = .chain),
           pars = c("alpha", "rel_hum_factor", "abs_hum_factor", "immunity_half_life"))




trajectory_data_imm <- fit |> 
  spread_draws(accumulated_immunity[ix_t, k], alpha) |> 
  mutate(date = training_counts$date[ix_t])

ggplot(trajectory_data_imm |> group_by(date, .draw) |> summarise(y = sum(accumulated_immunity) * alpha)) + 
  fastribbon::stat_fastribbon(aes(x = date, y = y)) +
  
  scale_fill_brewer()


ggplot(trajectory_data_imm) + 
  fastribbon::stat_fastribbon(aes(x = date, y = accumulated_immunity * alpha)) +
  
  facet_wrap(~k)



trajectory_data <- fit |> 
  spread_draws(log_R[ix_t], R_t[ix_t], inf[ix_t], inf_delayed[ix_t],
               cases[ix_t], f_annual[ix_t], f_noise[ix_t],
               latent_cases_log[ix_t]) |> 
  mutate(date = training_counts$date[ix_t],
         latent_cases = exp(latent_cases_log),
         
         cases = replace_na(cases, 0))


school_holidays <- get_school_holidays() |> 
  filter(location == x_location)

public_holidays <- get_public_holidays() |> 
  filter(location == x_location)


ggplot(trajectory_data) + 
  fastribbon::stat_fastribbon(aes(x = date, y = R_t))  + 
  fastribbon::stat_fastribbon(aes(x = date, y = exp(log_R)),
                              alpha = 0.2) +
  
  scale_fill_brewer()


ggplot(trajectory_data) + 
  fastribbon::stat_fastribbon(aes(x = date, y = latent_cases),
                              linewidth = 1) +
  geom_line(aes(x = notification_date, y = cases_trend), eval_counts,
            colour = "orangered") + 
  
  scale_fill_brewer()


p_reff <- plot_reff_components(trajectory_data, origin_date, school_holidays, public_holidays) +
  xlim(c(ymd(min(trajectory_data$date)), ymd(origin_date + days(365))))

p_latent_inf <- plot_latent_inf(trajectory_data, training_counts, origin_date, school_holidays, public_holidays) +
  geom_line(aes(x = notification_date, y = cases_trend), eval_counts,
            colour = "orangered") +
  
  xlim(c(ymd(min(trajectory_data$date)), ymd(origin_date + days(365))))

p_reff / p_latent_inf




