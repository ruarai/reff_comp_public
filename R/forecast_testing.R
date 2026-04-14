


source("R/functions.R")


round_id <- "2025-05-08"

x_location <- "QLD"
x_pathogen <- "RSV"

year <- as.character(year(round_id))
combined_counts <- read_csv(str_c("../forecast-hub/processed-data/", year, "/", round_id, "/combined-count-", round_id, ".csv"))

origin_date <- get_origin_date(combined_counts, x_pathogen, x_location)

eval_counts <- read_csv("data/combined-count-smoothed-2026-01-15.csv") |> 
  filter(notification_date > origin_date, location == x_location, pathogen == x_pathogen)



n_days_ahead <- 365

training_counts <- process_counts(combined_counts, x_location, x_pathogen, n_days_ahead) |> 
  # Save some extra data in case we need it later
  mutate(origin_date = origin_date)

data_list <- create_stan_data(training_counts, n_days_ahead, x_pathogen, x_location)


library(cmdstanr)
model <- cmdstan_model("src/gp_renewal_immunity.stan")

fit <- model$optimize(
  data = data_list, seed = 3, iter = 1e4, init = 0.1
)
fit <- model$pathfinder(
  data = data_list, seed = 6, init = 0.5
)
fit <- model$variational(
  data = data_list, seed = 10, iter = 5e4,
  init = fit
)

fit <- model$sample(
  data = data_list,
  seed = 1,
  chains = 4,
  parallel_chains = 4,
  refresh = 50,
  iter_warmup = 1000, iter_sampling = 1000
  # init = fit
)

stan_pairplot(fit, "scaling_humidity", c(pars, "r_0", "initial_immunity_prop[1]", "alpha", "immunity_half_life"))


draws <- fit |>
  tidybayes::tidy_draws() |> 
  select(-starts_with("cases"))

pars <- c(
  "scaling_noise", "ls_noise", "phi", "scaling_humidity", 
  "delay_mean", "delay_sd"
)

p_hist <- mcmc_hist(draws |> mutate(chain = .chain), pars = pars) +
  theme_bw()

p_hist

p_trace <- mcmc_trace(draws |> mutate(chain = .chain), pars = pars) +
  theme_bw()
p_trace

mcmc_trace(draws |> mutate(chain = .chain),
           pars = c("alpha", "immunity_half_life")) +
  theme_bw()

mcmc_pairs(draws |> mutate(chain = .chain),
           pars = c("alpha", "immunity_half_life"))


# mcmc_pairs(draws |> mutate(chain = .chain), pars = pars)



trajectory_data_imm <- fit|> 
  spread_draws(immune_effect[ix_t], alpha, humidity_effect[ix_t]) |> 
  mutate(date = training_counts$date[ix_t],
         )

ggplot(trajectory_data_imm) + 
  fastribbon::stat_fastribbon(aes(x = date, y = immune_effect)) +
  
  scale_fill_brewer()

ggplot(trajectory_data_imm) + 
  fastribbon::stat_fastribbon(aes(x = date, y = humidity_effect)) +
  
  scale_fill_brewer()

data_list$humidity


trajectory_data <- fit |> 
  spread_draws(log_R[ix_t], R_t[ix_t], inf[ix_t], inf_delayed[ix_t],
               cases[ix_t], f_noise[ix_t],
               latent_cases_log[ix_t],
               # log_R_0
               ) |> 
  mutate(date = training_counts$date[ix_t],
         latent_cases = exp(latent_cases_log),
         # log_R_0 = cut(log_R_0, 4),
         
         cases = replace_na(cases, 0))


school_holidays <- get_school_holidays() |> 
  filter(location == x_location)

public_holidays <- get_public_holidays() |> 
  filter(location == x_location)


ggplot(trajectory_data) + 
  fastribbon::stat_fastribbon(aes(x = date, y = R_t))  + 
  fastribbon::stat_fastribbon(aes(x = date, y = exp(log_R)),
                              alpha = 0.2) +
  
  # facet_wrap(~log_R_0) +
  
  scale_fill_brewer()

trajectory_data |> 
  filter(.draw %% 222 == 0) |>
  ggplot() +
  geom_line(aes(x = date, y = latent_cases, group = .draw)) +
  geom_line(aes(x = notification_date, y = cases_trend), eval_counts,
            colour = "orangered") +
  geom_line(aes(x = date, y = y_smooth), training_counts |> mutate(y_smooth = zoo::rollmean(y, 7, fill = NA)),
            colour = "blue")

ggplot(trajectory_data) + 
  fastribbon::stat_fastribbon(aes(x = date, y = latent_cases),
                              linewidth = 1) +
  geom_line(aes(x = date, y = y_smooth), training_counts |> mutate(y_smooth = zoo::rollmean(y, 7, fill = NA)),
            colour = "blue") +
  geom_line(aes(x = notification_date, y = cases_trend), eval_counts,
            colour = "orangered") + 
  
  # facet_wrap(~log_R_0) +
  
  coord_cartesian(ylim = c(5, 300)) +
  
  scale_y_log10() +
  
  scale_fill_brewer()

