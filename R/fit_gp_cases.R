library(targets)
library(tidyverse)

source("R/functions.R")




combined_counts <- read_csv("../autogpcast/data/combined_counts.csv", show_col_types = FALSE) |> 
  filter(round_id == max(round_id)) |> 
  select(-round_id)

training_counts <- make_training_counts(combined_counts, "VIC", "RSV")

training_date <- ymd("2025-03-01")
n_days_training <- which(training_counts$date == training_date)

ggplot(training_counts) +
  geom_point(aes(x = date, y = y), size = 0.5, stroke = 0.5) +
  geom_line(aes(x = date, y = y_smooth)) +
  
  geom_vline(xintercept = training_date, linetype = "44") +
  
  theme_bw()

library(cmdstanr)

model <- cmdstan_model("src/gp_cases.stan")

# TODO: how to pick these M values?
# see paper on HSGP
data_list <- list(
  x = training_counts$x, y = training_counts$y[1:n_days_training],
  N = nrow(training_counts), N_data = n_days_training,
  M_annual = 20, M_weekly = 20, M_noise = 30, c_factor = 1.2
)

fit <- model$sample(
  data = data_list,
  seed = 1,
  chains = 4,
  parallel_chains = 4,
  refresh = 50
)

library(bayesplot)
library(tidybayes)

draws <- fit |> 
  tidybayes::tidy_draws()

mcmc_hist(draws, pars = c("var_annual", "var_weekly", "var_noise", "ls_noise", "ls_annual", "ls_weekly", "mu", "phi"))
mcmc_pairs(draws, pars = c("var_annual", "var_weekly", "var_noise", "ls_noise", "ls_annual", "ls_weekly", "mu", "phi"))


plot_data <- fit |> 
  spread_draws(f_total[ix_x], f_annual[ix_x], f_noise[ix_x], mu) |> 
  mutate(date = training_counts$date[ix_x])


plot_data |> 
  filter(.draw < 10) |> 
  ggplot() +
  geom_line(aes(x = date, y = f_total, group = .draw)) +
  geom_point(aes(x = date, y = log(y)),
             training_counts)


plot_data |> 
  filter(.draw %% 100 == 0) |>
  ggplot() +
  geom_line(aes(x = date, y = exp(f_total), group = .draw),
            alpha = 0.1) +
  geom_point(aes(x = date, y = y),
             training_counts) +
  
  coord_cartesian(xlim = c(training_date - days(7 * 52), NA)) +
  
  geom_vline(xintercept = training_date, linetype = "44") +
  
  theme_bw()

plot_data |> 
  filter(date >= training_date - days(7 * 4),
         date < training_date + days(7 * 8)) |>
  ggplot() +
  ggdist::stat_lineribbon(aes(x = date, y = exp(f_total)), linewidth = 0.0) +
  geom_point(aes(x = date, y = y),
             training_counts) +
  
  coord_cartesian(xlim = c(training_date - days(7 * 4), training_date + days(7 * 12))) +
  
  scale_fill_brewer() +
  
  geom_vline(xintercept = training_date, linetype = "44") +
  
  theme_bw()

plot_data |> 
  filter(.draw %% 100 == 0) |>
  ggplot() +
  geom_line(aes(x = date, y = exp(f_annual + mu), group = .draw),
            alpha = 0.1, colour = ggokabeito::palette_okabe_ito(5)) +
  
  geom_line(aes(x = date, y = exp(f_annual + f_noise + mu), group = .draw),
            alpha = 0.1) +
  
  geom_point(aes(x = date, y = y_smooth),
             size = 0.5, stroke = 0.5,
             training_counts, colour = "red") +
  
  
  geom_vline(xintercept = training_date, linetype = "44") +
  
  theme_bw()
