library(targets)
library(tidyverse)

source("R/functions.R")


combined_counts <- read_csv("../autogpcast/data/combined_counts.csv", show_col_types = FALSE) |> 
  filter(round_id == max(round_id)) |> 
  select(-round_id)

training_counts <- make_training_counts(combined_counts, "VIC", "SARSCOV2", days_ahead = 365)


eval_counts <- read_csv("../autogpcast/data/combined-count-2025-12-25.csv") |> 
  group_by(location, pathogen) |> 
  mutate(cases_trend = zoo::rollmean(cases, k = 7, fill = NA),
         outlier = FALSE) |> 
  make_training_counts(x_location = "VIC", x_pathogen = "SARSCOV2")

training_date <- ymd("2025-10-06")
n_days_training <- which(training_counts$date == training_date)

ggplot(training_counts) +
  geom_point(aes(x = date, y = y), size = 0.5, stroke = 0.5) +
  geom_line(aes(x = date, y = y_smooth)) +
  
  geom_vline(xintercept = training_date, linetype = "44") +
  
  theme_bw()

# gi <- make_gi(12.76, 1.7, 0.995) # RSV
gi <- make_gi(0.89, rate = 0.27, 0.995) # SARSCOV2

library(cmdstanr)

model <- cmdstan_model("src/gp_renewal.stan")

# TODO: how to pick these M values?
# see paper on HSGP
data_list <- list(
  x = training_counts$x, y = training_counts$y[1:n_days_training],
  N = nrow(training_counts), N_data = n_days_training,
  M_annual = 20, M_noise = 30, c_factor = 1.2,
  
  gi = gi, N_gi = length(gi)
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

mcmc_hist(draws, pars = c("scaling_annual", "scaling_weekly", "scaling_noise", "ls_noise", "ls_annual", "r_0", "phi"))


fit |> 
  spread_draws(obs_offset_weekly[ix_x]) |> 
  ggplot() +
  geom_point(aes(x = ix_x, y = obs_offset_weekly))

plot_data <- fit |> 
  spread_draws(log_R[ix_x], inf[ix_x], inf_weekly[ix_x], f_annual[ix_x], f_noise[ix_x]) |> 
  mutate(date = training_counts$date[ix_x])

plot_data |> 
  filter(.draw %% 100 == 0) |> 
  ggplot() +
  geom_line(aes(x = date, y = inf_weekly, group = .draw), alpha = 0.1) +
  geom_point(aes(x = date, y = y), size = 0.5, stroke = 0.5,
             training_counts) +
  
  geom_line(aes(x = date, y = y), colour = "orangered",
            eval_counts |> filter(date > training_date)) +
  
  theme_bw() +
  
  
  geom_vline(xintercept = training_date, linetype = "44") +
  
  coord_cartesian(xlim = c(training_date - days(7 * 12), training_date + days(7 * 4)),
                  ylim = c(0, 350))

plot_data |> 
  ungroup() |> 
  filter(.draw %% 100 == 0) |>
  select(date, .draw, f_annual, f_noise, log_R) |> 
  pivot_longer(-c(date, .draw)) |> 
  ggplot() +
  geom_line(aes(x = date, y = exp(value), group = .draw, colour = name),
            alpha = 0.2) +
  
  facet_wrap(~name, ncol = 1) +
  
  ggokabeito::scale_colour_okabe_ito(order = c(3, 5, 9)) +
  
  theme_bw()
  
plot_data |> 
  filter(.draw %% 100 == 0) |>  
  ggplot() +
  
  geom_line(aes(x = date, y = exp(f_annual), group = .draw),
            alpha = 0.1, colour = ggokabeito::palette_okabe_ito(3)) +
  
  geom_line(aes(x = date, y = exp(f_noise), group = .draw),
            alpha = 0.1, colour = ggokabeito::palette_okabe_ito(5)) +
  
  geom_line(aes(x = date, y = exp(log_R), group = .draw),
            alpha = 0.2, colour = ggokabeito::palette_okabe_ito(9)) +
  
  
  geom_vline(xintercept = training_date, linetype = "44") +
  
  theme_bw()


plot_quants <- plot_data |> 
  ungroup() |> 
  select(.draw, date, inf) |> 
  pivot_wider(names_from = .draw,
              names_prefix = "sample_",
              values_from = inf) |> 
  make_intervals()


alpha_vals <- scales::rescale(rev(1/1.7^(1:8)), to = c(0.05, 0.99))


ggplot(plot_quants) +
  geom_ribbon(aes(x = date, ymin = lower, ymax = upper,
                  group = quant, alpha = quant),
              fill = ggokabeito::palette_okabe_ito(5)) +
  
  geom_vline(xintercept = training_date) +
  
  geom_line(aes(x = date, y = y_smooth),
            training_counts) +
  
  geom_line(aes(x = date, y = y_smooth), colour = "orangered",
            eval_counts |> filter(date > training_date)) +
  
  
  xlab("Notification date") + ylab("Count") +
  
  scale_alpha_manual(values = alpha_vals) +
  
  scale_y_log10() +
  # coord_cartesian(ylim = c(NA, 500)) +
  
  theme_bw() +
  
  theme(strip.text = element_blank(),
        legend.position = "none")

plot_quants_weekly <- plot_data |> 
  ungroup() |> 
  select(.draw, date, inf_weekly) |> 
  pivot_wider(names_from = .draw,
              names_prefix = "sample_",
              values_from = inf_weekly) |> 
  make_intervals()

ggplot(plot_quants_weekly) +
  geom_ribbon(aes(x = date, ymin = lower, ymax = upper,
                  group = quant, alpha = quant),
              fill = ggokabeito::palette_okabe_ito(5)) +
  
  geom_vline(xintercept = training_date) +
  
  geom_line(aes(x = date, y = y),
            training_counts) +
  
  geom_line(aes(x = date, y = y), colour = "orangered",
            eval_counts |> filter(date > training_date)) +
  
  
  xlab("Notification date") + ylab("Count") +
  
  scale_alpha_manual(values = alpha_vals) +
  
  coord_cartesian(ylim = c(0, 500),
                  xlim = c(training_date - days(7 * 12), training_date + days(7 * 20))
                  ) +
  
  theme_bw() +
  
  theme(strip.text = element_blank(),
        legend.position = "none")
