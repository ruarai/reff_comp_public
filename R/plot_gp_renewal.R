library(bayesplot)
library(tidybayes)
library(tidyverse)


source("R/functions.R")

x_location <- "NSW"
x_pathogen <- "RSV"

fit <- read_rds(str_c("results/2025-06-19/", x_pathogen, "_", x_location, ".rds"))
training_counts <- arrow::read_parquet(str_c("results/2025-06-19/", x_pathogen, "_", x_location, "_training_counts.parquet"))



eval_counts <- read_csv("data/combined-count-2025-12-25.csv") |> 
  group_by(location, pathogen) |> 
  mutate(cases_trend = zoo::rollmean(cases, k = 7, fill = NA),
         outlier = FALSE) |> 
  filter(location == x_location, pathogen == x_pathogen) |> 
  make_training_counts(x_location, x_pathogen)

origin_date <- training_counts |> 
  drop_na(y) |> 
  pull(date) |> 
  max()

draws <- fit |> 
  tidybayes::tidy_draws()

mcmc_hist(draws, pars = c("scaling_annual", "scaling_weekly", "scaling_noise", "ls_noise", "ls_annual", "r_0", "phi"))


plot_data <- fit |> 
  spread_draws(log_R[ix_x], inf[ix_x], inf_weekly[ix_x], f_annual[ix_x], f_noise[ix_x]) |> 
  mutate(date = training_counts$date[ix_x])

plot_data |> 
  filter(.draw %% 100 == 0) |> 
  ggplot() +
  geom_line(aes(x = date, y = inf, group = .draw), alpha = 0.1) +
  geom_point(aes(x = date, y = y_smooth), size = 0.5, stroke = 0.5,
             training_counts) +
  
  geom_line(aes(x = date, y = y_smooth), colour = "orangered",
            eval_counts |> filter(date >= ymd("2025-06-12"))) +
  
  theme_bw() +
  
  coord_cartesian(ylim = c(0, 500))

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
  
  geom_vline(xintercept = origin_date) +
  
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
  
  geom_vline(xintercept = origin_date) +
  
  geom_line(aes(x = date, y = y_smooth),
            training_counts) +
  
  geom_line(aes(x = date, y = y_smooth), colour = "orangered",
            eval_counts |> filter(date > origin_date)) +
  
  
  xlab("Notification date") + ylab("Count") +
  
  scale_alpha_manual(values = alpha_vals) +
  
  # scale_y_log10() +
  coord_cartesian(ylim = c(NA, 500)) +
  
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
  
  geom_vline(xintercept = origin_date) +
  
  geom_line(aes(x = date, y = y),
            training_counts) +
  
  geom_line(aes(x = date, y = y), colour = "orangered",
            eval_counts |> filter(date > origin_date)) +
  
  
  xlab("Notification date") + ylab("Count") +
  
  scale_alpha_manual(values = alpha_vals) +
  
  coord_cartesian(ylim = c(0, 500),
                  xlim = c(origin_date - days(7 * 12), origin_date + days(7 * 20))
  ) +
  
  theme_bw() +
  
  theme(strip.text = element_blank(),
        legend.position = "none")

