

library(tidyverse)
library(tidybayes)

x_location <- "VIC"
x_pathogen <- "RSV"

jobs <- expand_grid(
  x_location = c("VIC", "NSW", "QLD", "NT", "TAS", "SA", "WA", "ACT"),
  x_pathogen = c("RSV", "SARSCOV2", "flu")
)

plot_data <- jobs |> 
  pmap(
    function(x_location, x_pathogen) {
      fit_file <- fs::path("results", "2025-10-09", str_c(x_pathogen, "_", x_location, ".rds"))
      fit <- read_rds(fit_file)
      
      fit |> 
        spread_draws(school_holiday_Rt_effect[ix_holiday], school_holiday_asc_effect[ix_holiday],
                     delay_mean, delay_sd, scaling_annual, scaling_noise, scaling_weekly_asc, ls_annual,
                     ls_noise, phi) |> 
        filter(ix_holiday > 1) |> 
        summarise_draws() |> 
        mutate(location = x_location, pathogen = x_pathogen)
    },
    .progress = TRUE
  ) |> 
  bind_rows()



plot_data |> 
  bind_rows() |> 
  filter(variable == "school_holiday_asc_effect") |> 
  mutate(holiday_name = case_when(
    ix_holiday == 2 ~ "Summer holidays",
    ix_holiday == 3 ~ "Mid-year holidays"
  )) |> 
  ungroup()  |> 
  mutate(holiday_name = factor(holiday_name, c("Summer holidays", "Mid-year holidays"))) |> 
  ggplot() +
  geom_vline(xintercept = 1, linetype = "22",
             alpha = 0.5) +
  
  geom_linerange(aes(y = location, xmin = exp(q5), xmax = exp(q95), colour = holiday_name),
                 linewidth = 2.0,
                 position = position_dodge2(width = 0.5)) +
  facet_wrap(~pathogen) +
  
  scale_colour_brewer(palette = 6, type = "qual") +
  
  coord_cartesian(xlim = exp(c(-0.3, 0.3))) +
  scale_x_log10(breaks = c(0.5, 0.66, 0.8, 0.9, 1, 1.10, 1.25, 1.5, 2)) +
  
  theme_bw()


plot_data |> 
  bind_rows() |> 
  filter(variable == "delay_sd", ix_holiday == 2) |> 
  ggplot() +
  
  geom_linerange(aes(y = location, xmin = q5, xmax = q95),
                 linewidth = 2.0) +
  facet_wrap(~pathogen, ncol = 1) +
  
  scale_colour_brewer(palette = 6, type = "qual") +
  
  theme_bw()





plot_data_delay <- jobs |> 
  pmap(
    function(x_location, x_pathogen) {
      fit_file <- fs::path("results", "2025-10-09", str_c(x_pathogen, "_", x_location, ".rds"))
      fit <- read_rds(fit_file)
      
      fit |> 
        spread_draws(delay_mean, delay_sd)  |> 
        mutate(sigma = sqrt(log(1 + (delay_sd ^ 2) / (delay_mean ^ 2))),
               mu = log(delay_mean) - 0.5 * sigma ^ 2) |> 
        
        expand_grid(t = seq(1, 14, by = 1)) |> 
        mutate(p = plnorm(t, mu, sigma) - plnorm(t - 1, mu, sigma)) |> 
        mutate(location = x_location, pathogen = x_pathogen)
    },
    .progress = TRUE
  ) |> 
  bind_rows()


plot_data_delay |> 
  filter(.draw %% 200 == 0) |> 
  ggplot() +
  geom_step(aes(x = t - 1, y = p, group = .draw),
            alpha = 0.1) +
  
  facet_grid(rows = vars(location), cols = vars(pathogen)) +
  
  theme_bw()




fit |> 
  spread_draws(delay_mean, delay_sd)  |> 
  mutate(sigma = sqrt(log(1 + (delay_sd ^ 2) / (delay_mean ^ 2))),
         mu = log(delay_mean) - 0.5 * sigma ^ 2) |> 
  
  expand_grid(t = seq(1, 14, by = 1)) |> 
  mutate(p = plnorm(t, mu, sigma) - plnorm(t - 1, mu, sigma)) |> 
  mutate(location = x_location, pathogen = x_pathogen) |> 
  filter(.draw %% 200 == 0) |> 
  ggplot() +
  geom_step(aes(x = t - 1, y = p, group = .draw),
            alpha = 0.1) +
  
  facet_grid(rows = vars(location), cols = vars(pathogen)) +
  
  theme_bw()

