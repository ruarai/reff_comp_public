

library(tidyverse)
library(arrow)

source("R/functions.R")

forecasts_A <- open_dataset("forecast_outputs/baseline/") |> 
  mutate(model = "baseline")
forecasts_B <- open_dataset("forecast_outputs/no_annual/") |> 
  mutate(model = "no_annual")



eval_counts <- read_csv("data/combined-count-2026-01-15.csv")

locations <- forecasts_A |> distinct(location) |> collect()

calculate_crps <- function(x, y) {
  x <- sort(x)
  m <- length(x)
  (2 / m) * mean((x - y) * (m * (y < x) - seq_len(m) + 0.5))
}


perf_data <- locations |> 
  arrange(location) |> 
  rename(x_location = location) |> 
  pmap(
    function(x_location) {
      
      forecast_data <- bind_rows(
        forecasts_A |> 
          filter(target == "case_incidence", location == x_location) |> 
          collect(),
        forecasts_B |> 
          filter(target == "case_incidence", location == x_location) |> 
          collect()
      )
      
      round_data <- forecast_data |> 
        mutate(notification_date = origin_date + horizon) |> 
        left_join(eval_counts)
      
      perf_data <- round_data |> 
        mutate(log_val = log(value + 1), log_truth = log(cases + 1)) |> 
        group_by(model, location, pathogen, round_id, origin_date, notification_date) |> 
        summarise(crps_log = calculate_crps(log_val, log_truth[[1]]),
                  crps = calculate_crps(value, cases[[1]]), .groups = "drop")
      
    },
    .progress = TRUE
  ) |> 
  bind_rows()


perf_data |> 
  mutate(horizon = as.numeric(notification_date - origin_date)) |> 
  filter(location != "NZ", location != "SA", horizon <= 28) |> 
  group_by(model, pathogen, location, horizon) |> 
  summarise(crps_log = mean(crps_log),
            crps = mean(crps)) |> 
  
  ggplot() +
  geom_line(aes(x = horizon, y = crps_log, colour = model)) +
  
  facet_grid(rows = vars(location), cols = vars(pathogen), scales = "free_y") +
  
  # scale_y_log10() +
  
  theme_bw()


weeks_ahead <- 3



perf_data |> 
  select(-crps_log) |> 
  pivot_wider(names_from = model, values_from = crps) |> 
  mutate(skill_score = (baseline - no_annual) / (0 -no_annual)) |> 
  mutate(horizon = as.numeric(notification_date - origin_date)) |> 
  filter(horizon > 7 * weeks_ahead, horizon <= 7 * (weeks_ahead + 1),
         location != "NZ", location != "SA") |> 
  ungroup() |> 
  slice_sample(prop = 1.0) |> 
  ggplot() +
  
  geom_linerange(aes(x = notification_date, ymin = 0, ymax = skill_score,
                     colour = skill_score > 0),
                 size = 0.5, stroke = 0.5) +
  
  scale_colour_manual(values = c("TRUE" = ggokabeito::palette_okabe_ito(3),
                                 "FALSE" = ggokabeito::palette_okabe_ito(1))) +
  
  coord_cartesian(ylim = c(-1, 1)) +
  
  facet_grid(rows = vars(location), cols = vars(pathogen)) +
  
  # scale_y_log10() +
  
  theme_bw()




perf_data |> 
  select(-crps) |> 
  pivot_wider(names_from = model, values_from = crps_log) |>
  mutate(horizon = as.numeric(notification_date - origin_date)) |> 
  filter(horizon <= 21, location != "NZ", location != "SA") |> 
  group_by(location, pathogen) |> 
  
  drop_na(no_annual) |> 
  
  summarise(skill_score = 1 - sum(baseline) / sum(no_annual)) |> 
  ggplot() +
  geom_vline(xintercept = 0) +
  geom_linerange(aes(y = location, xmin = 0, xmax = skill_score, colour = skill_score > 0),
                 linewidth = 2.0) +
  
  scale_colour_manual(values = c("TRUE" = ggokabeito::palette_okabe_ito(3),
                                 "FALSE" = ggokabeito::palette_okabe_ito(1))) +
  
  facet_wrap(~pathogen, ncol = 3) +
  
  ylab(NULL) + xlab("Skill score") +
  
  coord_cartesian(xlim = c(-0.2, 0.2)) +
  
  theme_bw() +
  theme(legend.position = "none")
  

