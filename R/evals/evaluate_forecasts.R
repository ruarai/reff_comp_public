

library(tidyverse)
library(arrow)

source("R/functions.R")

forecasts_new <- open_dataset("forecast_outputs_raw/") |> 
  mutate(model = "new")
forecasts_old <- open_dataset("forecast_outputs_old/uom-reff-composure-raw/") |> 
  mutate(model = "old")

eval_counts <- read_csv("data/combined-count-2026-01-15.csv")


location_pathogens <- forecasts_new |> distinct(location, pathogen) |> collect()

calculate_crps <- function(x, y) {
  x <- sort(x)
  m <- length(x)
  (2 / m) * mean((x - y) * (m * (y < x) - seq_len(m) + 0.5))
}


perf_data <- location_pathogens |> 
  arrange(pathogen, location) |> 
  rename(x_location = location, x_pathogen = pathogen) |> 
  pmap(
    function(x_location, x_pathogen) {
      
      forecast_data <- bind_rows(
        forecasts_new |> 
          filter(target == "case_incidence",
                 pathogen == x_pathogen, location == x_location) |> 
          collect(),
        forecasts_old |> 
          filter(target == "case_incidence",
                 pathogen == x_pathogen, location == x_location) |> 
          collect()
      )
      
      round_data <- forecast_data |> 
        mutate(notification_date = origin_date + horizon) |> 
        left_join(eval_counts)
      
      perf_data <- round_data |> 
        mutate(log_val = log(value + 1), log_truth = log(cases + 1)) |> 
        group_by(model, location, pathogen, round_id, origin_date, notification_date) |> 
        summarise(crps_log = calculate_crps(log_val, log_truth[[1]]), .groups = "drop")
      
    },
    .progress = TRUE
  ) |> 
  bind_rows()

perf_data |> 
  filter(notification_date <= ymd("2025-10-01")) |> 
  group_by(model, location, pathogen, horizon = as.numeric(notification_date - origin_date)) |> 
  summarise(crps_log = mean(crps_log)) |> 

  ggplot() +
  geom_line(aes(x = horizon, y = crps_log, colour = model)) +
  facet_grid(rows = vars(location), cols = vars(pathogen)) +
  
  theme_bw()


weeks_ahead <- 11

perf_data |> 
  pivot_wider(names_from = model, values_from = crps_log) |> 
  mutate(skill_score = (new - old) / (0 -old)) |> 
  mutate(horizon = as.numeric(notification_date - origin_date)) |> 
  filter(horizon > 7 * weeks_ahead, horizon <= 7 * (weeks_ahead + 1)) |> 
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

