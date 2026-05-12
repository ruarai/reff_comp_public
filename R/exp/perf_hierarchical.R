
source("R/functions.R")

source("R/presentation/plot_theme.R")

perf_data <- bind_rows(
  read_parquet("saved_results/results_baseline_performance.parquet") |> 
    mutate(.model = "baseline"),
  read_parquet("saved_results/hierarchical_performance.parquet") |> 
    mutate(.model = "hierarchical"),
) |> 
  filter(location %in% c("ACT", "NSW", "WA", "QLD", "VIC", "TAS")) |> 
  mutate(horizon = as.numeric(date - origin_date))

perf_data |> 
  group_by(round_id, pathogen, location, date) |> 
  filter(n() == 2) |> 
  ungroup() |> 
  filter(horizon > 21, horizon <= 7*4) |> 
  group_by(.model, pathogen) |> 
  summarise(crps_log = mean(crps_log))


perf_data |> 
  filter(horizon == 7) |> 
  ggplot() + 
  geom_tile(aes(x = ymd(round_id), y = location, fill = log(crps_log))) +
  
  facet_grid(cols = vars(.model), rows = vars(pathogen)) +
  
  scale_fill_viridis_c()

perf_data |> 
  group_by(round_id, pathogen, location, date) |> 
  filter(n() == 2) |> 
  ungroup() |> 
  filter(horizon > 0, horizon <= 7*8) |> 
  group_by(pathogen, location, .model, horizon) |> 
  summarise(crps_log = mean(crps_log)) |> 
  
  ggplot() +
  geom_line(aes(x = horizon, y = crps_log, colour = .model, group = .model)) +
  
  ggokabeito::scale_colour_okabe_ito(order = c(9, 2)) +
  
  facet_grid(rows = vars(location), cols = vars(pathogen), scales = "free_y") +
  
  theme_bw()

perf_data |> 
  group_by(round_id, pathogen, location, date) |> 
  filter(n() == 2) |> 
  ungroup() |> 
  filter(horizon > 0, horizon <= 7*8) |> 
  group_by(pathogen, location, .model) |> 
  summarise(crps_log = mean(crps_log)) |> 
  
  ggplot() +
  geom_linerange(aes(x = .model, ymin = 0, ymax = crps_log),
                 linewidth = 4.0) +
  
  facet_grid(rows = vars(location), cols = vars(pathogen), scales = "free_y") +
  
  theme_bw()


perf_data |> 
  group_by(round_id, pathogen, location, date) |> 
  filter(n() == 2) |> 
  ungroup() |> 
  filter(horizon > 21, horizon <= 7*4) |> 
  group_by(pathogen, location, .model, round_id) |> 
  summarise(crps_log = mean(crps_log)) |> 
  
  ggplot() +
  geom_line(aes(x = ymd(round_id), y = crps_log, colour = .model, group = .model)) +
  
  ggokabeito::scale_colour_okabe_ito(order = c(2, 9)) +
  
  facet_grid(rows = vars(location), cols = vars(pathogen), scales = "free_y") +
  
  theme_bw()

