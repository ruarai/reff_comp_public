

perf_files <- fs::dir_ls("results", glob = "*_performance.parquet", recurse = TRUE)

length(perf_files)

perf_data <- arrow::open_dataset(perf_files)


perf_data |> 
  collect() |> 
  mutate(horizon = as.numeric(date - origin_date)) |> 
  filter(horizon > 21, horizon <= 28) |> 
  group_by(location, pathogen, round_id) |> 
  summarise(crps_log = mean(crps_log)) |> 
  
  ggplot() +
  geom_point(aes(y = location, x = crps_log)) +
  
  facet_wrap(~pathogen, ncol = 1) +
  
  theme_bw()

perf_data |> 
  collect() |> 
  mutate(horizon = as.numeric(date - origin_date)) |> 
  filter(horizon > 21, horizon <= 28) |> 
  group_by(location, pathogen, round_id) |>
  summarise(elpd = mean(elpd)) |>
  
  ggplot() +
  geom_point(aes(y = location, x = elpd)) +
  
  facet_wrap(~pathogen, ncol = 1) +
  
  theme_bw()


perf_data |> 
  collect() |> 
  mutate(horizon = as.numeric(date - origin_date)) |> 
  filter(horizon > 0,  horizon <= 7 * 8) |> 
  
  ggplot() +
  geom_point(aes(x = horizon, y = crps_log),
             size = 0.5, stroke = 0.5) +
  
  facet_grid(rows = vars(location), cols = vars(pathogen)) +
  
  theme_bw()



perf_data |> 
  collect() |> 
  mutate(horizon = as.numeric(date - origin_date)) |> 
  filter(horizon > 21, horizon <= 28) |> 
  
  ggplot() +
  geom_linerange(aes(x = date, ymin = 0, ymax = crps_log),
                 size = 0.5, linewidth = 0.5) +
  
  facet_grid(rows = vars(location), cols = vars(pathogen)) +
  
  coord_cartesian(ylim = c(0, 1)) +
  
  theme_bw()
