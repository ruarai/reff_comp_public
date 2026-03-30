
source("R/functions.R")

scores <- open_dataset(fs::dir_ls("saved_results/results_baseline/", glob = "*energy_score.parquet", recurse = TRUE))

eval_data <- read_csv("../forecast-hub/processed-data/2026/2026-01-15/combined-count-2026-01-15.csv")

plot_data_cases <- eval_data |> 
  mutate(location = str_c(location, "_data"))

locations <- c("VIC", "NSW", "QLD")

scores |>
  collect() |>
  filter(location %in% locations,
         pathogen == "RSV") |> 
  mutate(horizon = as.numeric(date - origin_date),
         round_id = as.numeric((ymd(round_id) - ymd("2025-05-08"))) / as.numeric(weeks(1))) |> 
  filter(horizon > 0, horizon <= 28) |> 
  
  ggplot() +
  geom_tile(aes(x = date, y = round_id, fill = energy_dist)) +
  
  geom_line(aes(x = notification_date, y = cases), plot_data_cases |> 
              filter(notification_date >= ymd("2025-05-08"),
                     pathogen == "RSV",
                     location %in% str_c(locations, "_data"))) +
  
  scale_fill_viridis_c(trans = "sqrt") +
  
  coord_cartesian(xlim = c(ymd("2025-05-08"), ymd("2025-09-01"))) +
  
  facet_grid(location~pathogen, scales = "free_y") +
  
  theme_bw()


scores |> 
  collect() |> 
  filter(location %in% locations) |> 
  mutate(horizon = as.numeric(date - origin_date),
         round_id = as.numeric((ymd(round_id) - ymd("2025-05-08"))) / as.numeric(weeks(1))) |> 
  filter(horizon > 21, horizon <= 28) |> 
  
  ggplot() +
  
  geom_line(aes(x = date, y = energy_dist, group = round_id)) +
  
  
  geom_line(aes(x = notification_date, y = log(cases + 1)), plot_data_cases |> 
              filter(notification_date >= ymd("2025-05-08"),
                     location %in% str_c(locations, "_data"))) +
  
  coord_cartesian(xlim = c(ymd("2025-05-08"), ymd("2025-08-01"))) +
  
  facet_grid(location~pathogen, scales = "free_y") +
  
  theme_bw()

  
