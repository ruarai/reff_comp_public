
source("R/functions.R")

source("R/presentation/plot_theme.R")

perf_data <- bind_rows(
  read_parquet("saved_results/results_baseline_performance.parquet") |> 
    mutate(.model = "baseline"),
  read_parquet("saved_results/normalised_incidence_0.1_performance.parquet") |> 
    mutate(.model = "normalised_incidence_0.1"),
  read_parquet("saved_results/normalised_incidence_0.2_performance.parquet") |> 
    mutate(.model = "normalised_incidence_0.2"),
  read_parquet("saved_results/normalised_incidence_0.3_performance.parquet") |> 
    mutate(.model = "normalised_incidence_0.3"),
  read_parquet("saved_results/normalised_incidence_0.4_performance.parquet") |> 
    mutate(.model = "normalised_incidence_0.4"),
  collect(open_dataset(fs::dir_ls("saved_results/normalised_incidence_0.5/", glob = "*_performance.parquet", recurse = TRUE))) |>
    mutate(.model = "normalised_incidence_0.5"),
  collect(open_dataset(fs::dir_ls("saved_results/normalised_incidence_0.75/", glob = "*_performance.parquet", recurse = TRUE))) |>
    mutate(.model = "normalised_incidence_0.75"),
  collect(open_dataset(fs::dir_ls("saved_results/normalised_incidence_1/", glob = "*_performance.parquet", recurse = TRUE))) |>
    mutate(.model = "normalised_incidence_1"),
  collect(open_dataset(fs::dir_ls("saved_results/normalised_incidence_1.5/", glob = "*_performance.parquet", recurse = TRUE))) |>
    mutate(.model = "normalised_incidence_1.5")
) |> 
  filter(location != "NZ", location != "SA") |> 
  mutate(horizon = as.numeric(date - origin_date))


perf_data |> 
  
  filter(horizon > 0, horizon <= 28) |> 
  
  group_by(pathogen, location, date, round_id) |>
  filter(n() != 9) 

perf_data |> 
  
  filter(horizon > 0, horizon <= 28) |> 
  
  group_by(pathogen, location, date, round_id) |>
  filter(n() == 9) |>
  
  group_by(pathogen, location) |> 
  mutate(crps_log_baseline = mean(crps_log[.model == "baseline"])) |> 
  filter(.model != "baseline") |> 
  group_by(pathogen, location, .model) |> 
  summarise(skill_score = 1 - mean(crps_log) / first(crps_log_baseline)) |> 
  mutate(sigma_total_inf = parse_number(.model)) |> 
  
  ggplot() +
  geom_line(aes(x = sigma_total_inf, y = skill_score, colour = location),
            linewidth = 0.5) +
  
  geom_hline(yintercept = 0, linetype = "42") +
  
  ggokabeito::scale_colour_okabe_ito(name = "Location") +
  
  xlab("sigma_log10_total_inf") +
  
  scale_y_continuous(trans = log1minus_trans()) +
  
  facet_wrap(~pathogen, ncol = 3) +
  
  theme_bw()

perf_data |> 
  
  filter(horizon > 0, horizon <= 28) |> 
  
  group_by(pathogen, location, date, round_id) |>
  filter(n() == 9) |> 
  
  group_by(pathogen) |> 
  mutate(crps_log_baseline = mean(crps_log[.model == "baseline"])) |> 
  filter(.model != "baseline") |> 
  group_by(pathogen, .model) |> 
  summarise(skill_score = 1 - mean(crps_log) / first(crps_log_baseline)) |> 
  mutate(sigma_total_inf = parse_number(.model)) |> 
  
  ggplot() +
  geom_line(aes(x = sigma_total_inf, y = skill_score),
            linewidth = 0.9) +
  
  geom_hline(yintercept = 0, linetype = "42") +
  
  ggokabeito::scale_colour_okabe_ito(name = "Location") +
  
  xlab("sigma_log10_total_inf") +
  
  scale_y_continuous(trans = log1minus_trans()) +
  
  facet_wrap(~pathogen, ncol = 3) +
  
  theme_bw()


perf_data |> 
  filter(horizon > 21, horizon <= 28) |> 
  
  ggplot() +
  geom_point(aes(x = date, y = crps_log, colour = .model),
             size = 0.5) +
  
  facet_grid(rows = vars(location), cols = vars(pathogen), scales = "free_y") +
  
  theme_bw()


perf_data |> 
  filter(horizon > 21, horizon <= 28) |> 
  group_by(pathogen, location, date, round_id) |>
  filter(n() == 9) |> 
  group_by(pathogen, location, round_id, .model) |> 
  
  summarise(crps_log = mean(crps_log),
            date = first(date)) |> 
  
  ggplot() +
  
  geom_line(aes(x = date, y = crps_log, colour = .model)) +
  
  geom_hline(yintercept = 0, linetype = "22") +
  
  facet_grid(rows = vars(location), cols = vars(pathogen), scales = "free_y") +
  
  theme_bw()




perf_data |> 
  filter(horizon > 21, horizon <= 28) |> 
  group_by(pathogen, location, date, round_id) |>
  filter(n() == 9) |> 
  group_by(pathogen, location, round_id) |> 
  mutate(crps_log_mean_baseline = mean(crps_log[.model == "baseline"])) |> 
  filter(.model != "baseline") |> 
  group_by(pathogen, location, round_id, .model) |> 
  
  summarise(skill_score = 1 - mean(crps_log) / first(crps_log_mean_baseline),
            date = first(date)) |> 
  
  ggplot() +
  
  geom_line(aes(x = date, y = skill_score, colour = .model)) +
  
  scale_colour_viridis_d(end = 0.9) +
  
  geom_hline(yintercept = 0.0, linetype = "22") +
  
  scale_y_continuous(trans = log1minus_trans()) +
  
  facet_grid(rows = vars(location), cols = vars(pathogen), scales = "free_y") +
  
  theme_bw() +
  theme(legend.position = "bottom")



expand_grid(
  i = 1:10000,
  # sigma_total_inf = c(0.1, 0.2, 0.3, 0.4)
  sigma_total_inf = c(0.5, 0.75, 1.0, 1.5)
) |> 
  group_by(sigma_total_inf) |> 
  mutate(x = 10^rnorm(10000, sd = sigma_total_inf)) |> 
  
  ggplot() +
  geom_histogram(aes(x = x), bins = 50) +
  
  facet_wrap(~sigma_total_inf) +
  
  scale_x_log10(labels = scales::label_number())


