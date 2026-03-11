
source("R/functions.R")

source("R/presentation/plot_theme.R")

perf_data <- bind_rows(
  collect(open_dataset(fs::dir_ls("results_baseline", glob = "*_performance.parquet", recurse = TRUE))) |> 
    mutate(.model = "baseline"),
  collect(open_dataset(fs::dir_ls("results_no_holidays", glob = "*_performance.parquet", recurse = TRUE))) |> 
    mutate(.model = "no_holidays"),
  collect(open_dataset(fs::dir_ls("results_no_annual", glob = "*_performance.parquet", recurse = TRUE))) |> 
    mutate(.model = "no_annual"),
  collect(open_dataset(fs::dir_ls("results_neither", glob = "*_performance.parquet", recurse = TRUE))) |> 
    mutate(.model = "neither")
) |> 
  filter(location != "NZ") |> 
  mutate(horizon = as.numeric(date - origin_date))


perf_data |> 
  filter(horizon > 21, horizon <= 28) |> 
  # group_by(.model, location, pathogen, date) |> 
  # summarise(crps_log = mean(crps_log)) |> 
  
  ggplot() +
  geom_point(aes(x = date, y = crps_log, colour = .model),
             size = 0.5) +
  
  facet_grid(rows = vars(location), cols = vars(pathogen), scales = "free_y") +
  
  theme_bw()

plot_data <- perf_data |> 
  filter(horizon > 21, horizon <= 28) |> 
  group_by(round_id, location, pathogen, horizon) |> 
  mutate(crps_log_null = crps_log[[which(.model == "neither")]]) |> 
  
  filter(.model != "neither") |> 
  group_by(pathogen, .model) |> 
  summarise(skill_score = 1 - sum(crps_log) / sum(crps_log_null)) |> 
  label_pathogens() |> 
  
  mutate(
    .model = factor(
      .model,
      c("baseline", "no_holidays", "no_annual"),
      c("With both", "With annual periodic variation", "With holiday effects")
    )
  )
  
ggplot(plot_data) +
  
  geom_vline(xintercept = 0, linetype = "44", linewidth = 0.75) +
  
  geom_linerange(aes(xmin = 0, xmax = skill_score, y = .model, colour = skill_score > 0),
                 position = position_dodge(width = 0.75),
                 linewidth = 8.0) +
  
  ggokabeito::scale_colour_okabe_ito(order = c(6, 3)) +
  ylab(NULL) + xlab("Skill score") +
  
  plot_theme_pres +
  
  facet_wrap(~pathogen, ncol = 1) +
  
  theme(legend.position = "none") +
  
  ggtitle("Forecasting skill (against null model)")



perf_data |> 
  filter(horizon > 21, horizon <= 28) |> 
  group_by(.model, pathogen) |> 
  summarise(crps_log = mean(crps_log)) |> 
  
  mutate(
    .model = factor(
      .model,
      c("baseline", "no_holidays", "no_annual", "neither"),
      c("With both", "With annual periodic term", "With holiday effects", "Null model")
    )
  ) |> 
  
  ggplot() +
  
  geom_linerange(aes(xmin = 0, xmax = crps_log, y = .model, colour = .model == "With both"),
                 linewidth = 8.0) +
  # geom_point(aes(y = .model, x = crps_log, colour = .model == "With both"),
  #            size = 8) +
  # geom_point(aes(y = .model, x = crps_log),
  #            size = 3, colour = "white") +
  
  ylab(NULL) + xlab("CRPS") +
  
  ggokabeito::scale_colour_okabe_ito(order = c(9, 2)) +
  coord_cartesian(xlim = c(0.0, 0.45)) +
  # scale_x_continuous(breaks = c(0.25, 0.3, 0.35)) +
  
  ggtitle("Mean 21–28 day CRPS\n(lower is better)") +
  
  facet_wrap(~pathogen, ncol = 1) +
  
  plot_theme_pres +
  
  theme(legend.position = "none")

