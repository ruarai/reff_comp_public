
source("R/functions.R")

source("R/presentation/plot_theme.R")

perf_files <- fs::dir_ls("saved_results/results_baseline/", glob = "*_performance.parquet", recurse = TRUE)
perf_data <- collect(open_dataset(perf_files)) |> 
  mutate(.model = "uom-reff-comp", .before = 1)


acefa_perf_data <- read_rds("data/case_forecast_measures.rds") |> 
  filter(measure_name == "crps", target == "case_incidence") |> 
  mutate(pathogen = as.character(if_else(pathogen == "Influenza", "flu", pathogen)),
         round_id = as.character(round_id)) |> 
  select(-c(measure_name, target)) |> 
  pivot_wider(names_from = measure_scale) |> 
  rename(crps = natural, crps_log = log)



perf_data_full <- bind_rows(
  perf_data, acefa_perf_data
) |> 
  select(.model, location, pathogen, round_id, origin_date, date, crps, crps_log) |> 
  mutate(horizon = as.numeric(date - origin_date)) |> 
  filter(horizon > 0, horizon <= 28, location != "NZ") |> 
  
  mutate(
    .model = factor(
      .model, 
      c("snaive", "dstg-RV_model", "UQ-DD", "UQ-renewal_global", 
        "UQ-renewal_i", "uom-multiwalk", "ensemble-mixture", "ensemble-quantile", 
        "uom-reff-comp")
    )
  )


perf_data_full |> 
  filter(horizon > 0, !str_detect(.model, "ensemble")) |> 
  group_by(pathogen, location, .model, horizon) |> 
  summarise(crps_log = mean(crps_log)) |> 
  
  ggplot() +
  geom_line(aes(x = horizon, y = crps_log, colour = .model == "uom-reff-comp", group = .model)) +
  
  ggokabeito::scale_colour_okabe_ito(order = c(2, 9)) +
  
  facet_grid(rows = vars(location), cols = vars(pathogen), scales = "free_y") +
  
  theme_bw()

perf_data_ranking <- perf_data_full |> 
  filter(!str_detect(.model, "ensemble"), #.model != "snaive",
         horizon > 21, horizon <= 28) |> 
  group_by(pathogen, location, round_id, horizon) |> 
  mutate(n_models = n(),
         rank = min_rank(crps_log)) |> 
  ungroup() |> 
  group_by(pathogen, location) |> 
  filter(n_models == max(n_models)) |> 
  group_by(.model) |> 
  count(rank) |> 
  mutate(n_cumulative = cumsum(n),
         p = n_cumulative / sum(n))




perf_data_ranking |>
  ggplot() +
  geom_linerange(aes(x = rank, ymin = 0, ymax = p, colour = .model),
                 position = position_dodge(width = 0.5),
                 linewidth = 2.0) +
  
  ggokabeito::scale_colour_okabe_ito(
    order = c(9, 9, 9, 9, 9, 9, 3),
    name = "Model"
  ) +
  
  coord_cartesian(xlim = c(0.5, 6)) +
  scale_x_continuous(breaks = 1:6) +
  
  xlab("Rank no worse than") +
  ylab("Proportion") +
  
  ggtitle("Performance ranking") +
  
  plot_theme_pres


perf_data_full |> 
  filter(horizon > 21, horizon <= 28) |> 
  group_by(location, pathogen, round_id, .model) |> 
  summarise(crps_log = mean(crps_log),
            date = first(date)) |> 
  
  ggplot() +
  geom_line(aes(x = date, y = crps_log, colour = .model == "uom-reff-comp", group = .model)) +
  
  ggokabeito::scale_colour_okabe_ito(order = c(2, 9)) +
  
  facet_grid(rows = vars(location), cols = vars(pathogen), scales = "free_y") +
  
  theme_bw()

perf_data_full |> 
  filter(horizon > 21, horizon <= 28,
         !str_detect(.model, "ensemble")) |> 
  group_by(.model) |> 
  summarise(crps_log = mean(crps_log)) |> 
  
  ggplot() +
  geom_linerange(aes(y = .model, xmin = 0, xmax = crps_log, colour = .model == "uom-reff-comp"),
                 linewidth = 14.0) +
  
  ylab(NULL) + xlab("CRPS") +
  
  ggokabeito::scale_colour_okabe_ito(order = c(9, 2)) +
  
  ggtitle("Mean 22–28 day CRPS\n(lower is better)") +
  
  plot_theme_pres +
  
  theme(legend.position = "none")

perf_data_full |> 
  filter(horizon > 21, horizon <= 28,
         !str_detect(.model, "ensemble")) |> 
  group_by(.model, location) |> 
  summarise(crps_log = mean(crps_log)) |> 
  
  ggplot() +
  geom_linerange(aes(y = .model, xmin = 0, xmax = crps_log, colour = .model == "uom-reff-comp"),
                 linewidth = 1) +
  
  ylab(NULL) + xlab("CRPS") +
  
  ggokabeito::scale_colour_okabe_ito(order = c(9, 2)) +
  
  facet_grid(rows = vars(location)) +
  
  ggtitle("Mean 22–28 day CRPS\n(lower is better)") +
  
  plot_theme_pres +
  
  theme(legend.position = "none")

perf_data_full |> 
  filter(horizon > 21, horizon <= 28,
         !str_detect(.model, "ensemble"), .model != "snaive") |> 
  group_by(.model) |> 
  summarise(crps_log = mean(crps_log)) |> 
  
  ggplot() +
  
  geom_linerange(aes(xmin = 0, xmax = crps_log, y = .model, colour = .model == "uom-reff-comp"),
                 linewidth = 2.0) +
  geom_point(aes(y = .model, x = crps_log, colour = .model == "uom-reff-comp"),
             size = 8) +
  geom_point(aes(y = .model, x = crps_log),
             size = 3, colour = "white") +
  
  ylab(NULL) + xlab("CRPS") +
  
  ggokabeito::scale_colour_okabe_ito(order = c(9, 2)) +
  
  ggtitle("Mean 22–28 day CRPS\n(lower is better)") +
  
  plot_theme_pres +
  
  theme(legend.position = "none")

perf_data_full |> 
  filter(horizon > 21, horizon <= 28,
         !str_detect(.model, "ensemble"), .model != "snaive") |> 
  
  ggplot() +
  ggdist::stat_histinterval(aes(y = .model, x = crps_log, colour = .model == "uom-reff-comp"),
                            breaks = seq(0, 10, by = 0.01)) +
  
  ylab(NULL) + xlab("CRPS") +
  
  ggokabeito::scale_colour_okabe_ito(order = c(9, 2)) +
  
  coord_cartesian(xlim = c(0, 1)) +
  
  ggtitle("21–28 day CRPS (lower is better)") +
  
  plot_theme_pres +
  
  theme(legend.position = "none")

  
perf_data_full |> 
  filter(horizon > 21, horizon <= 28, .model %in% c("snaive", "uom-reff-comp")) |> 
  group_by(location, pathogen, round_id, .model) |> 
  summarise(crps_log = mean(crps_log),
            date = first(date)) |> 
  
  ggplot() +
  geom_line(aes(x = date, y = crps_log, colour = .model, group = .model)) +
  
  ggokabeito::scale_colour_okabe_ito(order = c(2, 9)) +
  
  facet_grid(rows = vars(location), cols = vars(pathogen), scales = "free_y") +
  
  theme_bw()


perf_data_full |> 
  filter(horizon > 21, horizon <= 28) |>
  
  group_by(location, pathogen, round_id, horizon) |> 
  mutate(crps_log_snaive = crps_log[.model == "snaive"]) |> 
  filter(.model != "snaive", !str_detect(.model, "ensemble")) |> 
  
  group_by(location, pathogen, .model) |> 
  
  summarise(skill_score = 1 - sum(crps_log) / sum(crps_log_snaive)) |> 
  
  ggplot() +
  geom_linerange(aes(xmin = 0, xmax = skill_score, y = .model, colour = skill_score > 0),
                 linewidth = 2.0) +
  
  facet_grid(rows = vars(location), cols = vars(pathogen)) +
  
  ggokabeito::scale_colour_okabe_ito(order = c(6, 3)) +
  ylab(NULL) + xlab("Skill score") +
  
  theme_bw() +
  
  theme(legend.position = "none")


perf_data_full |> 
  filter(horizon > 21, horizon <= 28) |>
  
  group_by(location, pathogen, round_id, horizon) |> 
  mutate(crps_log_snaive = crps_log[.model == "snaive"]) |> 
  filter(.model != "snaive") |> 
  
  group_by(.model) |> 
  
  summarise(skill_score = 1 - sum(crps_log) / sum(crps_log_snaive)) |> 
  
  ggplot() +
  geom_linerange(aes(xmin = 0, xmax = skill_score, y = .model, colour = skill_score > 0),
                 linewidth = 2.0) +
  
  ggokabeito::scale_colour_okabe_ito(order = c(6, 3)) +
  ylab("Model") + xlab("Skill score") +
  
  plot_theme_pres +
  
  theme(legend.position = "none") +
  
  ggtitle("Forecasting skill (against naive model)")

perf_data_full |> 
  filter(horizon > 21, horizon <= 28) |>
  
  group_by(location, pathogen, round_id, horizon) |> 
  mutate(crps_log_snaive = crps_log[.model == "snaive"]) |> 
  filter(.model != "snaive") |> 
  
  group_by(.model) |> 
  
  summarise(skill_score = 1 - sum(crps_log) / sum(crps_log_snaive)) |> 
  
  ggplot() +
  geom_linerange(aes(xmin = 0, xmax = skill_score, y = .model, colour = skill_score > 0),
                 linewidth = 2.0) +
  
  ggokabeito::scale_colour_okabe_ito(order = c(6, 3)) +
  ylab("Model") + xlab("Skill score") +
  
  plot_theme_pres +
  
  theme(legend.position = "none") +
  
  ggtitle("Forecasting skill (against naive model)")


perf_data_full |> 
  filter(horizon > 21, horizon <= 28) |>
  
  group_by(location, pathogen, round_id, horizon) |> 
  mutate(crps_log_snaive = crps_log[.model == "snaive"]) |> 
  filter(.model == "uom-reff-comp") |> 
  
  group_by(pathogen) |> 
  
  summarise(skill_score = 1 - sum(crps_log) / sum(crps_log_snaive)) |> 
  
  label_pathogens() |> 
  
  ggplot() +
  geom_vline(xintercept = 0, linetype = "24", linewidth = 1.0) +
  geom_linerange(aes(xmin = 0, xmax = skill_score, y = pathogen),
                 linewidth = 12.0) +
  ylab(NULL) + xlab("Skill score") +
  
  coord_cartesian(xlim = c(-0.2, 0.6)) +
  
  plot_theme_pres +
  
  theme(legend.position = "none") +
  
  ggtitle("Forecasting skill (higher is better)")

