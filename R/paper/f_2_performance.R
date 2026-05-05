source("R/functions.R")
source("R/paper/plot_theme.R")

snaive_perf_data <- read_rds("data/case_forecast_measures.rds") |> 
  filter(measure_name == "crps", target == "case_incidence") |> 
  mutate(pathogen = as.character(if_else(pathogen == "Influenza", "flu", pathogen)),
         round_id = as.character(round_id)) |> 
  filter(.model == "ensemble-mixture") |> 
  select(-c(measure_name, target)) |> 
  pivot_wider(names_from = measure_scale) |> 
  rename(crps = natural, crps_log = log)



perf_data <- read_parquet("saved_results/results_baseline_performance.parquet") |> 
  mutate(.model = "baseline") |> 
  bind_rows(snaive_perf_data) |> 
  mutate(horizon = as.numeric(date - origin_date)) |> 
  
  # Drop dates with known issues
  anti_join(
    tribble(
      ~location, ~pathogen, ~date,
      "WA", "flu", ymd("2025-08-16"),
      "WA", "flu", ymd("2025-08-23"),
      "WA", "flu", ymd("2025-08-30"),
      "WA", "RSV", ymd("2025-08-16"),
      "WA", "RSV", ymd("2025-08-23"),
      "WA", "RSV", ymd("2025-08-30")
    )
  ) |> 
  filter(location != "NZ", location != "SA")

perf_data_paired <- perf_data |> 
  group_by(pathogen, location, round_id, date) |> 
  filter(n() == 2) |> 
  ungroup()


perf_data |> 
  filter(horizon == 7) |> 
  ggplot() + 
  geom_tile(aes(x = ymd(round_id), y = location, fill = crps_log)) +
  
  facet_grid(cols = vars(.model), rows = vars(pathogen))

n_days_ahead <- 7 * 4

p_skill <- perf_data_paired |> 
  filter(horizon > 0, horizon <= n_days_ahead) |> 
  group_by(pathogen, location, round_id, date) |> 
  mutate(crps_log_ensemble = crps_log[.model == "ensemble-mixture"]) |> 
  filter(.model == "baseline") |> 
  group_by(pathogen, location) |> 
  summarise(skill_score = 1 - mean(crps_log) / mean(crps_log_ensemble),
            crps_log = mean(crps_log)) |> 
  label_pathogen() |> 
  
  ggplot() +
  
  geom_linerange(aes(x = location, ymin = 0, ymax = skill_score),
                 linewidth = 6.0) +
  
  geom_hline(yintercept = 0) +
  
  # ggokabeito::scale_colour_okabe_ito(order = c(6, 3)) +
  
  
  scale_y_continuous(breaks = c(-0.3, -0.15, 0, 0.15, 0.3)) +
  
  facet_wrap(~pathogen_label, ncol = 3) +
  
  xlab("Region") + ylab("Skill score") +
  
  plot_theme_paper +
  
  theme(legend.position = "none",
        panel.grid.major.y = element_gridline,
        panel.grid.major.x = element_gridline,
        panel.background = element_facet_background)


p_skill + 
  plot_annotation("<b>A</b>. Forecasting skill (against ACEFA 2025 ensemble forecast)",
                  theme = theme(plot.title = element_markdown(family = "Helvetica"),
                                plot.margin = margin(t = 0.2, l = 0.2, b = 0.2, r = 0.2, unit = "cm")))

perf_data_paired |> 
  filter(horizon > 21, horizon <= n_days_ahead) |> 
  group_by(pathogen, location, round_id, date) |> 
  mutate(crps_log_ensemble = crps_log[.model == "ensemble-mixture"]) |> 
  filter(.model == "baseline") |> 
  group_by(pathogen, location, round_id) |> 
  summarise(skill_score = 1 - mean(crps_log) / mean(crps_log_ensemble),
            crps_log = mean(crps_log)) |> 
  label_pathogen() |> 
  mutate(round_id = ymd(round_id),
         skill_score = pmax(-1, skill_score)) |> 
  
  ungroup() |> 
  
  complete(location, pathogen_label, round_id) |> 
  
  ggplot() +
  
  geom_tile(aes(x = round_id, y = location, fill = skill_score)) +
  
  scale_fill_distiller(type = "div", direction = 1, palette = "RdBu",
                       limits = c(-1, 1), name = "Skill score",
                       breaks = c(-1, 0, 1), labels = c("<1", "0", "1")) +
  
  facet_wrap(~pathogen_label, ncol = 1) +
  
  xlab("Region") + ylab("Skill score") +
  
  plot_theme_paper +
  
  theme(legend.position = "bottom",
        panel.background = element_facet_background)


perf_data_paired |> 
  filter(horizon > 0, horizon <= n_days_ahead) |> 
  group_by(pathogen, location, round_id, date) |> 
  mutate(crps_log_ensemble = crps_log[.model == "ensemble-mixture"]) |> 
  group_by(pathogen, location, horizon = floor((horizon - 1) / 7)) |> 
  summarise(skill_score = 1 - mean(crps_log) / mean(crps_log_ensemble),
            crps_log = mean(crps_log)) |> 
  label_pathogen() |> 
  ggplot() +
  geom_line(aes(x = horizon, y = skill_score, colour = location)) +
  
  ggokabeito::scale_colour_okabe_ito(order = c(1, 2, 3, 5, 6, 7, 8),
                                     name = "Region") +
  
  geom_hline(yintercept = 0, linetype = "44") +
  
  plot_theme_paper +
  
  scale_x_continuous(breaks = c(0, 1, 2, 3),
                     labels = c("1–7", "8–14", "15–21", "22–28")) +
  
  xlab("Horizon (days)") + ylab("Skill score") +
  
  facet_wrap(~pathogen_label) +
  
  theme(panel.spacing.x = unit(0.5, "cm"))

perf_data |> 
  filter(horizon > 0, horizon <= n_days_ahead) |> 
  group_by(pathogen, location, horizon = floor((horizon - 1) / 7)) |> 
  summarise(crps_log = mean(crps_log)) |> 
  label_pathogen() |> 
  ggplot() +
  geom_line(aes(x = horizon, y = crps_log, colour = location)) +
  geom_point(aes(x = horizon, y = crps_log, colour = location), size = 1) +
  
  ggokabeito::scale_colour_okabe_ito(order = c(1, 2, 3, 5, 6, 7, 8),
                                     name = "Region") +
  
  plot_theme_paper +
  
  scale_x_continuous(breaks = c(0, 1, 2, 3),
                     labels = c("1–7", "8–14", "15–21", "22–28")) +
  
  coord_cartesian(ylim = c(0, NA)) +
  
  xlab("Horizon (days)") + ylab("CRPS") +
  
  facet_wrap(~pathogen_label) +
  
  theme(panel.spacing.x = unit(0.5, "cm"))


combined_counts <- read_csv("data/combined-count-smoothed-2026-04-22.csv")

growth_estimates <- combined_counts |> 
  filter(target == "case_incidence") |> 
  group_by(pathogen, location) |>
  mutate(g = log(value_trend + 1) - lag(log(value_trend + 1), 21))

perf_data |> 
  filter(horizon > 0, horizon <= 28) |> 
  left_join(growth_estimates |> select(date = notification_date, location, pathogen, g)) |> 
  
  mutate(g = cut(g, c(-10, -0.1, 0.1, 10))) |> 
  
  group_by(g, location, pathogen) |> 
  summarise(crps_log = mean(crps_log)) |> 
  
  label_pathogen() |> 
  
  
  ggplot() +
  geom_point(aes(x = g, y = crps_log, colour = location), size = 0.5, stroke = 0.5) +
  
  facet_wrap(~pathogen_label, ncol = 3) +
  plot_theme_paper



