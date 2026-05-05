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


perf_data |> 
  filter(horizon == 7) |> 
  ggplot() + 
  geom_tile(aes(x = round_id, y = location, fill = crps_log)) +
  
  facet_wrap(~.model)


n_days_ahead <- 7 * 4

p_skill <- perf_data |> 
  filter(horizon > 0, horizon <= n_days_ahead) |> 
  group_by(pathogen, location, round_id, date) |> 
  filter(n() == 2) |>
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
  plot_annotation("<b>A</b>. Forecasting skill",
                  theme = theme(plot.title = element_markdown(family = "Helvetica"),
                                plot.margin = margin(t = 0.2, l = 0.2, b = 0.2, r = 0.2, unit = "cm")))




perf_data |> 
  filter(horizon > 21, horizon <= n_days_ahead) |> 
  group_by(pathogen, location, round_id, date) |> 
  filter(n() == 2) |>
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
