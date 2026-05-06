source("R/functions.R")
source("R/paper/plot_theme.R")

snaive_perf_data <- read_rds("data/case_forecast_measures.rds") |> 
  filter(measure_name == "crps", target == "case_incidence") |> 
  mutate(pathogen = as.character(if_else(pathogen == "Influenza", "flu", pathogen)),
         round_id = as.character(round_id)) |> 
  filter(.model == "snaive") |> 
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



location_y <- tibble(
  location = unique(perf_data_paired$location)
) |> 
  mutate(y_offset = row_number())

combined_counts <- read_csv("data/combined-count-smoothed-2026-04-22.csv")

plot_data_counts <- combined_counts |> 
  filter(notification_date >= ymd("2025-05-01"),
         notification_date <= ymd("2025-12-01")) |> 
  inner_join(location_y, by = join_by(location)) |> 
  group_by(location, pathogen) |> 
  mutate(y = y_offset + scales::rescale(value_trend, to = c(0.1, 0.9)),
         y_full = scales::rescale(value_trend)) |> 
  label_pathogen()

perf_data_paired |> 
  filter(horizon > 21, horizon <= 28) |> 
  group_by(pathogen, location, round_id, date) |> 
  mutate(crps_log_snaive = crps_log[.model == "snaive"]) |> 
  filter(.model == "baseline") |> 
  group_by(pathogen, location, round_id) |> 
  summarise(skill_score = 1 - mean(crps_log) / mean(crps_log_snaive),
            crps_log = mean(crps_log),
            date = first(date)) |> 
  label_pathogen() |> 
  mutate(skill_score = pmax(-1, skill_score)) |> 
  
  ungroup() |> 
  
  left_join(location_y, by = join_by(location)) |> 
  arrange(location, desc(date)) |> 
  
  ggplot() +
  
  geom_rect(aes(xmin = date, xmax = date + days(7),
                ymin = y_offset - 0.5, ymax = y_offset + 0.5, fill = skill_score),
            colour = "white", linewidth = 0.25) +
  
  geom_hline(aes(yintercept = y_offset + 0.5), location_y |> filter(location != "WA"),
             colour = "white", linewidth = 4.0) +
  
  geom_line(aes(x = notification_date, y = y - 0.5, group = location), plot_data_counts,
            colour = "white", linewidth = 2.0, alpha = 0.8) +
  geom_line(aes(x = notification_date, y = y - 0.5, group = location), plot_data_counts) +
  
  
  scale_fill_distiller(type = "div", direction = 1, palette = "RdBu",
                       limits = c(-1, 1), name = "Skill score",
                       breaks = c(-1, 0, 1), labels = c("<1", "0", "1")) +
  
  facet_wrap(~pathogen_label, ncol = 3) +
  
  scale_y_continuous(breaks = location_y$y_offset, labels = location_y$location) +
  
  ylab("Region") +
  xlab("Date") +
  
  plot_theme_paper +
  
  theme(legend.position = "bottom",
        panel.background = element_facet_background)



p_skill <- perf_data_paired |> 
  filter(horizon > 21, horizon <= 28) |> 
  group_by(pathogen, location, round_id, date) |> 
  mutate(crps_log_snaive = crps_log[.model == "snaive"]) |> 
  filter(.model == "baseline") |> 
  group_by(pathogen, location, round_id) |> 
  summarise(skill_score = 1 - mean(crps_log) / mean(crps_log_snaive),
            crps_log = mean(crps_log),
            date = first(date)) |> 
  label_pathogen() |> 
  mutate(skill_score = pmax(-1, skill_score)) |> 
  
  ungroup() |> 
  
  ggplot() +
  
  geom_tile(aes(x = ymd(round_id) + days(21), y = location, fill = skill_score)) +
  
  scale_fill_distiller(type = "div", direction = 1, palette = "RdBu",
                       limits = c(-1, 1), name = "Skill score",
                       breaks = c(-1, 0, 1), labels = c("<1", "0", "1")) +
  
  facet_wrap(~pathogen_label, ncol = 3) +
  
  ylab("Region") +
  xlab("Date") +
  
  plot_theme_paper +
  
  theme(legend.position = "none",
        panel.background = element_facet_background)

p_inc <- ggplot(plot_data_counts) +
  geom_line(aes(x = notification_date, y = y_full, colour = location)) +
  
  facet_wrap(~pathogen_label, ncol = 3) +
  
  plot_theme_paper +
  
  ylab("Relative case incidence") +
  xlab("Date") +
  
  theme(legend.position = "bottom",
        panel.background = element_facet_background)


wrap_plots(p_skill, p_inc, heights = c(2, 1))
