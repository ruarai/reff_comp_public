source("R/functions.R")
source("R/paper/common.R")

snaive_perf_data <- read_parquet("data/snaive_perf_extended.parquet")


perf_data <- read_parquet("saved_results/baseline_2_performance.parquet") |> 
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
  
  filter(date >= ymd("2025-01-01"),
         ymd(round_id) < ymd("2025-11-27"),
         location != "SA", location != "NZ")

combined_counts <- read_csv("data/combined-count-smoothed-2026-04-22.csv")


perf_data_paired <- perf_data |> 
  group_by(pathogen, location, round_id, date) |> 
  filter(n() == 2) |> 
  ungroup()


perf_data |> 
  filter(horizon == 7) |> 
  ggplot() + 
  geom_tile(aes(x = ymd(round_id), y = location, fill = .model),
            alpha = 0.5) +
  
  facet_grid(rows = vars(pathogen))

location_y <- tibble(
  location = unique(perf_data_paired$location)
) |> 
  mutate(y_offset = row_number())

plot_data_counts <- combined_counts |> 
  filter(notification_date >= ymd("2025-04-01"),
         notification_date <= ymd("2026-01-01")) |> 
  inner_join(location_y, by = join_by(location)) |> 
  group_by(location, pathogen) |> 
  mutate(y = y_offset + scales::rescale(value_trend, to = c(0.15, 0.9)),
         y_full = scales::rescale(value_trend)) |> 
  label_pathogen()
  


skill_days_ahead <- 7 * 3

plot_data_skill <- perf_data_paired |> 
  filter(horizon > skill_days_ahead, horizon <= skill_days_ahead + 7) |> 
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
  arrange(location, desc(date))

p_heatmap <- ggplot(plot_data_skill) +
  
  geom_rect(aes(xmin = date, xmax = date + days(7),
                ymin = y_offset - 0.5, ymax = y_offset + 0.5, fill = skill_score),
            colour = "white", linewidth = 0.25) +
  
  geom_rect(aes(xmin = date, xmax = date + days(7),
                ymin = y_offset - 0.5, ymax = y_offset - 0.4),
            fill = "black",
            colour = "white", linewidth = 0.25) +
  
  geom_linerange(aes(y = y_offset + 0.5,
                     xmin = min(plot_data_skill$date), xmax = max(plot_data_skill$date) + 14),
                 location_y |> filter(location != "VIC"),
             colour = "white", linewidth = 1.6) +
  
  geom_linerange(aes(y = y_offset - 0.5,
                     xmin = min(plot_data_skill$date), xmax = max(plot_data_skill$date) + 14),
                 location_y,
                 colour = "white", linewidth = 1.6) +
  
  geom_line(aes(x = notification_date, y = y - 0.5, group = location), plot_data_counts,
            colour = "white", linewidth = 1.5, alpha = 0.8) +
  geom_line(aes(x = notification_date, y = y - 0.5, group = location), plot_data_counts,
            linewidth = 0.4) +
  
  
  scale_fill_distiller(type = "div", direction = 1, palette = "RdBu",
                       limits = c(-1, 1), name = "Skill score",
                       breaks = c(-1, 0, 1), labels = c("≤-1", "0", "1")) +
  
  facet_wrap(~pathogen_label, ncol = 3) +
  
  scale_y_continuous(breaks = location_y$y_offset, labels = location_y$location) +
  
  scale_x_date(date_breaks = "2 months", labels = scales::label_date_short(sep = "\n")) +
  
  ylab(NULL) + xlab(NULL) +
  
  plot_theme_paper +
  
  theme(legend.position = "right",
        panel.background = element_facet_background)
p_heatmap


p_skill <- perf_data_paired |> 
  group_by(pathogen, location, round_id, date) |> 
  mutate(crps_log_snaive = crps_log[.model == "snaive"]) |> 
  filter(.model == "baseline") |> 
  group_by(pathogen, location, horizon = floor((horizon - 1) / 7)) |> 
  summarise(skill_score = 1 - mean(crps_log) / mean(crps_log_snaive)) |> 
  label_pathogen() |> 
  
  ggplot() +
  geom_line(aes(x = horizon + 1, y = skill_score, colour = location),
            linewidth = 0.9) +
  geom_point(aes(x = horizon + 1, y = skill_score, colour = location),
             size = 0.8) +
  
  geom_hline(yintercept = 0.0, linetype = "44") +
  
  ggokabeito::scale_colour_okabe_ito(order = c(1, 2, 3, 5, 6, 7, 8, 9, 4),
                                     name = "Region") +
  
  coord_cartesian(ylim = c(-0.1, NA),
                  xlim = c(0, 8)) +
  
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8)) +
  
  facet_wrap(~pathogen_label) +
  
  plot_theme_paper +
  
  xlab("Horizon (<i>n</i>th week ahead)") +
  ylab("Skill score") +
  
  theme(legend.position = "right",
        panel.background = element_facet_background,
        panel.grid.major = element_gridline)


p_top <- p_skill + 
  plot_annotation("<b>A</b>. Forecast skill by horizon",
                  theme = plot_theme_annotation)


p_bottom <- p_heatmap + 
  plot_annotation("<b>B</b>. 4th week forecast skill by epidemic period",
                  theme = plot_theme_annotation)


(wrap_elements(full = p_top) / wrap_elements(full = p_bottom)) +
  plot_layout(heights = c(1, 1.5))

ggsave(
  "paper_figures/f_2_performance.png",
  width = 10, height = 8,
  bg = "white"
)



plot_data_bias <- perf_data_paired |> 
  filter(horizon > skill_days_ahead, horizon <= skill_days_ahead + 7) |> 
  filter(.model == "baseline") |> 
  group_by(pathogen, location, round_id) |> 
  summarise(bias = mean(bias), date = first(date)) |> 
  
  ungroup() |> 
  
  left_join(location_y, by = join_by(location)) |> 
  arrange(location, desc(date)) |> 
  label_pathogen()


ggplot(plot_data_bias) +
  
  geom_rect(aes(xmin = date, xmax = date + days(7),
                ymin = y_offset - 0.5, ymax = y_offset + 0.5, fill = bias),
            colour = "white", linewidth = 0.25) +
  
  geom_rect(aes(xmin = date, xmax = date + days(7),
                ymin = y_offset - 0.5, ymax = y_offset - 0.4),
            fill = "black",
            colour = "white", linewidth = 0.25) +
  
  geom_linerange(aes(y = y_offset + 0.5,
                     xmin = min(plot_data_skill$date), xmax = max(plot_data_skill$date) + 14),
                 location_y |> filter(location != "VIC"),
                 colour = "white", linewidth = 4.0) +
  
  geom_linerange(aes(y = y_offset - 0.5,
                     xmin = min(plot_data_skill$date), xmax = max(plot_data_skill$date) + 14),
                 location_y,
                 colour = "white", linewidth = 4.0) +
  
  geom_line(aes(x = notification_date, y = y - 0.5, group = location), plot_data_counts,
            colour = "white", linewidth = 1.5, alpha = 0.8) +
  geom_line(aes(x = notification_date, y = y - 0.5, group = location), plot_data_counts,
            linewidth = 0.4) +
  
  
  scale_fill_distiller(type = "div", direction = 1, palette = "RdBu",
                       limits = c(0, 1), name = "Skill score",
                       breaks = c(0, 1), labels = c("0", "1")) +
  
  facet_wrap(~pathogen_label, ncol = 3) +
  
  scale_y_continuous(breaks = location_y$y_offset, labels = location_y$location) +
  
  ylab(NULL) +
  xlab("Date") +
  
  plot_theme_paper +
  
  theme(legend.position = "right",
        panel.background = element_facet_background)
