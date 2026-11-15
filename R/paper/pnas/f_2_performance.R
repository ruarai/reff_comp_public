source("R/functions.R")
source("R/paper/common.R")

perf_summary <- read_perf_summary(c("baseline", "snaive"))

combined_counts <- read_csv("data/combined-count-smoothed_final_anonymised.csv")




perf_summary |> 
  filter(horizon_week == 0) |> 
  ggplot() + 
  geom_tile(aes(x = ymd(round_id), y = location, fill = .model),
            alpha = 0.5) +
  
  facet_grid(rows = vars(pathogen))




plot_data_counts <- combined_counts |> 
  filter(notification_date >= ymd("2025-04-01") | (pathogen == "RSV" & notification_date >= ymd("2025-02-01")),
         notification_date <= ymd("2026-01-01")) |> 
  inner_join(location_y, by = join_by(location)) |> 
  group_by(location, pathogen) |> 
  mutate(y_clipped = pmin(value_trend_and_remainder, max(value_trend) * 1.2),
         y = y_offset + scales::rescale(y_clipped, to = c(0.15, 0.9), from = c(0, max(value_trend) * 1.2)),
         y_trend = y_offset + scales::rescale(value_trend, to = c(0.15, 0.9), from = c(0, max(value_trend) * 1.2)),
         y_full = scales::rescale(value_trend)) |> 
  label_pathogen()

example_points <- tribble(
  ~round_id, ~location, ~pathogen, ~label,
  "2025-04-04", "QLD", "RSV", "i",
  # "2025-06-05", "VIC", "flu", "ii",
  # "2025-04-24", "ACT", "SARSCOV2", "iii",
  "2025-05-08", "VIC", "SARSCOV2", "ii"
) |> 
  label_pathogen() |> 
  mutate(notification_date = ymd(round_id) + days(23),
         approx_date = ymd(round_id) + 23.5) |> 
  
  left_join(location_y, by = join_by(location)) |> 
  left_join(plot_data_counts)


plot_data_skill <- perf_summary |> 
  filter(.model == "baseline", horizon_week == 3) |> 
  select(pathogen, location, round_id, crps_log, date = first_date) |> 
  label_pathogen() |> 
  mutate(crps_log = pmin(0.75, crps_log)) |> 
  
  left_join(location_y, by = join_by(location)) |> 
  arrange(location, desc(date))


school_holidays <- get_school_holidays() |> 
  filter(location %in% locations_paper) |> 
  left_join(location_y, by = join_by(location)) |> 
  filter(year(start_date) == 2025) |> 
  expand_grid(pathogen = c("RSV", "flu", "SARSCOV2")) |> 
  label_pathogen()


# Tight titles, strips and margins for the compact PNAS layout
title_theme <- theme(
  plot.title.position = "plot",
  plot.title = element_markdown(family = "Helvetica", margin = margin(b = 2)),
  strip.text.x = element_text(size = 11, hjust = 0, margin = margin(t = 0, b = 2)),
  plot.margin = margin(0.1, 0.05, 0.1, 0.05, unit = "cm")
)


p_heatmap <- ggplot(plot_data_skill) +
  
  geom_rect(aes(xmin = date, xmax = date + days(7),
                ymin = y_offset - 0.4, ymax = y_offset + 0.4, fill = crps_log,
                label = round_id),
            colour = "white", linewidth = 0.25) +
  
  geom_rect(aes(xmin = date, xmax = date + days(7),
                ymin = y_offset - 0.44, ymax = y_offset - 0.4),
            fill = "black",
            colour = "white", linewidth = 0.25) +
  
  geom_rect(aes(xmin = start_date, xmax = end_date,
                ymin = y_offset - 0.44, ymax = y_offset - 0.48),
            fill = "#017A4AFF", school_holidays,
            linewidth = 0.25) +
  
  geom_line(aes(x = notification_date, y = y_trend - 0.5, group = location), plot_data_counts,
            colour = "white", linewidth = 1.5, alpha = 0.8) +
  geom_line(aes(x = notification_date, y = y_trend - 0.5, group = location), plot_data_counts,
            linewidth = 0.4) +
  
  geom_point(aes(x = approx_date, y = y_trend - 0.5),
             example_points, size = 2) +
  geom_point(aes(x = approx_date, y = y_trend - 0.5),
             colour = "white", size = 0.5, stroke = 0.5,
             example_points) +
  
  geom_label(aes(x = approx_date, y = y_offset - 0.3, label = label),
             size = 3, alpha = 0.8,
             example_points) +
  
  
  scale_fill_gradientn(
    name = "CRPS", limits = c(0, 0.75),
    breaks = c(0, 0.25, 0.5, 0.75), labels = c("0", "0.25", "0.5", "≥0.75"),
    colours = RColorBrewer::brewer.pal(10, "RdYlBu") |> rev()
  ) +
  
  facet_wrap(~pathogen_label, ncol = 3,
             scales = "free_x") +
  
  scale_y_continuous(breaks = location_y$y_offset, labels = location_y$location) +
  
  scale_x_date(date_breaks = "2 months", labels = scales::label_date_short(sep = "\n"),
               date_minor_breaks = "1 months",
               guide = guide_axis(minor.ticks = TRUE)) +
  
  ylab("Jurisdiction") + xlab(NULL) +
  
  plot_theme_paper +
  
  coord_cartesian(xlim = c(NA_Date_, ymd("2026-01-01")), ylim = c(0.5, 5.5)) +
  
  theme(legend.position = "right",
        panel.background = element_facet_background)
p_heatmap



p_crps <- perf_summary |> 
  filter(.model == "baseline") |> 
  group_by(pathogen, location, horizon = horizon_week) |> 
  summarise(crps_log = mean(crps_log)) |> 
  label_pathogen() |> 
  
  ggplot() +
  geom_line(aes(x = horizon + 1, y = crps_log, colour = location),
            linewidth = 0.9) +
  geom_point(aes(x = horizon + 1, y = crps_log, colour = location),
             size = 0.8) +
  
  ggokabeito::scale_colour_okabe_ito(order = c(1, 2, 3, 5, 6, 7, 8, 9, 4),
                                     name = "Jurisdiction") +
  
  coord_cartesian(ylim = c(-0, NA), xlim = c(1, 4)) +
  scale_x_continuous(breaks = c(1, 2, 3, 4)) +
  
  facet_wrap(~pathogen_label) +
  
  plot_theme_paper +
  
  xlab("Horizon (<i>n</i>th week ahead)") +
  ylab("CRPS") +
  
  theme(legend.position = "right",
        panel.background = element_facet_background,
        panel.grid.major = element_gridline)


p_skill <- perf_summary |> 
  group_by(pathogen, location, round_id, horizon_week) |> 
  mutate(crps_log_snaive = crps_log[.model == "snaive"]) |> 
  filter(.model == "baseline") |> 
  group_by(pathogen, location, horizon = horizon_week) |> 
  summarise(skill_score = 1 - mean(crps_log) / mean(crps_log_snaive)) |> 
  label_pathogen() |> 
  
  ggplot() +
  geom_line(aes(x = horizon + 1, y = skill_score, colour = location),
            linewidth = 0.9) +
  geom_point(aes(x = horizon + 1, y = skill_score, colour = location),
             size = 0.8) +
  
  geom_hline(yintercept = 0.0, linetype = "44") +
  
  ggokabeito::scale_colour_okabe_ito(order = c(1, 2, 3, 5, 6, 7, 8, 9, 4),
                                     name = "Jurisdiction") +
  
  coord_cartesian(ylim = c(-0.1, 0.7), xlim = c(1, 4)) +
  scale_x_continuous(breaks = c(1, 2, 3, 4)) +
  
  facet_wrap(~pathogen_label) +
  
  plot_theme_paper +
  
  xlab("Horizon (<i>n</i>th week ahead)") +
  ylab("Skill score") +
  
  theme(legend.position = "right",
        panel.background = element_facet_background,
        panel.grid.major = element_gridline)


p_skill <- perf_summary |> 
  group_by(pathogen, location, round_id, horizon_week) |> 
  mutate(crps_log_snaive = crps_log[.model == "snaive"]) |> 
  filter(.model == "baseline") |> 
  group_by(pathogen, location) |> 
  summarise(skill_score = 1 - mean(crps_log) / mean(crps_log_snaive)) |> 
  label_pathogen() |> 
  
  ggplot() +
  geom_linerange(aes(x = location, ymin = 0, ymax = skill_score, colour = location),
                 linewidth = 5) +
  
  geom_hline(yintercept = 0.0, linetype = "44") +
  
  ggokabeito::scale_colour_okabe_ito(order = c(1, 2, 3, 5, 6, 7, 8, 9, 4),
                                     name = "Jurisdiction") +
  
  coord_cartesian(ylim = c(-0.3, 0.7)) +
  
  facet_wrap(~pathogen_label) +
  
  plot_theme_paper +
  
  xlab("Jurisdiction") +
  ylab("Skill score") +
  
  theme(legend.position = "none",
        panel.background = element_facet_background,
        panel.grid.major = element_gridline)


p_top <- p_crps +
  ggtitle("<b>A</b>. Forecast performance by horizon week (lower value is better)") +
  title_theme

p_mid <- p_heatmap +
  ggtitle("<b>B</b>. 4th week forecast performance by epidemic period") +
  title_theme

p_bottom <- p_skill +
  ggtitle("<b>C</b>. Forecast skill against baseline by jurisdiction (higher value is better)") +
  xlab(NULL) +
  title_theme

(p_top / p_mid / p_bottom) +
  plot_layout(heights = c(1, 1.5, 1))

ggsave(
  "pnas_figures/f_2_performance.png",
  width = 0.81 * 9, height = 6.3,
  bg = "white"
)

