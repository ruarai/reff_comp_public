source("R/functions.R")
source("R/paper/common.R")

perf_summary <- read_perf_summary(c("baseline", "no_holidays_no_annual", "snaive"))

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
  mutate(y = y_offset + scales::rescale(value_trend, to = c(0.15, 0.9)),
         y_full = scales::rescale(value_trend)) |> 
  label_pathogen()



school_holidays <- get_school_holidays() |> 
  filter(location %in% locations_paper) |> 
  left_join(location_y, by = join_by(location)) |> 
  filter(year(start_date) == 2025) |> 
  expand_grid(pathogen = c("RSV", "flu", "SARSCOV2")) |> 
  label_pathogen()


plot_data_skill_snaive <- perf_summary |> 
  filter(horizon_week == 3) |> 
  group_by(pathogen, location, round_id) |> 
  mutate(crps_log_snaive = crps_log[.model == "snaive"]) |> 
  filter(.model == "baseline") |> 
  mutate(skill_score = pmax(1 - crps_log / crps_log_snaive, -1)) |> 
  select(pathogen, location, round_id, skill_score, date = first_date) |> 
  label_pathogen() |> 
  left_join(location_y, by = join_by(location)) |> 
  arrange(location, desc(date))

p_heatmap_snaive <- ggplot(plot_data_skill_snaive) +
  
  geom_rect(aes(xmin = date, xmax = date + days(7),
                ymin = y_offset - 0.4, ymax = y_offset + 0.4, fill = skill_score,
                label = round_id),
            colour = "white", linewidth = 0.25) +
  
  geom_rect(aes(xmin = date, xmax = date + days(7),
                ymin = y_offset - 0.46, ymax = y_offset - 0.4),
            fill = "black",
            colour = "white", linewidth = 0.25) +
  
  geom_rect(aes(xmin = start_date, xmax = end_date,
                ymin = y_offset - 0.44, ymax = y_offset - 0.48),
            fill = "#017A4AFF", school_holidays,
            linewidth = 0.25) +
  
  geom_line(aes(x = notification_date, y = y - 0.5, group = location), plot_data_counts,
            colour = "white", linewidth = 1.5, alpha = 0.8) +
  geom_line(aes(x = notification_date, y = y - 0.5, group = location), plot_data_counts,
            linewidth = 0.4) +
  
  scale_fill_gradientn(
    name = "Skill score", limits = c(-1, 1),
    breaks = c(-1, -0.5, 0, 0.5, 1), labels = c("<-1", "-0.5", "0.0", "0.5", "1.0"),
    colours = RColorBrewer::brewer.pal(10, "RdBu")
  ) +
  
  facet_wrap(~pathogen_label, ncol = 3,
             scales = "free_x") +
  
  scale_y_continuous(breaks = location_y$y_offset, labels = location_y$location) +
  
  scale_x_date(date_breaks = "2 months", labels = scales::label_date_short(sep = "\n"),
               date_minor_breaks = "1 months",
               guide = guide_axis(minor.ticks = TRUE)) +
  
  ylab("Jurisdiction") + xlab(NULL) +
  
  plot_theme_paper +
  
  ggtitle("<b>A.</b> Fourth week forecast skill against random-walk baseline model") +
  
  theme(legend.position = "right",
        panel.background = element_facet_background)

plot_data_skill_null <- perf_summary |> 
  filter(horizon_week == 3) |> 
  group_by(pathogen, location, round_id) |> 
  mutate(crps_log_snaive = crps_log[.model == "no_holidays_no_annual"]) |> 
  filter(.model == "baseline") |> 
  mutate(skill_score = pmax(1 - crps_log / crps_log_snaive, -1)) |> 
  select(pathogen, location, round_id, skill_score, date = first_date) |> 
  label_pathogen() |> 
  left_join(location_y, by = join_by(location)) |> 
  arrange(location, desc(date))


p_heatmap_null <- ggplot(plot_data_skill_null) +
  
  geom_rect(aes(xmin = date, xmax = date + days(7),
                ymin = y_offset - 0.4, ymax = y_offset + 0.4, fill = skill_score,
                label = round_id),
            colour = "white", linewidth = 0.25) +
  
  geom_rect(aes(xmin = date, xmax = date + days(7),
                ymin = y_offset - 0.46, ymax = y_offset - 0.4),
            fill = "black",
            colour = "white", linewidth = 0.25) +
  
  geom_rect(aes(xmin = start_date, xmax = end_date,
                ymin = y_offset - 0.44, ymax = y_offset - 0.48),
            fill = "#017A4AFF", school_holidays,
            linewidth = 0.25) +
  
  geom_line(aes(x = notification_date, y = y - 0.5, group = location), plot_data_counts,
            colour = "white", linewidth = 1.5, alpha = 0.8) +
  geom_line(aes(x = notification_date, y = y - 0.5, group = location), plot_data_counts,
            linewidth = 0.4) +
  
  scale_fill_gradientn(
    name = "Skill score", limits = c(-1, 1),
    breaks = c(-1, -0.5, 0, 0.5, 1), labels = c("<-1", "-0.5", "0.0", "0.5", "1.0"),
    colours = RColorBrewer::brewer.pal(10, "RdBu")
  ) +
  
  facet_wrap(~pathogen_label, ncol = 3,
             scales = "free_x") +
  
  scale_y_continuous(breaks = location_y$y_offset, labels = location_y$location) +
  
  scale_x_date(date_breaks = "2 months", labels = scales::label_date_short(sep = "\n"),
               date_minor_breaks = "1 months",
               guide = guide_axis(minor.ticks = TRUE)) +
  
  ylab("Jurisdiction") + xlab(NULL) +
  
  plot_theme_paper +
  
  ggtitle("<b>B.</b> Fourth week forecast skill against model without holiday or seasonal effects") +
  
  theme(legend.position = "right", panel.background = element_facet_background)


p_heatmap_snaive / p_heatmap_null


ggsave(
  "paper_figures/supp/f_skill_heatmap.png",
  width = 10, height = 9,
  bg = "white"
)

