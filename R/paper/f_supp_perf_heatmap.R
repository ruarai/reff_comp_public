source("R/functions.R")
source("R/paper/common.R")

perf_summary <- read_perf_summary("baseline")

combined_counts <- read_csv("data/combined-count-smoothed_final_anonymised.csv")




plot_data_counts <- combined_counts |> 
  filter(notification_date >= ymd("2025-04-01") | (pathogen == "RSV" & notification_date >= ymd("2025-02-01")),
         notification_date <= ymd("2026-01-01")) |> 
  inner_join(location_y, by = join_by(location)) |> 
  group_by(location, pathogen) |> 
  mutate(y = scales::rescale(value_trend, to = c(0.15 * 4, 0.9 * 4)),
         y_full = scales::rescale(value_trend)) |> 
  label_pathogen() |> 
mutate(location = factor(location, c("NSW", "QLD", "ACT", "NT", "VIC") |> rev()))


plot_data_skill <- perf_summary |> 
  mutate(days_ahead_group = horizon_week * 7) |> 
  select(days_ahead_group, pathogen, location, round_id, crps_log, date = first_date) |> 
  label_pathogen() |> 
  mutate(crps_log = pmin(0.75, crps_log)) |> 
  
  left_join(location_y, by = join_by(location)) |> 
  arrange(location, desc(date)) |> 
  mutate(location = factor(location, c("NSW", "QLD", "ACT", "NT", "VIC") |> rev())) |> 
  
  mutate(days_ahead_group_label = str_c(days_ahead_group + 1, "-", days_ahead_group + 7, " days") |> 
           factor(c("1-7 days", "8-14 days", "15-21 days", "22-28 days")))

p_heatmap <- ggplot(plot_data_skill) +
  
  geom_rect(aes(xmin = date, xmax = date + days(7),
                ymin = days_ahead_group / 7 - 0.4, ymax = days_ahead_group / 7 + 0.4, fill = crps_log,
                label = round_id),
            colour = "white", linewidth = 0.25) +
  
  geom_rect(aes(xmin = date, xmax = date + days(7),
                ymin = days_ahead_group / 7 - 0.46, ymax = days_ahead_group / 7 - 0.4),
            fill = "black",
            colour = "white", linewidth = 0.25) +
  
  geom_line(aes(x = notification_date, y = y - 0.5, group = location), plot_data_counts,
            colour = "white", linewidth = 1.5, alpha = 0.8) +
  geom_line(aes(x = notification_date, y = y - 0.5, group = location), plot_data_counts,
            linewidth = 0.4) +
  
  
  scale_fill_gradientn(
    name = "CRPS", limits = c(0, 0.75),
    breaks = c(0, 0.25, 0.5, 0.75), labels = c("0", "0.25", "0.5", "≥0.75"),
    colours = RColorBrewer::brewer.pal(10, "RdYlBu") |> rev()
  ) +
  
  facet_grid(rows = vars(location),
             cols = vars(pathogen_label),
             scales = "free_x") +
  
  scale_y_continuous(breaks = 0:3, labels = c("1-7 days", "8-14 days", "15-21 days", "22-28 days")) +
  
  scale_x_date(date_breaks = "2 months", labels = scales::label_date_short(sep = "\n"),
               date_minor_breaks = "1 months",
               guide = guide_axis(minor.ticks = TRUE)) +
  
  ylab("Horizon (days ahead)") + xlab(NULL) +
  
  plot_theme_paper +
  
  theme(legend.position = "right",
        panel.background = element_facet_background) +
  
  ggtitle("Forecast performance by horizon (days ahead)")
p_heatmap


ggsave(
  "paper_figures/supp/f_perf_heatmap.png",
  width = 10, height = 8,
  bg = "white"
)


plot_data_skill |> 
  mutate(x_2 = date,
         x_1 = date + days(7),
         x_3 = date + days(7),
         x_4 = date + days(14),
         y_2 = days_ahead_group / 7 - 0.5,
         y_1 = days_ahead_group / 7 - 0.5,
         y_3 = days_ahead_group / 7 + 0.5,
         y_4 = days_ahead_group / 7 + 0.5) |> 
  expand_grid(corner = 1:4) |> 
  group_by(days_ahead_group, pathogen, location, round_id) |> 
  mutate(
    x = c(first(x_1), first(x_2), first(x_3), first(x_4))[corner],
    y = c(first(y_1), first(y_2), first(y_3), first(y_4))[corner],
  ) |> 
  
  filter(location == "VIC", pathogen == "RSV") |> 
  
  ggplot() +
  
  geom_polygon(aes(x, y, fill = crps_log, group = interaction(days_ahead_group, round_id))) +
  
  # geom_rect(aes(xmin = date, xmax = date + days(7),
  #               ymin = days_ahead_group / 7 - 0.4, ymax = days_ahead_group / 7 + 0.4, fill = crps_log,
  #               label = round_id),
  #           colour = "white", linewidth = 0.25) +
  
  
  
  # geom_rect(aes(xmin = date, xmax = date + days(7),
  #               ymin = days_ahead_group / 7 - 0.46, ymax = days_ahead_group / 7 - 0.4),
  #           fill = "black",
  #           colour = "white", linewidth = 0.25) +
  
  geom_line(aes(x = notification_date, y = y - 0.5, group = location), plot_data_counts |> filter(location == "VIC", pathogen == "RSV"),
            colour = "white", linewidth = 1.5, alpha = 0.8) +
  geom_line(aes(x = notification_date, y = y - 0.5, group = location), plot_data_counts |> filter(location == "VIC", pathogen == "RSV"),
            linewidth = 0.4) +
  
  
  scale_fill_gradientn(
    name = "CRPS", limits = c(0, 0.75),
    breaks = c(0, 0.25, 0.5, 0.75), labels = c("0", "0.25", "0.5", "≥0.75"),
    colours = RColorBrewer::brewer.pal(10, "RdYlBu") |> rev()
  ) +
  
  facet_grid(rows = vars(location),
             cols = vars(pathogen_label),
             scales = "free_x") +
  
  scale_y_continuous(breaks = 0:3, labels = c("1-7 days", "8-14 days", "15-21 days", "22-28 days")) +
  
  scale_x_date(date_breaks = "2 months", labels = scales::label_date_short(sep = "\n"),
               date_minor_breaks = "1 months",
               guide = guide_axis(minor.ticks = TRUE)) +
  
  ylab("Horizon (days ahead)") + xlab(NULL) +
  
  plot_theme_paper +
  
  theme(legend.position = "right",
        panel.background = element_facet_background) +
  
  ggtitle("Forecast performance by horizon (days ahead)")




