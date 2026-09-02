source("R/functions.R")
source("R/paper/common.R")

snaive_perf_data <- read_parquet("data/snaive_perf_extended.parquet")




perf_data <- read_perf("baseline") |> 
  bind_rows(snaive_perf_data) |> 
  process_perf_data()

combined_counts <- read_csv("data/combined-count-smoothed_paper.csv")



perf_data |> 
  filter(horizon == 7) |> 
  ggplot() + 
  geom_tile(aes(x = ymd(round_id), y = location, fill = .model),
            alpha = 0.5) +
  
  facet_grid(rows = vars(pathogen))

perf_data_paired <- perf_data |> 
  group_by(pathogen, location, round_id, date) |> 
  filter(n() == 2) |> 
  ungroup() |> 
  filter(horizon <= 7 * 4)

location_y <- tibble(
  location = unique(perf_data_paired$location)
) |> 
  mutate(y_offset = row_number())




plot_data_counts <- combined_counts |> 
  filter(notification_date >= ymd("2025-04-01") | (pathogen == "RSV" & notification_date >= ymd("2025-02-01")),
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
            skill_score = pmax(skill_score, -1),
            date = first(date)) |> 
  label_pathogen() |> 
  left_join(location_y, by = join_by(location)) |> 
  arrange(location, desc(date))

p_heatmap <- ggplot(plot_data_skill) +
  
  geom_rect(aes(xmin = date, xmax = date + days(7),
                ymin = y_offset - 0.4, ymax = y_offset + 0.4, fill = skill_score,
                label = round_id),
            colour = "white", linewidth = 0.25) +
  
  geom_rect(aes(xmin = date, xmax = date + days(7),
                ymin = y_offset - 0.46, ymax = y_offset - 0.4),
            fill = "black",
            colour = "white", linewidth = 0.25) +
  
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
  
  ylab("Region") + xlab(NULL) +
  
  plot_theme_paper +
  
  theme(legend.position = "right",
        panel.background = element_facet_background)

p_heatmap


ggsave(
  "paper_figures/supp/f_skill_heatmap.png",
  width = 10, height = 5,
  bg = "white"
)

