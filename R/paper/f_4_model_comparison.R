source("R/functions.R")
source("R/paper/common.R")


perf_data <- bind_rows(
  read_perf("baseline"),
  read_perf("no_holidays"),
  read_perf("no_annual"),
  read_perf("no_holidays_no_annual")
) |> 
  process_perf_data()

combined_counts <- read_csv("data/combined-count-smoothed-2026-04-22.csv")




perf_data |> 
  filter(horizon == 7) |> 
  count(round_id = ymd(round_id), location, pathogen) |> 
  ggplot() + 
  geom_tile(aes(x = ymd(round_id), y = location, fill = factor(n))) +
  
  facet_grid(rows = vars(pathogen))

model_labels <- c(
  "baseline" = "With both holiday<br>effects and annual GP",
  "no_holidays" = "With annual GP",
  "no_annual" = "With holiday effects",
  "no_holidays_no_annual" = "With neither"
)

n_models <- length(unique(perf_data$.model))

perf_data_paired <- perf_data |> 
  group_by(location, pathogen, round_id, horizon, date) |> 
  filter(n() == n_models)

plot_data_performance <- perf_data_paired |> 
  filter(horizon > 21, horizon <= 28) |> 
  group_by(.model, pathogen) |> 
  summarise(crps_log = mean(crps_log)) |> 
  label_pathogen() |> 
  mutate(model_label = factor(.model, names(model_labels), model_labels)) |> 
  group_by(pathogen) |> 
  mutate(skill_score = 1 - crps_log / crps_log[.model == "no_holidays_no_annual"]) |> 
  filter(.model != "no_holidays_no_annual")

p_skill <- ggplot(plot_data_performance) +
  
  geom_hline(aes(yintercept = model_label),
             colour = "grey70", linetype = "28", linewidth = 0.3) +
  
  geom_linerange(aes(y = model_label, xmin = 0, xmax = skill_score,
                     colour = skill_score < 0),
                 linewidth = 6.0) +
  
  geom_vline(xintercept = 0, linetype = "22") +
  
  scale_colour_brewer(type = "qual", palette = 2) +
  
  scale_x_continuous(breaks = c(-0.15, -0.1, -0.05, 0, 0.05, 0.1, 0.15)) +
  
  coord_cartesian(xlim = c(-0.15, 0.15)) +
  
  ylab(NULL) + xlab("Skill-score") +
  
  plot_theme_paper + 
  
  facet_wrap(~pathogen_label, ncol = 3) +
  
  theme(legend.position = "none", panel.grid.major.x = element_gridline,
        axis.text.y.left = element_markdown(size = 11, hjust = 0),
        panel.background = element_facet_background,
        panel.spacing.x = unit(0.5, "cm"))

p_skill


location_y <- tibble(
  location = unique(perf_data_paired$location)
) |> 
  mutate(y_offset = row_number())


plot_data_heatmap <- perf_data_paired |> 
  filter(horizon > 21, horizon <= 28) |> 
  group_by(.model, pathogen, location, round_id) |> 
  summarise(crps_log = mean(crps_log),
            date = first(date)) |> 
  label_pathogen() |> 
  mutate(model_label = factor(.model, names(model_labels), model_labels)) |> 
  group_by(pathogen, location, date) |> 
  mutate(skill_score = 1 - crps_log / crps_log[.model == "no_holidays_no_annual"]) |> 
  filter(.model == "baseline") |> 
  mutate(skill_score = pmax(-1, skill_score)) |> 
  
  left_join(location_y) |> 
  arrange(location, desc(date))


plot_data_counts <- combined_counts |> 
  filter(notification_date >= ymd("2025-04-01") | (pathogen == "RSV" & notification_date >= ymd("2025-02-01")),
         notification_date <= ymd("2026-01-01")) |> 
  inner_join(location_y, by = join_by(location)) |> 
  group_by(location, pathogen) |> 
  mutate(y = y_offset + scales::rescale(value_trend, to = c(0.15, 0.9)),
         y_full = scales::rescale(value_trend)) |> 
  label_pathogen()



p_heatmap <- ggplot(plot_data_heatmap) +
  
  
  geom_rect(aes(xmin = date, xmax = date + days(7),
                ymin = y_offset - 0.4, ymax = y_offset + 0.4, fill = skill_score),
            colour = "white", linewidth = 0.25) +
  
  
  geom_rect(aes(xmin = date, xmax = date + days(7),
                ymin = y_offset - 0.46, ymax = y_offset - 0.4),
            fill = "black",
            colour = "white", linewidth = 0.25) +
  
  
  geom_line(aes(x = notification_date, y = y - 0.5, group = location), plot_data_counts,
            colour = "white", linewidth = 1.5, alpha = 0.8) +
  geom_line(aes(x = notification_date, y = y - 0.5, group = location), plot_data_counts,
            linewidth = 0.4) +
  
  facet_wrap(~pathogen_label, ncol = 3, scales = "free_x") +
  
  scale_y_continuous(breaks = location_y$y_offset, labels = location_y$location) +
  
  scale_x_date(date_breaks = "2 months", labels = scales::label_date_short(sep = "\n"),
               date_minor_breaks = "1 months",
               guide = guide_axis(minor.ticks = TRUE)) +
  
  ylab("Region") + xlab(NULL) +
  
  scale_fill_distiller(type = "div", direction = 1, palette = "RdBu",
                       limits = c(-1, 1), name = "Skill score",
                       breaks = c(-1, 0, 1), labels = c("≤-1", "0", "1")) +
  
  
  plot_theme_paper +
  
  theme(legend.position = "left")

# Skill score (relative to null model)
# "Forecast skill (model with both holiday effects and annual GP relative to null model)"


p_top <- p_skill + 
  plot_annotation("<b>A</b>. 4th week forecast skill relative to null model",
                  theme = plot_theme_annotation)

p_bottom <- p_heatmap + 
  plot_annotation("<b>B</b>. 4th week forecast skill (model with both holiday effects and annual<br>GP relative to null model) by epidemic period",
                  theme = plot_theme_annotation)


(wrap_elements(full = p_top) / 
    wrap_elements(full = p_bottom)) +
  plot_layout(heights = c(1, 1.5))


p_top

ggsave(
  "paper_figures/f_4_model_comparison.png",
  width = 9, height = 3,
  bg = "white"
)


