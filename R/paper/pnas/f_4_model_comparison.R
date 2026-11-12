source("R/functions.R")
source("R/paper/common.R")


perf_data <- bind_rows(
  read_perf("baseline"),
  read_perf("no_holidays"),
  read_perf("no_annual"),
  read_perf("no_holidays_no_annual")
) |> 
  process_perf_data()

combined_counts <- read_csv("data/combined-count-smoothed_final_anonymised.csv")




perf_data |> 
  filter(horizon == 7) |> 
  count(round_id = ymd(round_id), location, pathogen) |> 
  ggplot() + 
  geom_tile(aes(x = ymd(round_id), y = location, fill = factor(n))) +
  
  facet_grid(rows = vars(pathogen))

model_labels <- c(
  "baseline" = "With both holiday<br>effects and seasonal GP",
  "no_holidays" = "With seasonal GP",
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
  
  facet_wrap(~pathogen_label, ncol = 1) +
  
  theme(legend.position = "none", panel.grid.major.x = element_gridline,
        axis.text.y.left = element_markdown(size = 11, hjust = 0),
        panel.background = element_facet_background,
        panel.spacing.x = unit(0.5, "cm"))

p_skill


p_skill

ggsave(
  "pnas_figures/f_4_model_comparison.png",
  width = 8.7 * (9 / 12.5), height = 5,
  bg = "white"
)


