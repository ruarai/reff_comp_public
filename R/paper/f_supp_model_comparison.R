source("R/functions.R")
source("R/paper/common.R")


perf_summary <- read_perf_summary(c(
  "baseline", "no_public_holidays", "no_school_holidays", "no_holidays",
  "no_asc_holidays", "no_Rt_holidays", "no_delay", "no_annual", "no_holidays_no_annual"
))


perf_summary |> 
  filter(horizon_week == 0) |> 
  count(round_id = ymd(round_id), location, pathogen) |> 
  ggplot() + 
  geom_tile(aes(x = ymd(round_id), y = location, fill = factor(n))) +
  
  facet_grid(rows = vars(pathogen))

model_labels <- c(
  # "baseline" = "Full model",
  
  "label_hol_component" = "<b>Holiday effect component</b>",
  
  "no_asc_holidays" = "No holiday reporting effects",
  "no_Rt_holidays" = "No holiday R(t) effects",
  
  "label_hol_type" = "<b>Type of holiday</b>",
  
  "no_school_holidays" = "No school holiday effects",
  "no_public_holidays" = "No public holiday effects",
  
  "label_delay" = "<b>Infection-reporting delay</b>",
  
  "no_delay" = "No delay"
  
  # "no_holidays" = "No holiday reporting or R(t) effects" # TODO remove
  
  # "label_periodic" = "<b>Periodic GP</b>",
  
  # "no_annual" = "No periodic GP"
  
  # "no_holidays_no_annual" = "No holiday effects or periodic GP"
) |> rev()

plot_data_label <- tibble(model_label = factor(model_labels, model_labels)) |> 
  expand_grid(pathogen = c("RSV", "flu", "SARSCOV2")) |> 
  label_pathogen()

plot_data_performance <- perf_summary |> 
  filter(.model %in% c("baseline", names(model_labels))) |> 
  filter(horizon_week == 3) |> 
  group_by(.model, pathogen) |> 
  summarise(crps_log = mean(crps_log)) |> 
  label_pathogen() |> 
  mutate(model_label = factor(.model, names(model_labels), model_labels)) |> 
  group_by(pathogen) |> 
  mutate(skill_score = 1 - crps_log / crps_log[.model == "baseline"]) |> 
  filter(.model != "baseline")

p_skill <- ggplot(plot_data_performance) +
  
  geom_blank(aes(x = 0, y = model_label),
             plot_data_label) +
  
  geom_hline(aes(yintercept = model_label),
             colour = "grey70", linetype = "28", linewidth = 0.3) +
  geom_rug(aes(y = model_label), length = unit(0.1, "cm")) +
  
  geom_linerange(aes(y = model_label, xmin = 0, xmax = skill_score,
                     colour = skill_score < 0),
                 linewidth = 4.0) +
  
  geom_vline(xintercept = 0, linetype = "22") +
  
  scale_colour_brewer(type = "qual", palette = 2) +
  
  coord_cartesian(xlim = c(-0.05, 0.05)) +
  
  ylab(NULL) + xlab("Skill-score") +
  
  plot_theme_paper + 
  
  facet_wrap(~pathogen_label, ncol = 1, scales = "free_y") +
  
  theme(legend.position = "none", panel.grid.major.x = element_gridline,
        axis.text.y.left = element_markdown(size = 11, hjust = 0), axis.ticks.y = element_blank()) +
  
  ggtitle("Skill score (relative to full model)")

p_skill


ggsave(
  "paper_figures/supp/f_model_comparison.png",
  width = 9, height = 6,
  bg = "white"
)


