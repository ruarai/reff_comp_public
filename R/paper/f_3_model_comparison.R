source("R/functions.R")
source("R/paper/common.R")


perf_data <- bind_rows(
  read_parquet("saved_results/baseline_2_performance.parquet") |> 
    mutate(.model = "baseline"),
  read_parquet("saved_results/no_annual_performance.parquet") |> 
    collect() |> mutate(.model = "no_annual"),
  read_parquet("saved_results/no_holidays_performance.parquet") |> 
    collect() |> mutate(.model = "no_holidays"),
  read_parquet("saved_results/no_holidays_no_annual_performance.parquet") |>
    collect() |> mutate(.model = "no_holidays_no_annual"),
  # open_dataset(fs::dir_ls("saved_results/no_delay/", recurse = TRUE, glob = "*performance.parquet")) |> 
  #   collect() |> mutate(.model = "no_delay"),
  open_dataset(fs::dir_ls("saved_results/no_school_holidays/", recurse = TRUE, glob = "*performance.parquet")) |> 
    collect() |> mutate(.model = "no_school_holidays"),
  open_dataset(fs::dir_ls("saved_results/no_public_holidays/", recurse = TRUE, glob = "*performance.parquet")) |> 
    collect() |> mutate(.model = "no_public_holidays")
) |> 
  process_perf_data()


perf_data |> 
  filter(horizon == 7) |> 
  count(round_id = ymd(round_id), location, pathogen) |> 
  ggplot() + 
  geom_tile(aes(x = ymd(round_id), y = location, fill = n)) +
  
  facet_grid(rows = vars(pathogen))


model_labels <- tribble(
  ~.model,                ~eff_school_hol, ~eff_public_hol, ~eff_ann, 
  "baseline",             TRUE,            TRUE,            TRUE,     
  "no_school_holidays",   FALSE,           TRUE,            TRUE,     
  "no_public_holidays",   TRUE,            FALSE,           TRUE,     
  "no_holidays",          FALSE,           FALSE,           TRUE,     
  "no_annual",            TRUE,            TRUE,            FALSE,    
  "no_holidays_no_annual",FALSE,           FALSE,           FALSE,    
  # "no_delay",             TRUE,            TRUE,            TRUE,   
)

effect_labels <- c(
  "delay" = "Reporting delay",
  "ann" = "Periodic GP",
  "school_hol" = "School holidays",
  "public_hol" = "Public holidays"
)

plot_data_labels <- model_labels |> 
  pivot_longer(starts_with("eff_"), names_prefix = "eff_") |> 
  mutate(name = factor(name, names(effect_labels), effect_labels),
         .model = factor(.model, model_labels$.model)) |>
  
  expand_grid(pathogen = c("RSV", "flu", "SARSCOV2")) |> 
  label_pathogen()

plot_data_models <- plot_data_labels |> 
  distinct(.model) |> 
  mutate(colour = row_number() %% 2 == 0)

p_labels <- ggplot() +
  
  geom_linerange(aes(x = .model, ymin = -Inf, ymax = Inf, alpha = colour), plot_data_models,
                 linewidth = 12.0) +
  scale_alpha_manual(values = c(0, 0.05)) +
  
  geom_point(aes(y = name, x = .model),
             size = 5,
             plot_data_labels |> filter(value)) +
  
  geom_point(aes(y = name, x = .model),
             size = 3, pch = 4, stroke = 1.5,
             plot_data_labels |> filter(!value)) +
  
  
  
  geom_blank(aes(y = name, x = .model), plot_data_labels) +
  
  xlab(NULL) + ylab("Feature") +
  
  facet_wrap(~pathogen, ncol = 3) +
  
  plot_theme_paper +
  theme(legend.position = "none",
        strip.text = element_markdown(colour = "white"),
        panel.grid.major.x = element_gridline,
        axis.text.x = element_blank())

p_labels


n_models <- length(unique(perf_data$.model))

perf_data_paired <- perf_data |> 
  group_by(location, pathogen, round_id, horizon, date) |> 
  filter(n() == n_models)

p_crps <- perf_data_paired |> 
  filter(horizon > 21, horizon <= 28) |> 
  group_by(.model, pathogen) |> 
  summarise(crps_log = mean(crps_log)) |> 
  label_pathogen() |> 
  mutate(.model = factor(.model, model_labels$.model)) |> 
  ggplot() +
  
  
  geom_linerange(aes(x = .model, ymin = -Inf, ymax = Inf, alpha = colour), plot_data_models,
                 linewidth = 12.0) +
  scale_alpha_manual(values = c(0, 0.05)) +
  
  geom_linerange(aes(x = .model, ymin = 0, ymax = crps_log,
                     colour = .model == "Baseline"),
                 linewidth = 6.0) +
  
  ggokabeito::scale_colour_okabe_ito(order = c(9, 5)) +
  
  xlab(NULL) + ylab("CRPS (log-counts)") +
  
  plot_theme_paper + 
  
  facet_wrap(~pathogen_label, ncol = 3) +
  
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        panel.grid.major.x = element_gridline) +
  
  ggtitle("Fourth week forecasting performance")

p_skill <- perf_data_paired |> 
  filter(horizon > 21, horizon <= 28) |> 
  group_by(.model, pathogen) |> 
  summarise(crps_log = mean(crps_log)) |> 
  group_by(pathogen) |> 
  mutate(skill_score = 1 - crps_log / crps_log[.model == "baseline"]) |> 
  
  label_pathogen() |> 
  mutate(.model = factor(.model, model_labels$.model)) |> 
  ggplot() +
  
  geom_linerange(aes(x = .model, ymin = -Inf, ymax = Inf, alpha = colour), plot_data_models,
                 linewidth = 12.0) +
  scale_alpha_manual(values = c(0, 0.05)) +
  
  geom_linerange(aes(x = .model, ymin = 0, ymax = skill_score,
                     colour = skill_score < 0),
                 linewidth = 6.0) +
  
  geom_hline(yintercept = 0, linetype = "22") +
  
  scale_colour_brewer(type = "qual", palette = 2) +
  
  xlab(NULL) + ylab("Skill-score") +
  
  plot_theme_paper + 
  
  facet_wrap(~pathogen_label, ncol = 3) +
  
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        panel.grid.major.x = element_gridline) +
  
  ggtitle("Skill score (relative to full model)")

p_skill

(p_skill / p_labels) +
  plot_layout(heights = c(2, 1))

