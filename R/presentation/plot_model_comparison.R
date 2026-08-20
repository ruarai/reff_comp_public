
source("R/functions.R")
source("R/presentation/plot_theme.R")

perf_data <- bind_rows(
  read_perf("baseline"),
  read_perf("no_holidays"),
  read_perf("no_annual"),
  read_perf("no_holidays_no_annual")
) |> 
  process_perf_data()


perf_data |> 
  filter(horizon > 21, horizon <= 28) |> 
  # group_by(.model, location, pathogen, date) |> 
  # summarise(crps_log = mean(crps_log)) |> 
  
  ggplot() +
  geom_point(aes(x = date, y = crps_log, colour = .model),
             size = 0.5) +
  
  facet_grid(rows = vars(location), cols = vars(pathogen), scales = "free_y") +
  
  theme_bw()

plot_data <- perf_data |> 
  filter(horizon > 21, horizon <= 28) |> 
  group_by(round_id, location, pathogen, horizon) |> 
  mutate(crps_log_null = crps_log[[which(.model == "no_holidays_no_annual")]]) |> 
  
  filter(.model != "no_holidays_no_annual") |> 
  group_by(pathogen, .model) |> 
  summarise(skill_score = 1 - sum(crps_log) / sum(crps_log_null)) |> 
  label_pathogen() |> 
  
  mutate(
    .model = factor(
      .model,
      c("baseline", "no_holidays", "no_annual"),
      c("With both", "With annual periodic GP", "With holiday effects")
    )
  )
  
ggplot(plot_data) +
  
  geom_vline(xintercept = 0, linetype = "44", linewidth = 0.75) +
  
  geom_linerange(aes(xmin = 0, xmax = skill_score, y = .model, colour = skill_score > 0),
                 position = position_dodge(width = 0.75),
                 linewidth = 8.0) +
  
  ggokabeito::scale_colour_okabe_ito(order = c(3, 6)) +
  ylab(NULL) + xlab("Skill score") +
  
  plot_theme_pres +
  
  facet_wrap(~pathogen, ncol = 1) +
  
  coord_cartesian(xlim = c(-0.15, 0.15)) +
  
  theme(legend.position = "none",
        plot.title = element_text(size = 18)) +
  
  ggtitle("Forecasting skill\n(against a model without either component)")



perf_data |> 
  filter(horizon > 21, horizon <= 28) |> 
  group_by(.model, pathogen) |> 
  summarise(crps_log = mean(crps_log)) |> 
  
  mutate(
    .model = factor(
      .model,
      c("baseline", "no_holidays", "no_annual", "no_holidays_no_annual"),
      c("With both", "With annual periodic term", "With holiday effects", "Without both")
    )
  ) |> 
  
  ggplot() +
  
  geom_linerange(aes(xmin = 0, xmax = crps_log, y = .model, colour = .model == "With both"),
                 linewidth = 8.0) +
  # geom_point(aes(y = .model, x = crps_log, colour = .model == "With both"),
  #            size = 8) +
  # geom_point(aes(y = .model, x = crps_log),
  #            size = 3, colour = "white") +
  
  ylab(NULL) + xlab("CRPS") +
  
  ggokabeito::scale_colour_okabe_ito(order = c(9, 2)) +
  coord_cartesian(xlim = c(0.0, 0.45)) +
  # scale_x_continuous(breaks = c(0.25, 0.3, 0.35)) +
  
  ggtitle("Mean 21–28 day CRPS\n(lower is better)") +
  
  facet_wrap(~pathogen, ncol = 1) +
  
  plot_theme_pres +
  
  theme(legend.position = "none")

