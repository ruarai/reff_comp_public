

source("R/functions.R")

source("R/presentation/plot_theme.R")

snaive_perf_data <- read_parquet("data/snaive_perf_extended.parquet")



perf_data <- read_perf("baseline") |> 
  bind_rows(snaive_perf_data) |> 
  process_perf_data()


acefa_perf_data <- read_rds("data/case_forecast_measures.rds") |> 
  filter(measure_name == "crps", target == "case_incidence") |> 
  mutate(pathogen = as.character(if_else(pathogen == "Influenza", "flu", pathogen)),
         round_id = as.character(round_id)) |> 
  select(-c(measure_name, target)) |> 
  pivot_wider(names_from = measure_scale) |> 
  rename(crps = natural, crps_log = log)


perf_data_full <- bind_rows(
  perf_data, acefa_perf_data
) |> 
  select(.model, location, pathogen, round_id, origin_date, date, crps, crps_log) |> 
  mutate(horizon = as.numeric(date - origin_date)) |> 
  filter(horizon > 0, horizon <= 28, location != "NZ") |> 
  
  mutate(
    .model = factor(
      .model, 
      c("snaive", "dstg-RV_model", "UQ-DD", "UQ-renewal_global", 
        "UQ-renewal_i", "uom-multiwalk", "ensemble-mixture", "ensemble-quantile", 
        "uom-reff-comp")
    )
  )



eval_counts <- read_csv("data/combined-count-smoothed-2026-01-15.csv")

make_two_panel_plot <- function(x_location, x_pathogen) {
  
  p_skill <- perf_data_full |> 
    filter(horizon > 21, horizon <= 28) |>
    filter(location == x_location, pathogen == x_pathogen) |> 
    
    group_by(location, pathogen, round_id, date) |> 
    mutate(crps_log_snaive = crps_log[.model == "snaive"]) |> 
    filter(.model == "uom-reff-comp") |> 
    
    group_by(location, pathogen, round_id) |> 
    
    summarise(skill_score = 1 - sum(crps_log) / sum(crps_log_snaive),
              date = first(date), .groups = "drop") |> 
    
    ggplot() +
    geom_hline(yintercept = 0, linetype = 22, linewidth = 0.5) +
    geom_linerange(aes(x = date, ymin = 0, ymax = skill_score, colour = skill_score > 0),
                   linewidth = 2.0) +
    
    ggokabeito::scale_colour_okabe_ito(order = c(6, 3)) +
    
    scale_y_continuous(trans = log1minus_trans(),
                       breaks = c(-3, -1, 0, 0.5, 0.75)) +
    
    scale_x_date(date_breaks = "month",
                 label = scales::label_date_short()) +
    
    coord_cartesian(ylim = c(-3.1, 0.8),
                    xlim = c(ymd("2025-05-01"), ymd("2025-11-01"))) +
    
    plot_theme_pres +
    xlab(NULL) + ylab("Skill score") +
    
    theme(legend.position = "none",
          panel.grid.major.y = element_blank())
  
  
  
  p_trend <- ggplot(eval_counts |> filter(location == x_location, pathogen == x_pathogen, notification_date >= ymd("2025-01-01"))) +
    geom_line(aes(x = notification_date, y = cases_trend)) +
    
    coord_cartesian(xlim = c(ymd("2025-05-01"), ymd("2025-11-01"))) +
    
    xlab(NULL) + ylab("Cases (trend)") +
    ggtitle(str_c(x_location, " ", x_pathogen)) +
    
    scale_x_date(date_breaks = "month",
                 label = scales::label_date_short()) +
    
    plot_theme_pres +
    theme(panel.grid.major.y = element_blank(),
          axis.text.x = element_blank())
  
  p_trend / p_skill
}

make_two_panel_plot("VIC", "flu") | 
  make_two_panel_plot("VIC", "RSV")


make_two_panel_plot("WA", "flu") | 
  make_two_panel_plot("WA", "RSV") |
  make_two_panel_plot("WA", "SARSCOV2")


make_two_panel_plot("QLD", "flu") | 
  make_two_panel_plot("QLD", "RSV") |
  make_two_panel_plot("QLD", "SARSCOV2")


perf_data_full |> 
  filter(horizon > 21, horizon <= 28) |>
  
  group_by(location, pathogen, round_id, date) |> 
  mutate(crps_log_snaive = crps_log[.model == "snaive"]) |> 
  filter(.model == "uom-reff-comp") |> 
  
  group_by(location, pathogen) |> 
  
  summarise(skill_score = 1 - sum(crps_log) / sum(crps_log_snaive),
            date = first(date), .groups = "drop") |> 
  
  ggplot() +
  
  geom_linerange(aes(y = location, xmin = 0, xmax = skill_score,
                     colour = skill_score > 0),
                 linewidth = 7.0) +
  
  ggokabeito::scale_colour_okabe_ito(order = c(3)) +
  
  scale_x_continuous(trans = log1minus_trans(),
                     breaks = c(-3, -1, 0, 0.5, 0.75)) +
  
  plot_theme_pres_2 +
  
  facet_wrap(~pathogen) +
  
  theme(legend.position = "none") +
  
  ylab(NULL) + xlab("Skill score")

perf_data_full |> 
  
  group_by(location, pathogen, round_id, date, horizon) |> 
  mutate(crps_log_snaive = crps_log[.model == "snaive"]) |> 
  filter(.model == "uom-reff-comp") |> 
  
  group_by(location, pathogen, horizon) |> 
  
  summarise(skill_score = 1 - sum(crps_log) / sum(crps_log_snaive),
            date = first(date), .groups = "drop") |> 
  
  ggplot() +
  
  geom_linerange(aes(x = horizon, ymin = 0, ymax = skill_score,
                     colour = skill_score > 0),
                 linewidth = 1.0) +
  
  geom_hline(yintercept = 0, linetype = "44") +
  
  ggokabeito::scale_colour_okabe_ito(order = c(6, 3)) +
  
  scale_y_continuous(trans = log1minus_trans(),
                     breaks = c(-3, -1, 0, 0.5, 0.75)) +
  
  plot_theme_pres_2 +
  
  facet_grid(rows = vars(pathogen), cols = vars(location)) +
  
  theme(legend.position = "none") +
  
  xlab("Horizon (days ahead)") + ylab("Skill score")

