Rt_draws <- read_parquet("data/paper/Rt_draws.parquet")

combined_counts <- read_csv("../forecast-hub/processed-data/2026/2026-08-19/combined-count-smoothed-2026-08-19.csv") |> 
  rename(date = notification_date)


get_pop_size <- function(location) { # 2025 ABS ERP
  case_when(
    location == "VIC" ~ 7.121 * 1e6,
    location == "NSW" ~ 8.641 * 1e6,
    location == "QLD" ~ 5.712 * 1e6,
    location == "SA" ~ 1.910 * 1e6,
    location == "WA" ~ 3.076 * 1e6,
    location == "TAS" ~ 0.579 * 1e6,
    location == "NT" ~ 0.267 * 1e6,
    location == "ACT" ~ 0.487 * 1e6
  )
}

plot_data <- combined_counts |> 
  expand_grid(lag = c(50, 100, 150, 200, 250, 300, 350)) |> 
  group_by(location, pathogen, lag) |> 
  mutate(cumulative_cases = cumsum(value_trend) - lag(cumsum(value_trend), lag[1], default = 0)) |> 
  ungroup() |> 
  filter(pathogen == "RSV", location == "VIC") |> 
  left_join(Rt_draws) |> 
  
  mutate(pop_size = get_pop_size(location)) |> 
  
  mutate(cumulative_frac = (cumulative_cases / pop_size),
         location = factor(location),
         log_R_holidayless = log(R_holidayless),
         day_of_year = yday(date)) |> 
  
  drop_na(log_R_holidayless)


ggplot(plot_data) +
  geom_path(aes(x = cumulative_frac, y = exp(log_R_holidayless), colour = day_of_year)) +
  geom_point(aes(x = cumulative_frac, y = exp(log_R_holidayless), colour = day_of_year),
             size = 0.1) +
  
  xlab("Fraction of population infected (and detected) over past N days") +
  
  ylab("R(t)") +
  
  scale_colour_viridis_c() +
  
  facet_wrap(~lag) +
  
  plot_theme_pres
