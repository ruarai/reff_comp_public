



public_holidays <- get_public_holidays()



combined_counts <- read_csv("data/combined-count-smoothed-2026-04-15.csv")


public_holidays |> 
  expand_grid(pathogen = c("RSV", "SARSCOV2", "flu")) |> 
  filter(location == "NSW") |> 
  
  filter(!(holiday_name == "Easter" & lag(holiday_name) == "Easter"),
         !str_detect(holiday_name, "Observed")) |> 
  
  expand_grid(days_around = -14:14) |> 
  mutate(date_offset = date + days_around) |> 
  
  left_join(combined_counts |> select(date_offset = notification_date, location, pathogen, value, value_trend, value_trend_and_remainder)) |> 
  
  group_by(location, pathogen, date, holiday_name) |> 
  mutate(y = value_trend_and_remainder / mean(value_trend_and_remainder)) |> 
  
  ggplot() +
  
  geom_line(aes(x = days_around, y = y, group = interaction(date, pathogen))) +
  
  geom_vline(xintercept = 0, linetype = "44") +
  
  scale_y_log10() +
  coord_cartesian(ylim = c(0.25, 4.0)) +
  
  facet_wrap(~holiday_name) +
  
  theme_bw()


public_holidays |> 
  expand_grid(pathogen = c("RSV", "SARSCOV2", "flu")) |> 
  
  filter(!(holiday_name == "Easter" & lag(holiday_name) == "Easter"),
         !str_detect(holiday_name, "Observed"),
         location != "NZ") |> 
  
  expand_grid(days_around = -14:14) |> 
  mutate(date_offset = date + days_around) |> 
  
  left_join(combined_counts |> select(date_offset = notification_date, location, pathogen, value, value_trend, value_trend_and_remainder)) |> 
  
  group_by(location, pathogen, date, holiday_name) |> 
  mutate(y = value_trend_and_remainder / mean(value_trend_and_remainder)) |> 
  group_by(holiday_name, days_around) |> 
  summarise(y = median(y, na.rm = TRUE)) |> 
  
  ggplot() +
  
  geom_line(aes(x = days_around, y = y)) +
  
  geom_vline(xintercept = 0, linetype = "44") +
  
  xlab(NULL) + ylab("Relative case incidence (trend + remainder)") +
  
  scale_y_log10() +
  coord_cartesian(ylim = c(0.25, 4.0)) +
  
  facet_wrap(~holiday_name) +
  
  theme_bw()
