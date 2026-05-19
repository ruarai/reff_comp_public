

school_holidays <- get_school_holidays() |> 
  mutate(start_date_minus = end_date - days(14),
         end_date_plus = end_date + days(14)) |> 
  
  rowwise() |> 
  mutate(date = list(seq(start_date_minus, end_date_plus, by = "days"))) |> 
  unnest(date) |> 
  select(location, holiday_name, date) |> 
  mutate(is_near_holiday = TRUE)


perf_data_joined <- perf_data_paired |> 
  
  left_join(school_holidays, by = join_by(location, date)) |> 
  
  mutate(is_near_holiday = replace_na(is_near_holiday, FALSE))


perf_data_joined |> 
  group_by(pathogen, .model, is_near_holiday) |>
  summarise(crps_log = mean(crps_log)) |> 
  
  ggplot() +
  geom_linerange(aes(y = .model, xmin = 0, xmax = crps_log, colour = is_near_holiday),
                 position = position_dodge2(width = 0.25)) +
  
  facet_wrap(~pathogen) +
  
  theme_bw()
  