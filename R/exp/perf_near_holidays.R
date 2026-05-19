
near_width <- 0

school_holidays <- get_school_holidays() |> 
  mutate(start_date_minus = end_date - days(near_width),
         end_date_plus = end_date + days(near_width)) |> 
  
  rowwise() |> 
  mutate(date = list(seq(start_date_minus, end_date_plus, by = "days"))) |> 
  unnest(date) |> 
  select(location, holiday_name, date) |> 
  mutate(is_near_holiday = TRUE)


perf_data_joined <- perf_data_paired |> 
  
  left_join(school_holidays, by = join_by(location, date)) |> 
  
  mutate(is_in_holiday = replace_na(is_near_holiday, FALSE))


perf_data_joined |> 
  filter(.model %in% c("baseline", "no_school_holidays", "no_Rt_holidays", "no_asc_holidays")) |> 
  group_by(pathogen, .model, is_in_holiday) |>
  summarise(crps_log = mean(crps_log)) |> 
  
  group_by(pathogen, is_in_holiday) |> 
  mutate(skill_score = 1 - crps_log / crps_log[.model == "baseline"]) |> 
  filter(.model != "baseline") |> 
  
  ggplot() +
  
  geom_vline(xintercept = 0, linetype = "44") +
  
  geom_linerange(aes(y = .model, xmin = 0, xmax = skill_score, colour = is_in_holiday),
                 position = position_dodge2(width = 0.5),
                 linewidth = 8.0) +
  
  facet_wrap(~pathogen) +
  
  theme_bw()
  