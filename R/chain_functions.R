
summarise_holiday_effects <- function(fit, public_holidays) {
  school_holiday_draws <- fit |> 
    spread_draws(school_holiday_asc_effect[school_holiday_type],
                 school_holiday_Rt_effect[school_holiday_type]) |> 
    mutate(school_holiday_type = school_holiday_type - 1) |> 
    filter(school_holiday_type > 0) |> 
    mutate(holiday_name = case_when(
      school_holiday_type == 1 ~ "Summer holidays",
      school_holiday_type == 2 ~ "Mid-year holidays"
    )) |> 
    ungroup() |> 
    mutate(holiday_name = factor(holiday_name, c("Summer holidays", "Mid-year holidays")))
  
  
  public_holiday_labels <- public_holidays |> 
    distinct(holiday_name) |> 
    mutate(public_holiday_type = row_number())
  
  public_holiday_draws <- fit |> 
    spread_draws(public_holiday_asc_effect[public_holiday_type],
                 public_holiday_Rt_effect[public_holiday_type]) |> 
    mutate(public_holiday_type = public_holiday_type - 1) |> 
    filter(public_holiday_type > 0) |> 
    
    left_join(public_holiday_labels) |> 
    
    ungroup() |> 
    mutate(holiday_name = factor(holiday_name, public_holiday_labels$holiday_name))
  
  
  school_holiday_draws |> 
    pivot_longer(c(school_holiday_asc_effect, school_holiday_Rt_effect)) |> 
    group_by(holiday_name, name) |> 
    summarise(q95_lower = quantile(exp(value), 0.01 / 2),
              median = median(exp(value)),
              q95_upper = quantile(exp(value), 1 - 0.01 / 2))
  
  
}