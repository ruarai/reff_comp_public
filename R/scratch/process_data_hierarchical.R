
process_counts_hierarchical <- function(
    combined_counts,
    x_pathogen, x_target,
    n_days_ahead,
    origin_dates
) {
  training_dates <- get_training_dates()
  
  # Filter the counts for our given location and pathogen
  # and drop data prior to the training (cutoff) date
  reduced_counts <- combined_counts |> 
    left_join(training_dates) |>
    
    filter(notification_date >= cutoff_date) |>
    filter(pathogen == x_pathogen, target == x_target) |> 
    filter(is_complete)
  
  
  # Drop the last day of data, but only for forecasts where we had
  # issues with right truncation often (i.e. pre-2026 season)
  if(any(origin_dates$origin_date < ymd("2026-01-01"))) {
    reduced_counts <- reduced_counts |> 
      filter(notification_date < max(notification_date))
  }
  
  min_origin_date <- reduced_counts |> group_by(location) |> summarise(d = max(notification_date)) |> pull(d) |> min()
  
  training_counts <- reduced_counts |> 
    select(location, pathogen, target, date = notification_date, 
           y = value, y_trend = value_trend, y_trend_and_remainder = value_trend_and_remainder) |> 
    filter(date <= min_origin_date) |> 
    group_by(location) |> 
    complete(date = seq(min(date), max(date), by = "days"),
             fill = list(pathogen = x_pathogen, target = x_target))
  
  # Make sure all values are complete
  n_na <- sum(is.na(training_counts$y))
  stopifnot(n_na == 0)
  
  locations <- unique(training_counts$location)
  
  # Bind extra rows that make up the forecasting horizon
  if(n_days_ahead > 0) {
    
    end_date <- max(origin_dates$origin_date) + days(n_days_ahead)
    
    training_counts <- training_counts |> 
      ungroup() |> 
      group_split(location) |> 
      map(function(x) {
        bind_rows(
          x,
          expand_grid(
            date = seq(max(x$date) + days(1), end_date, by = "days"),
            location = x$location[[1]],
            pathogen = x_pathogen, target = x_target
          )
        )
      }) |> 
      bind_rows()
  }
  
  school_holidays_location <- get_school_holidays_long() |> 
    filter(location %in% locations)
  
  public_holidays_location <- get_public_holidays() |> 
    filter(location %in% locations) |> 
    mutate(public_holiday_type = match(holiday_name, unique(holiday_name))) |> 
    select(date, location, public_holiday_type)
  
  # Add a bunch of the effects data (day of week, day of year, holidays, ...)
  training_counts_full <- training_counts |>
    arrange(location, date) |> 
    mutate(
      x = as.numeric(date - min(date)),
      day_of_week = wday(date, week_start = 7),
      day_of_year = yday(date)
    ) |> 
    left_join(school_holidays_location, by = join_by(date, location)) |> 
    left_join(public_holidays_location, by = join_by(date, location)) |> 
    mutate(school_holiday_type = replace_na(school_holiday_type, 0) + 1,
           public_holiday_type = replace_na(public_holiday_type, 0) + 1) |> 
    
    left_join(origin_dates, by = join_by(location))
  
  return(training_counts_full)
}


create_stan_data_hierarchical <- function(
    training_counts,
    n_days_ahead,
    x_pathogen
) {
  gi <- get_gi(x_pathogen)
  
  training_counts_split <- training_counts |> 
    ungroup() |> group_split(location)
  
  n_total <- nrow(training_counts_split[[1]])
  n_days_data <- nrow(training_counts_split[[1]] |> drop_na(y))
  
  N_public_holiday_types <- max(training_counts$public_holiday_type) - 1
  N_school_holiday_types <- max(training_counts$school_holiday_type) - 1
  
  N_loc <- length(unique(training_counts$location))
  
  data_list <- list(
    x = training_counts_split[[1]]$x,
    day_of_week = training_counts_split[[1]]$day_of_week,
    day_of_year = training_counts_split[[1]]$day_of_year,
    school_holiday_type = do.call(rbind, training_counts_split |> map(\(x) x$school_holiday_type)),
    public_holiday_type = do.call(rbind, training_counts_split |> map(\(x) x$public_holiday_type)),
    
    y = do.call(rbind, training_counts_split |> map(\(x) x$y[1:n_days_data])),
    N = n_total, N_data = n_days_data,
    N_public_holiday_types = N_public_holiday_types,
    N_school_holiday_types = N_school_holiday_types,
    
    N_loc = N_loc,
    
    M_annual = 40,
    M_noise = calculate_M_matern(c = 1.25, length_scale = 10.0, x = training_counts$x),
    c_factor = 1.25,
    
    gi_reversed = rev(gi)
  )
  
  return(data_list)
}