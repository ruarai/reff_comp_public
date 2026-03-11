
fit_model <- function(data_list) {
  require(cmdstanr)
  model <- cmdstan_model("src/gp_renewal.stan")
  
  fit <- model$sample(
    data = data_list,
    seed = 1,
    chains = 4,
    parallel_chains = 4,
    refresh = 50,
    iter_warmup = 1000, iter_sampling = 1000,
    init = 0.5
  )
  
  return(fit)
}


get_origin_date <- function(combined_counts, x_pathogen, x_location) {
  combined_counts |> 
    filter(location == x_location, pathogen == x_pathogen) |> 
    pull(notification_date) |> 
    max()
}

write_forecast_file <- function(round_id, x_pathogen, x_location, origin_date) {
  round_results <- fs::path("results", round_id)
  
  fit <- read_rds(str_c(round_results, "/", x_pathogen, "_", x_location, ".rds"))
  training_counts <- read_parquet(str_c(round_results, "/", x_pathogen, "_", x_location, "_training_counts.parquet"))
  
  case_trajectories <- fit |> 
    spread_draws(cases[ix_x], inf_delayed[ix_x]) |> 
    mutate(date = training_counts$date[ix_x])
  
  
  forecast_data <- case_trajectories |> 
    ungroup() |> 
    select(output_type_id = .draw,
           notification_date = date,
           case_incidence = cases,
           case_incidence_smoothed = inf_delayed) |> 
    pivot_longer(c(case_incidence, case_incidence_smoothed),
                 names_to = "target") |> 
    mutate(
      horizon = as.integer(notification_date - origin_date),
      origin_date = origin_date,
      round_id = ymd(round_id),
      output_type = "sample",
      location = x_location,
      pathogen = x_pathogen,
      
      # We get NA values sometimes where rng freaks out (high variance)
      # replace with zero as a stop-gap.
      value = as.integer(replace_na(value, 0)),
    ) |> 
    
    filter(horizon >= -7, horizon <= 7 * 12) |> 
    
    select(
      round_id, origin_date,
      target, horizon, location, pathogen,
      
      output_type, output_type_id,
      
      value
    )
  
  output_file <- str_c(round_results, "/", x_pathogen, "_", x_location, "_uom-gp-reff.parquet")
  
  arrow::write_parquet(forecast_data, output_file)
}




process_counts <- function(
    combined_counts,
    x_location, x_pathogen,
    n_days_ahead
) {
  training_dates <- get_training_dates()
  
  # Filter the counts for our given location and pathogen
  # and drop data prior to the training date
  reduced_counts <- combined_counts |> 
    left_join(training_dates) |> 
    filter(notification_date >= cutoff_date) |> 
    filter(location == x_location, pathogen == x_pathogen) |> 
    
    filter(
      # Drop the last day of data (often issues)
      # If this improves, remove this line
      notification_date < max(notification_date)
    )
  
  training_counts <- reduced_counts |> 
    make_training_counts(x_location, x_pathogen, n_days_ahead)
  
  return(training_counts)
}


get_gi <- function(pathogen) {
  # From ACEFA standards
  gis <- list(
    "SARSCOV2" = make_gi(0.89, rate = 0.27, 0.995),
    "RSV" = make_gi(12.76, 1.7, 0.995),
    "flu" = make_gi(4.00, 1.54, 0.995)
  )
  
  if(!(pathogen %in% names(gis))) {
    stop("Invalid pathogen for gi")
  }
  
  gi <- gis[[pathogen]]
}


make_training_counts <- function(
    reduced_counts,
    x_location, x_pathogen, 
    days_ahead = 0
) {
  school_holidays_location <- get_school_holidays_long() |> 
    filter(location == x_location)
  
  public_holidays_location <- get_public_holidays() |> 
    filter(location == x_location) |> 
    mutate(public_holiday_type = match(holiday_name, unique(holiday_name))) |> 
    select(date, location, public_holiday_type)
  
  training_counts <- reduced_counts |> 
    select(location, pathogen, date = notification_date, y = cases) |> 
    complete(date = seq(min(date), max(date), by = "days"), fill = list(location = x_location, pathogen = x_pathogen)) |>
    arrange(date)
  
  n_na <- sum(is.na(training_counts$y))
  
  print(str_c(n_na, " NA values to be imputed"))
  
  if(n_na > 0) {
    training_counts |> 
      filter(is.na(y)) |> 
      print()
    
    training_counts <- training_counts |> 
      fill(y, .direction = "down")
  }
  
  if(days_ahead > 0) {
    training_counts <- bind_rows(
      training_counts,
      tibble(
        date = seq(max(training_counts$date) + days(1), max(training_counts$date) + days(days_ahead), by = "days"),
        
        location = x_location, pathogen = x_pathogen
      )
    )
  }
  
  
  training_counts_full <- training_counts |>
    arrange(date) |> 
    mutate(
      x = as.numeric(date - min(date)),
      day_of_week = wday(date, week_start = 7),
      is_weekend = if_else(day_of_week == 1 | day_of_week == 7, 1, 0),
      year = as.integer(
        interval(min(floor_date(date, "year")), floor_date(date, "year")) / years(1)
      ) + 1
    ) |> 
    left_join(school_holidays_location, by = join_by(date, location)) |> 
    left_join(public_holidays_location, by = join_by(date, location)) |> 
    mutate(school_holiday_type = replace_na(school_holiday_type, 0) + 1,
           public_holiday_type = replace_na(public_holiday_type, 0) + 1)
  
  # Make sure we have a complete sequence of dates
  n_diff <- diff(training_counts_full$x) |> table() |> unique() |> length()
  
  stopifnot(n_diff == 1)
  
  return(training_counts_full)
}


make_gi <- function(shape, rate, p_threshold) {
  cdf_gi <- pgamma(0:21, shape = shape, rate = rate)
  
  discrete_gi <- diff(cdf_gi)
  discrete_gi[length(discrete_gi)] <- discrete_gi[length(discrete_gi)] + 1 - sum(discrete_gi)
  
  return(discrete_gi)
}



find_inv_gamma_params <- function(l, u, tail_prob = 0.01) {
  f <- function(x) {
    alpha <- exp(x[1])
    theta <- exp(x[2])
    
    cdf_l <- pgamma(1/l, shape = alpha, rate = theta, lower.tail = FALSE)
    cdf_u <- pgamma(1/u, shape = alpha, rate = theta, lower.tail = FALSE)
    
    y1 <- cdf_l - tail_prob
    y2 <- cdf_u - (1.0 - tail_prob)
    
    return(c(y1, y2))
  }
  
  initial_guess <- c(log(5.0), log(5.0))
  
  require(nleqslv)
  result <- nleqslv(initial_guess, f)
  
  if (result$termcd %in% c(1, 2)) {
    sol <- exp(result$x)
    print(str_c("inv_gamma(", round(sol[1], digits = 3), ", ", round(sol[2], digits = 3), ")"))
    return(list(alpha = sol[1], theta = sol[2]))
  } else {
    stop(sprintf("Solver failed to converge. Return code: %d, Message: %s", 
                 result$termcd, result$message))
  }
}


create_stan_data <- function(
    training_counts,
    days_ahead,
    x_pathogen,
    x_location
) {
  gi <- get_gi(x_pathogen)
  n_days_data <- nrow(training_counts) - n_days_ahead
  
  N_public_holiday_types <- max(training_counts$public_holiday_type) - 1
  N_school_holiday_types <- max(training_counts$school_holiday_type) - 1
  N_years <- max(training_counts$year)
  
  year_coverage <- training_counts |> count(year) |> mutate(year_covered = pmin(365, n) / 365) |> arrange(year) |> pull(year_covered)
  
  year_ranges <- training_counts |> group_by(year) |> summarise(ix_year_start = min(x) + 1, ix_year_end = max(x) + 1)
  
  data_list <- list(
    x = training_counts$x,
    day_of_week = training_counts$day_of_week,
    is_weekend = training_counts$is_weekend,
    school_holiday_type = training_counts$school_holiday_type,
    
    public_holiday_type = training_counts$public_holiday_type,
    
    y = training_counts$y[1:n_days_data],
    N = nrow(training_counts), N_data = n_days_data,
    N_public_holiday_types = N_public_holiday_types,
    N_school_holiday_types = N_school_holiday_types,
    
    M_annual = 40,
    M_noise = calculate_M_matern(c = 1.25, length_scale = 10.0, x = training_counts$x),
    c_factor = 1.25,
    
    gi_reversed = rev(gi), N_gi = length(gi)
  )
  
  return(data_list)
}