
fit_model <- function(data_list, model_file = "src/gp_renewal.stan") {
  
  print(str_c("Fitting model ", model_file))
  
  library(cmdstanr)
  model <- cmdstan_model(model_file)
  
  fit <- model$sample(
    data = data_list,
    seed = 1,
    chains = 4,
    parallel_chains = 4,
    refresh = 50,
    iter_warmup = 3000, iter_sampling = 1000,
    init = 0.5
  )
  
  return(fit)
}


get_origin_date <- function(combined_counts, x_pathogen, x_location, x_target) {
  combined_counts |> 
    filter(location == x_location, pathogen == x_pathogen, target == x_target,
           is_complete) |> 
    pull(notification_date) |> 
    max()
}

write_forecast_file <- function(round_id, x_pathogen, x_location, x_target) {
  round_results <- fs::path("results", round_id)
  
  fit <- read_rds(str_c(round_results, "/", x_pathogen, "_", x_location, "_", x_target, ".rds"))
  training_counts <- read_parquet(str_c(round_results, "/", x_pathogen, "_", x_location, "_", x_target, "_training_counts.parquet"))
  
  origin_date <- training_counts$origin_date[[1]]
  
  case_trajectories <- fit |> 
    spread_draws(cases[ix_x], latent_cases_log[ix_x]) |> 
    mutate(date = training_counts$date[ix_x],
           latent_cases = exp(latent_cases_log))
  
  peak_forecasts <- case_trajectories |> 
    filter(year(date) == year(round_id)) |> 
    group_by(.draw) |> 
    summarise(peak_case_incidence = max(latent_cases),
              peak_date = first(date[latent_cases == max(latent_cases)])) |> 
    ungroup() |> 
    transmute(
      output_type = "sample",
      output_type_id = .draw,
      peak_case_incidence, 
      peak_day_of_year = yday(peak_date),
      location = x_location,
      pathogen = x_pathogen,
      round_id = ymd(round_id),
      horizon = 0,
      origin_date = origin_date
    ) |> 
    pivot_longer(c(peak_case_incidence, peak_day_of_year),
                 names_to = "target", values_to = "value") |> 
    
    transmute(
      round_id, origin_date,
      target, horizon, location, pathogen,
      
      output_type, output_type_id,
      
      value = as.integer(round(value))
    )
  
  
  case_forecasts <- case_trajectories |>
    ungroup() |>
    transmute(
      round_id = ymd(round_id),
      origin_date = origin_date,
      target = x_target,
      horizon = as.integer(date - origin_date),
      location = x_location,
      pathogen = x_pathogen,
      output_type = "sample",
      output_type_id = .draw,
      value = as.integer(replace_na(cases, 0))
    ) |>
    filter(horizon >= -7, horizon <= 7 * 8)
  
  case_forecasts_smoothed <- case_trajectories |>
    add_forecast_noise(training_counts) |> 
    ungroup() |>
    transmute(
      round_id = ymd(round_id),
      origin_date = origin_date,
      target = "case_incidence_smoothed",
      horizon = as.integer(date - origin_date),
      location = x_location,
      pathogen = x_pathogen,
      output_type = "sample",
      output_type_id = .draw,
      value = as.integer(replace_na(count, 0))
    ) |>
    filter(horizon >= -7, horizon <= 7 * 8)
  
  forecast_data <- bind_rows(
    case_forecasts,
    
    # Only include derived targets if target is case_incidence
    case_forecasts_smoothed |> filter(x_target == "case_incidence"),
    peak_forecasts |> filter(x_target == "case_incidence")
  )
  
  output_file <- str_c(round_results, "/", x_pathogen, "_", x_location, "_forecast-trajectories.parquet")
  
  arrow::write_parquet(forecast_data, output_file)
}




process_counts <- function(
    combined_counts,
    x_location, x_pathogen, x_target,
    n_days_ahead,
    origin_date
) {
  training_dates <- get_training_dates()
  
  # Filter the counts for our given location and pathogen
  # and drop data prior to the training date
  reduced_counts <- combined_counts |> 
    left_join(training_dates) |>
    
    filter(notification_date >= cutoff_date) |>
    filter(location == x_location, pathogen == x_pathogen, target == x_target) |> 
    filter(is_complete) |> 
    
    filter(
      # Drop the last day of data (often issues)
      # If this improves, remove this line
      notification_date < max(notification_date)
    )
  
  training_counts <- reduced_counts |> 
    make_training_counts(x_location, x_pathogen, x_target, n_days_ahead) |> 
    # Save some extra data in case we need it later
    mutate(origin_date = origin_date)
  
  return(training_counts)
}


get_gi <- function(pathogen) {
  # From ACEFA standards
  gis <- list(
    "SARSCOV2" = make_gi(0.89, rate = 0.27),
    "RSV" = make_gi(12.76, 1.7),
    "flu" = make_gi(4.00, 1.54)
  )
  
  if(!(pathogen %in% names(gis))) {
    stop("Invalid pathogen for gi")
  }
  
  gi <- gis[[pathogen]]
  
  return(gi)
}


make_training_counts <- function(
    reduced_counts,
    x_location, x_pathogen, x_target,
    days_ahead = 0,
    holidays_cutoff = ymd("2023-01-01")
) {
  school_holidays_location <- get_school_holidays_long(date_cutoff = holidays_cutoff) |> 
    filter(location == x_location)
  
  public_holidays_location <- get_public_holidays(date_cutoff = holidays_cutoff) |> 
    filter(location == x_location) |> 
    mutate(public_holiday_type = match(holiday_name, unique(holiday_name))) |> 
    select(date, location, public_holiday_type)
  
  training_counts <- reduced_counts |> 
    select(location, pathogen, target, date = notification_date, 
           y = value, y_trend = value_trend, y_trend_and_remainder = value_trend_and_remainder) |> 
    complete(date = seq(min(date), max(date), by = "days"), fill = list(location = x_location, pathogen = x_pathogen, target = x_target)) |>
    arrange(date)
  
  n_na <- sum(is.na(training_counts$y))
  
  message(str_c(n_na, " NA values to be imputed"))
  
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
        
        location = x_location, pathogen = x_pathogen, target = x_target
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
      ) + 1,
      day_of_year = yday(date)
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


make_gi <- function(shape, rate) {
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
  
  library(nleqslv)
  result <- nleqslv(initial_guess, f)
  
  if (result$termcd %in% c(1, 2)) {
    sol <- exp(result$x)
    message(glue::glue("inv_gamma({round(sol[1], 3)}, {round(sol[2], 3)})"))
    return(list(alpha = sol[1], theta = sol[2]))
  } else {
    stop(sprintf("Solver failed to converge. Return code: %d, Message: %s", 
                 result$termcd, result$message))
  }
}


calculate_M_matern <- function(c, length_scale, x) {
  x_centered <- x - mean(x) 
  S <- max(abs(x_centered))
  
  # Formula of Riutort-Mayol 2023
  M <- 2.65 * (c * S) / length_scale
  
  M <- max(10, M)
  
  return(round(M))
}

get_population_size <- function(location) {
  population_sizes <- list(
    # December 2021 ERP - ABS
    "ACT" = 453.3 * 1000,
    "NSW" = 8095.4 * 1000,
    "NT" = 249.3 * 1000,
    "QLD" = 5265.0 * 1000,
    "SA" = 1806.6 * 1000,
    "TAS" = 569.8 * 1000,
    "VIC" = 6559.9 * 1000,
    "WA" = 2762.2 * 1000,
    
    # September 25 ERP - Stats NZ
    "NZ" = 5334.2 * 1000
  )
  
  if(!(location %in% names(population_sizes))) {
    stop("Invalid location for population size")
  }
  
  return(population_sizes[[location]])
}



create_stan_data <- function(
    training_counts,
    n_days_ahead,
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
  
  humidity_data <- read_csv("data/modelled_humidity.csv", show_col_types = FALSE) |> 
    filter(location == x_location) |> 
    arrange(day_of_year)
  
  data_list <- list(
    x = training_counts$x,
    day_of_week = training_counts$day_of_week,
    day_of_year = training_counts$day_of_year,
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
    
    population_size = get_population_size(x_location),
    
    gi_reversed = rev(gi), N_gi = length(gi),
    
    # Immunity model params
    k = 5, init_factor = 0.5e-2,
    humidity = humidity_data$humidity_scaled,
    sample_from_prior = 0, log_R_0 = 0.33
  )
  
  return(data_list)
}