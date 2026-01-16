

source("R/holidays.R")
source("R/plotting_functions.R")


get_training_dates <- function() {
  locations <- c("ACT", "NSW", "NT", "QLD", "TAS", "VIC", "WA", "SA", "NZ")
  training_dates <- bind_rows(
    expand_grid(location = locations, pathogen = "flu", cutoff_date = ymd("2023-01-01")),
    expand_grid(
      location = c("ACT", "NSW", "NT", "QLD", "TAS", "VIC", "WA"),
      pathogen = "SARSCOV2", cutoff_date = ymd("2023-03-01")
    ),
    expand_grid(
      location = c("NZ", "SA"),
      pathogen = "SARSCOV2", cutoff_date = ymd("2024-09-01")
    ),
    expand_grid(location = locations, pathogen = "RSV", cutoff_date = ymd("2023-01-01")),
  )
  
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
    select(location, pathogen, date = notification_date, y = cases, y_smooth = cases_trend_and_remainder) |> 
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
  
  
  training_counts <- training_counts |>
    arrange(date) |> 
    mutate(x = as.numeric(date - min(date)),
           day_of_week = wday(date)) |> 
    left_join(school_holidays_location, by = join_by(date, location)) |> 
    left_join(public_holidays_location, by = join_by(date, location)) |> 
    mutate(school_holiday_type = replace_na(school_holiday_type, 0) + 1,
           public_holiday_type = replace_na(public_holiday_type, 0) + 1)
  
  
  n_diff <- diff(training_counts$x) |> table() |> unique() |> length()
  
  stopifnot(n_diff == 1)
  
  return(training_counts)
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
    return(list(alpha = sol[1], theta = sol[2]))
  } else {
    stop(sprintf("Solver failed to converge. Return code: %d, Message: %s", 
                 result$termcd, result$message))
  }
}



make_intervals <- function(
    tbl,
    vars,
    grouping_var = rlang::sym(".draw"),
    date_var = rlang::sym("date"),
    probs = c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9),
    na.rm = FALSE
) {
  tbl |> 
    ungroup() |> 
    select(!!date_var, !!grouping_var, all_of(vars)) |> 
    pivot_longer(
      cols = c(-!!grouping_var, -!!date_var), 
      names_to = "variable", 
      values_to = "value"
    ) |> 
    pivot_wider(
      names_from = !!grouping_var, 
      names_prefix = "sample_", 
      values_from = value
    ) |> 
    make_intervals_internal(probs = probs, na.rm = na.rm)
}

make_intervals_internal <- function(
    tbl, 
    probs = c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), 
    prefix = "sample_", 
    na.rm = FALSE
) {
  data_matrix <- tbl |> 
    select(starts_with(prefix)) |> 
    as.matrix()
  
  id_tbl <- tbl |> 
    select(!starts_with(prefix))
  
  medians <- data_matrix |> 
    matrixStats::rowMedians(na.rm = na.rm)
  
  quant_probs <- c(rev(1 - probs) / 2, 0.5 + probs / 2)
  quant_labels <- as.character(probs * 100)
  
  quant_names <- c(
    str_c("lower_", rev(quant_labels)), 
    str_c("upper_", quant_labels)
  )
  
  quants <- data_matrix |> 
    matrixStats::rowQuantiles(probs = quant_probs, na.rm = na.rm)
  
  colnames(quants) <- quant_names
  
  bind_cols(id_tbl, median = medians, as_tibble(quants)) |> 
    pivot_longer(
      cols = all_of(quant_names),
      names_to = c(".value", "quant"),
      names_pattern = "(lower|upper)_(.*)"
    ) |> 
    mutate(
      quant = factor(quant, levels = quant_labels) |> fct_rev()
    )
}


get_gi <- function(pathogen) {
  gis <- list(
    "SARSCOV2" = make_gi(0.89, rate = 0.27, 0.995),
    "RSV" = make_gi(12.76, 1.7, 0.995),
    "flu" = make_gi(4.00, 1.54, 0.995)
  )
  
  gi <- gis[[pathogen]]
}

make_retrospective_jobs <- function() {
  hub_path <- "../forecast-hub/"
  
  rounds <- fs::dir_ls(fs::path(hub_path, "processed-data", "2025"))
  
  tibble(
    round_dir = rounds
  ) |> 
    mutate(round_id = str_extract(round_dir, "\\d{4}-\\d{2}-\\d{2}"),
           counts_file = fs::path(round_dir, str_c("combined-count-", round_id, ".csv")),
           smoothed_file = fs::path(round_dir, str_c("combined-count-smoothed-", round_id, ".csv"))) |> 
    
    filter(fs::file_exists(counts_file), fs::file_exists(smoothed_file)) |> 
    
    rowwise() |> 
    mutate(counts = list(read_csv(counts_file, show_col_types = FALSE)), 
           jobs = list(distinct(counts, location, pathogen))) |> 
    select(-counts) |> 
    unnest(jobs) |> 
    
    select(round_id, location, pathogen)
}

calculate_M_matern <- function(c, length_scale, x) {
  x_centered <- x - mean(x) 
  S <- max(abs(x_centered))
  
  # Formula of Riutort-Mayol 2023
  M <- 2.65 * (c * S) / length_scale
  
  M <- max(10, M)
  
  return(round(M))
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
  training_counts <- arrow::read_parquet(str_c(round_results, "/", x_pathogen, "_", x_location, "_training_counts.parquet"))
  
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
      value = as.integer(value)
    ) |> 
    
    filter(horizon >= -7, horizon <= 28) |> 
    
    select(
      round_id, origin_date,
      target, horizon, location, pathogen,
      
      output_type, output_type_id,
      
      value
    )
  
  output_file <- fs::path("forecast_outputs_raw", str_c(x_location, "-", round_id, "-uom-gp-reff.parquet"))
  
  arrow::write_parquet(forecast_data, output_file)
}




process_counts <- function(
    combined_counts, combined_counts_smoothed,
    x_location, x_pathogen,
    n_days_ahead
) {
  training_dates <- get_training_dates()
  
  all_counts <- combined_counts |> 
    left_join(combined_counts_smoothed)
  
  combined_counts_outliers <- all_counts |> 
    left_join(training_dates) |> 
    filter(notification_date >= cutoff_date) |> 
    filter(location == x_location, pathogen == x_pathogen) |> 
    mutate(outlier = abs(log(cases_trend_and_remainder + 0.5) - log(cases_trend + 0.5)) > 1.38 & cases_trend > 20) |> 
    
    filter(
      notification_date < max(notification_date)
    )
  
  print(str_c("Dropping ", sum(combined_counts_outliers$outlier), " outliers"))
  
  training_counts <- combined_counts_outliers |> 
    filter(!outlier) |> 
    make_training_counts(x_location, x_pathogen, n_days_ahead)
  
  return(training_counts)
}

create_stan_data <- function(
    training_counts,
    days_ahead,
    x_pathogen
) {
  gi <- get_gi(x_pathogen)
  n_days_data <- nrow(training_counts) - n_days_ahead
  
  N_public_holiday_types <- length(unique(training_counts$public_holiday_type))
  
  data_list <- list(
    x = training_counts$x,
    day_of_week = training_counts$day_of_week,
    school_holiday_type = training_counts$school_holiday_type,
    
    public_holiday_type = training_counts$public_holiday_type,
    
    y = training_counts$y[1:n_days_data],
    N = nrow(training_counts), N_data = n_days_data,
    N_public_holiday_types = N_public_holiday_types,
    
    M_annual = 40,
    M_noise = calculate_M_matern(c = 1.25, length_scale = 10.0, x = training_counts$x),
    c_factor = 1.25,
    
    gi_reversed = rev(gi), N_gi = length(gi)
  )
  
  return(data_list)
}

fit_model <- function(data_list) {
  require(cmdstanr)
  model <- cmdstan_model("src/gp_renewal.stan")
  
  model$sample(
    data = data_list,
    seed = 2,
    chains = 4,
    parallel_chains = 4,
    refresh = 50,
    iter_warmup = 1000, iter_sampling = 1000
  )
}
