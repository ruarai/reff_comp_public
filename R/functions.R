

source("R/school_holidays.R")
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
  
  training_counts <- reduced_counts |> 
    select(location, pathogen, date = notification_date, y = cases, y_smooth = cases_trend) |> 
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
    left_join(school_holidays_location |> mutate(school_holidays = 1), by = join_by(date)) |> 
    mutate(school_holidays = replace_na(school_holidays, 0))
  
  
  n_diff <- diff(training_counts$x) |> table() |> unique() |> length()
  
  stopifnot(n_diff == 1)
  
  return(training_counts)
}


make_gi <- function(shape, rate, p_threshold) {
  cdf_gi <- pgamma(0:30, shape = shape, rate = rate)
  
  discrete_gi <- diff(cdf_gi[cdf_gi < p_threshold])
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
    probs = c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9),
    prefix = "sample_",
    na.rm = FALSE
) {
  data_matrix <- tbl %>%
    select(starts_with(prefix)) %>%
    as.matrix()
  
  id_tbl <- tbl %>%
    select(!starts_with(prefix))
  
  medians <- data_matrix %>%
    matrixStats::rowMedians() %>%
    tibble(median = .)
  
  quant_probs <- c(rev(1 - probs) / 2, 0.5 + probs / 2)
  quant_names <- c(str_c("lower_", rev(probs) * 100), str_c("upper_", probs * 100))
  
  quants <- data_matrix %>%
    matrixStats::rowQuantiles(probs = quant_probs, na.rm = na.rm) %>%
    `colnames<-`(quant_names) %>%
    as_tibble() %>%
    bind_cols(id_tbl, medians, .) %>%
    pivot_longer(cols = -all_of(c(colnames(id_tbl), colnames(medians))),
                 names_to = c("type", "quant"),
                 names_sep = "_") %>%
    pivot_wider(names_from = "type",
                values_from = "value") %>%
    
    mutate(quant = factor(quant, levels = as.character(probs * 100)) %>% fct_rev())
  
  quants
}

get_gi <- function(pathogen) {
  gis <- list(
    "SARSCOV2" = make_gi(0.89, rate = 0.27, 0.995),
    "RSV" = make_gi(12.76, 1.7, 0.995),
    "flu" = make_gi(4.00, 1.54, 0.995)
  )
  
  gi <- gis[[pathogen]]
}
