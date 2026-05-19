

get_peak_priors <- function() {
  get_normal_params <- function(q, prob = 0.99) {
    p1 <- (1 - prob) / 2
    p2 <- 1 - (1 - prob) / 2
    
    z1 <- qnorm(p1)
    z2 <- qnorm(p2)
    
    sd_val <- (q[2] - q[1]) / (z2 - z1)
    mean_val <- q[1] - (sd_val * z1)
    
    return(c(mean = mean_val, sd = sd_val))
  }
  
  year_start <- ymd("2026-01-01")
  
  peak_priors <- tribble(
    ~location, ~pathogen, ~mean_peak_date, ~peak_yday_sd, ~peak_cases_q005, ~peak_cases_q995,
    "ACT", "RSV", ymd("2026-07-01"), 45,  10,  40,
    "NSW", "RSV", ymd("2026-06-01"), 30, 180, 600,
    "QLD", "RSV", ymd("2026-07-01"), 90,  75, 500,
    "SA",  "RSV", ymd("2026-07-01"), 60,   60,  200,
    "TAS", "RSV", ymd("2026-07-01"), 60,   5,  40,
    "VIC", "RSV", ymd("2026-07-01"), 30, 125, 400,
    "WA",  "RSV", ymd("2026-08-01"), 60, 40, 180,
    
    # Assume flu might be a bit lower this year due to big year in 2025
    "ACT", "flu", ymd("2026-07-01"), 45, 10,   150,
    "NSW", "flu", ymd("2026-07-01"), 30, 250, 5000,
    "QLD", "flu", ymd("2026-07-01"), 90, 200, 1500,
    "SA",  "flu", ymd("2026-07-01"), 45, 50,    400,
    "TAS", "flu", ymd("2026-07-01"), 45, 10,    60,
    "VIC", "flu", ymd("2026-07-01"), 45, 50, 1200,
    "WA",  "flu", ymd("2026-08-01"), 45, 50,  600,
    
    # More uncertain timing with SARS-CoV-2 (unclear epi)
    "ACT", "SARSCOV2", ymd("2026-07-01"), 90, 10,   45,
    "NSW", "SARSCOV2", ymd("2026-07-01"), 90, 250, 1000,
    "QLD", "SARSCOV2", ymd("2026-07-01"), 90, 100, 450,
    "SA",  "SARSCOV2", ymd("2026-07-01"), 90, 30,    120,
    "TAS", "SARSCOV2", ymd("2026-07-01"), 90, 5,    40,
    "VIC", "SARSCOV2", ymd("2026-07-01"), 90, 50, 400,
    "WA",  "SARSCOV2", ymd("2026-07-01"), 90, 25,  150,
  ) |> 
    
    rowwise() |> 
    mutate(peak_size_normal = list(get_normal_params(log(c(peak_cases_q005, peak_cases_q995)))),
           log_peak_size_mean = peak_size_normal[1], log_peak_size_sd = peak_size_normal[2],
           peak_yday_mean = as.numeric(mean_peak_date - year_start) + 1) |> 
    select(-peak_size_normal) |> 
    ungroup()
}

plot_peak_priors <- function() {
  
  round_id <- "2026-04-29"
  
  year <- as.character(year(round_id))
  combined_counts <- read_csv(str_c("../forecast-hub/processed-data/", year, "/", round_id, "/combined-count-smoothed-", round_id, ".csv"))
  
  
  plot_data_counts <- combined_counts |> 
    filter(notification_date >= ymd("2023-01-01"),
           target == "case_incidence",
           location %in% c("VIC", "ACT", "WA", "NSW", "TAS", "QLD", "SA"))
  
  n_samples <- 4000
  peak_samples <- get_peak_priors() |> 
    rowwise() |> 
    mutate(
      day_of_year = list(rnorm(n_samples, peak_yday_mean, peak_yday_sd)),
      peak_size = list(exp(rnorm(n_samples, log_peak_size_mean, log_peak_size_sd)))
    ) |> 
    unnest(c(day_of_year, peak_size))
  
  
  
  ggplot(plot_data_counts) +
    
    geom_line(aes(x = notification_date, y = value_trend)) +
    
    ggh4x::facet_grid2(rows = vars(location), cols = vars(pathogen),
                       scales = "free_y", independent = "y") +
    
    geom_density2d(aes(x = ymd("2026-01-01") + day_of_year - 1, y = peak_size),
                   contour_var = "ndensity", alpha = 0.6,
                   peak_samples, bins = 6) +
    
    geom_point(aes(x = ymd("2026-01-01") + day_of_year - 1, y = peak_size),
               size = 0.4, stroke = 0.4, alpha = 0.2,
                   peak_samples |> sample_n(1000)) +
    
    coord_cartesian(xlim = c(ymd("2023-01-01"), ymd("2027-01-01"))) +
    
    theme_bw()
}
