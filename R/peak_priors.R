

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
    "VIC", "RSV", ymd("2026-07-01"), 30, 150, 400,
    "WA",  "RSV", ymd("2026-08-01"), 60, 60, 180,
    
    "ACT", "flu", ymd("2026-07-01"), 45, 20,   150,
    "NSW", "flu", ymd("2026-07-01"), 30, 500, 5000,
    "QLD", "flu", ymd("2026-07-01"), 90, 300, 1500,
    "SA",  "flu", ymd("2026-07-01"), 45, 80,    400,
    "TAS", "flu", ymd("2026-07-01"), 45, 15,    60,
    "VIC", "flu", ymd("2026-07-01"), 45, 300, 1200,
    "WA",  "flu", ymd("2026-08-01"), 45, 100,  600,
    
    "ACT", "SARSCOV2", ymd("2026-07-01"), 45, 10,   45,
    "NSW", "SARSCOV2", ymd("2026-07-01"), 45, 250, 1000,
    "QLD", "SARSCOV2", ymd("2026-07-01"), 45, 100, 450,
    "SA",  "SARSCOV2", ymd("2026-07-01"), 45, 30,    120,
    "TAS", "SARSCOV2", ymd("2026-07-01"), 45, 5,    40,
    "VIC", "SARSCOV2", ymd("2026-07-01"), 45, 50, 400,
    "WA",  "SARSCOV2", ymd("2026-07-01"), 45, 25,  150,
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
  
  
  peak_priors <- get_peak_priors()
  
  n_samples <- 4000
  peak_samples <- peak_priors |> 
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
                   contour_var = "ndensity",
                   peak_samples, bins = 4) +
    
    coord_cartesian(xlim = c(ymd("2023-01-01"), ymd("2027-01-01"))) +
    
    theme_bw()
}
