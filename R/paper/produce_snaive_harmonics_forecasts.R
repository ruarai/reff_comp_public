

source("R/functions.R")

snaive_harmonics_forecast <- function(cases, horizon_days = 7 * 8) {
  library(tsibble)
  library(fabletools)
  library(distributional)
  
  cases |>
    left_join(get_training_dates()) |> 
    filter(notification_date >= cutoff_date) |> 
    rename(date = notification_date) |>
    mutate(log_counts = log(value + 1)) |> 
    as_tsibble(
      key = c(pathogen, location, target),
      index = date
    ) |>
    model(ARIMA(
      log_counts ~ PDQ(0, 0, 0) + pdq(d = 1) +
        fourier(period = "week", K = 3) +
        fourier(period = "year", K = 10)
    )) |>
    forecast(h = horizon_days) |>
    as_tibble() |>
    group_by(pathogen, location, target) |>
    mutate(
      origin_date = min(date) - 1,
      value = dist_truncated(exp(log_counts) - 1, lower = 0),
      value_log = log(dist_truncated(value, lower = 0) + 1),
    ) |>
    ungroup()
}


calculate_crps <- function(dist, truth) {
  dist_family <- family(dist)
  scores <- rep(NA_real_, length(dist_family))
  
  probs <- seq(0.01, 0.99, by = 0.01)
  x <- quantile(dist, p = probs)
  percentile <- list(probs * 100)
  
  scores <- purrr::pmap_dbl(
    list(
      x,
      percentile,
      truth
    ),
    function(x, p, y) {
      2 * mean(
        ifelse(y < x, (1 - p / 100), p / 100) * abs(x - y),
        na.rm = TRUE
      )
    },
    .progress = TRUE
  )
  
  scores
}


library(fable)
combined_count_latest <- read_csv(
  "data/combined-count-smoothed_final_anonymised.csv",
  show_col_types = FALSE
) |> 
  filter(target == "case_incidence", location == "SA", pathogen == "RSV")

test_forecast <- snaive_harmonics_forecast(combined_count_latest)

ggplot() +
  stat_lineribbon(aes(x = date, ydist = value),
                  linewidth = 0,
                  test_forecast) +
  
  geom_line(aes(x = notification_date, y = value),
            filter(combined_count_latest, location == "SA", pathogen == "RSV", notification_date >= ymd("2023-01-01"))) +
  
  
  coord_cartesian(ylim = c(0, 150)) +
  
  theme_bw() +
  
  scale_fill_brewer()

library(furrr)
library(future.callr)

plan(callr, workers = 4)

forecast_data <- tibble(
  round_id = c(seq(ymd("2025-04-10"), ymd("2025-12-25"), by = "1 week"),
               seq(ymd("2025-01-17"), ymd("2025-04-10"), by = "1 week"))
) |> 
  mutate(case_file = glue::glue("data/retro_rounds/{round_id}/combined-count-{round_id}.csv")) |> 
  
  future_pmap(function(round_id, case_file) {
    print(round_id)
    if(!fs::file_exists(case_file)) {
      print(glue::glue("No data found for {round_id}"))
      return(NULL)
    }
    
    cases <- read_csv(case_file, show_col_types = FALSE) |> 
      filter(target == "case_incidence", is_complete)
    
    case_forecast <- snaive_harmonics_forecast(cases) |> 
      mutate(round_id = round_id, .before = 1)
  }) |> 
  bind_rows() |> 
  left_join(combined_count_latest |> rename(date = notification_date, true_value = value))


forecast_data$crps <- calculate_crps(forecast_data$value, forecast_data$true_value)
forecast_data$crps_log <- calculate_crps(forecast_data$value_log, log(forecast_data$true_value + 1))



snaive_perf <- forecast_data |> 
  mutate(horizon = as.numeric(date - origin_date)) |> 
  select(round_id, pathogen, location, target, date, horizon, origin_date,
         crps, crps_log) |> 
  mutate(.model = "snaive_harmonics",
         round_id = as.character(round_id))

snaive_perf |> 
  write_parquet("data/snaive_harmonics_perf_extended.parquet")





