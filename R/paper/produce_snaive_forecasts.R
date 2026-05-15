
# Adapted from https://github.com/robmoss/paper-2025-winter-forecasting/blob/master/R/snaive.R
# itself adapted from code by Mitch O'Hara-Wild
snaive_forecast <- function(cases, start_date) {
  library(tsibble)
  library(fabletools)
  library(distributional)
  
  cases |>
    rename(date = notification_date) |>
    filter(date > !! start_date) |>
    as_tsibble(
      key = c(pathogen, location),
      index = date
    ) |>
    model(fable::SNAIVE(value)) |>
    forecast(h = 7 * 8) |>
    as_tibble() |>
    group_by(pathogen, location) |>
    mutate(
      origin_date = min(date) - 1,
      value = dist_truncated(value, lower = 0),
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

start_date <- ymd("2024-01-01")


combined_count_latest <- read_csv(
  "../forecast-hub/processed-data/2026/2026-05-06/combined-count-2026-05-06.csv",
  show_col_types = FALSE
) |> 
  filter(target == "case_incidence")

test_forecast <- snaive_forecast(combined_count_latest, start_date)

ggplot() +
  stat_lineribbon(aes(x = date, ydist = value_log),
                  test_forecast |> filter(location == "NSW", pathogen == "RSV"))


forecast_data <- tibble(
  round_id = c(seq(ymd("2025-04-10"), ymd("2025-12-25"), by = "1 week"),
               seq(ymd("2025-01-17"), ymd("2025-04-10"), by = "1 week"))
) |> 
  mutate(case_file = glue::glue("../forecast-hub/processed-data/2025/{round_id}/combined-count-{round_id}.csv")) |> 
  
  pmap(function(round_id, case_file) {
    if(!fs::file_exists(case_file)) {
      print(glue::glue("No data found for {round_id}"))
      return(NULL)
    }
    
    cases <- read_csv(case_file, show_col_types = FALSE) |> 
      filter(target == "case_incidence", is_complete)
    
    case_forecast <- snaive_forecast(cases, start_date) |> 
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
  mutate(.model = "snaive",
         round_id = as.character(round_id))

snaive_perf |> 
  write_parquet("data/snaive_perf_extended.parquet")
