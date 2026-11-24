# Minimal example: fit the baseline model and produce a 28-day forecast for a
# single location-pathogen-target, using the anonymised counts in this repository.
#
#   Rscript R/example_forecast.R
#
# Requires CmdStan (see README). Fitting takes roughly an hour on 4 cores.

source("R/functions.R")

x_location <- "ACT"
x_pathogen <- "RSV"
x_target <- "case_incidence"
origin_date <- ymd("2025-05-01")

n_days_ahead <- 365 # as in the study; only the first 28 days are shown below
n_days_forecast <- 28

all_counts <- read_csv("data/combined-count-smoothed_final_anonymised.csv", show_col_types = FALSE)

# Pretend we are forecasting at origin_date: keep only data available by then
combined_counts <- all_counts |> 
  filter(notification_date <= origin_date)

training_counts <- process_counts(combined_counts, x_location, x_pathogen, x_target, n_days_ahead, origin_date)

data_list <- create_stan_data(training_counts, n_days_ahead, x_pathogen, x_location)

fit <- fit_model(
  data_list, model_file = "src/gp_renewal.stan",
  n_warmup = 1000, n_samples = 1000
)

# Forecast quantiles of daily case counts
forecast_quants <- fit |>
  spread_draws(cases[ix_t]) |>
  mutate(date = training_counts$date[ix_t]) |>
  filter(date > origin_date - days(56), date <= origin_date + days(n_days_forecast)) |>
  group_by(date) |>
  median_qi(cases = replace_na(cases, 0), .width = c(0.5, 0.8, 0.95))

# Observed counts, including those after origin_date for comparison
observed <- all_counts |>
  filter(location == x_location, pathogen == x_pathogen, target == x_target,
         notification_date > origin_date - days(56),
         notification_date <= origin_date + days(n_days_forecast)) |>
  mutate(after_origin = notification_date > origin_date)

ggplot() +
  geom_ribbon(aes(x = date, ymin = .lower, ymax = .upper, group = .width),
              forecast_quants, fill = "steelblue", alpha = 0.25) +
  geom_line(aes(x = date, y = cases), filter(forecast_quants, .width == 0.5), colour = "steelblue4") +
  geom_point(aes(x = notification_date, y = value, colour = after_origin), observed, size = 0.8) +
  geom_vline(xintercept = origin_date, linetype = "dashed") +
  scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = "grey60"), guide = "none") +
  labs(x = NULL, y = "Daily cases",
       title = str_c(x_location, " ", x_pathogen, " forecast from ", origin_date)) +
  theme_minimal()

