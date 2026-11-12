# Summarise forecast scores for public release.
#
# The per-day scores (saved_results/<model_version>/**/*_performance.parquet) were computed
# against the restricted case counts, and per-day point-error scores can be inverted to
# recover those counts. This script reduces them to the level the manuscript figures need:
# mean CRPS (log scale) per model, pathogen, jurisdiction, round and horizon week.
#
# Requires the per-day scores for every model version and the random-walk baseline scores
# (data/snaive_perf_extended.parquet, from produce_snaive_forecasts.R).
#
# Writes saved_results/perf_summary.parquet.

library(tidyverse)
library(arrow)

source("R/paper/common.R")

model_versions <- c(
  "baseline", "no_holidays_no_annual", "no_holidays", "no_annual", "no_school_holidays",
  "no_public_holidays", "no_asc_holidays", "no_Rt_holidays", "no_delay"
)

perf_data <- map(model_versions, read_perf) |>
  map(\(x) select(x, .model, pathogen, location, round_id, origin_date, date, crps_log)) |>
  bind_rows() |>
  bind_rows(
    read_parquet("data/snaive_perf_extended.parquet") |>
      select(.model, pathogen, location, round_id, origin_date, date, crps_log)
  ) |>
  process_perf_data() |>
  filter(horizon >= 1, horizon <= 7 * 4)

# The figures compare models on the dates they all share. Summarising each model separately
# is only equivalent if every model has a score for exactly the same dates, so check that.
n_models <- n_distinct(perf_data$.model)
incomplete <- perf_data |>
  count(pathogen, location, round_id, date) |>
  filter(n != n_models)

if (nrow(incomplete) > 0) {
  print(incomplete)
  stop("Some forecast dates are missing scores for at least one model; see above.")
}

perf_summary <- perf_data |>
  mutate(horizon_week = (horizon - 1) %/% 7) |>
  group_by(.model, pathogen, location, round_id, horizon_week) |>
  summarise(
    first_date = min(date),
    n_days = n(),
    crps_log = mean(crps_log),
    .groups = "drop"
  )

# Every figure averages bin means with equal weight, which matches the per-day mean only
# when every bin is a full week.
stopifnot(all(perf_summary$n_days == 7))

perf_summary |>
  select(-n_days) |>
  arrange(.model, pathogen, location, round_id, horizon_week) |>
  write_parquet("saved_results/perf_summary.parquet")
