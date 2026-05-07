
source("R/functions.R")
source("R/exp/process_data_hierarchical.R")
source("R/exp/results_hierarchical.R")


args <- commandArgs(trailingOnly = TRUE)
ix_job <- as.integer(args[1])


results_dir <- "saved_results/hierarchical"

jobs <- make_retrospective_jobs() |> 
  filter(round_id >= "2025-05-08", round_id <= "2025-10-09") |> 
  
  filter(location != "SA", location != "NZ") |> 
  filter(target == "case_incidence") |> 
  distinct(round_id, pathogen)


print(str_c("Generated ", nrow(jobs), " jobs"))

if(ix_job > nrow(jobs)) {
  print("Exiting without completing a job (ix_job > nrow(jobs))")
  q('no')
}


job <- jobs[ix_job, ]

x_target <- "case_incidence"
x_pathogen <- job$pathogen[[1]]
round_id <- job$round_id[[1]]


year <- as.character(year(round_id))
combined_counts <- read_csv(str_c("../forecast-hub/processed-data/", year, "/", round_id, "/combined-count-smoothed-", round_id, ".csv")) |> 
  filter(location %in% c("ACT", "NSW", "WA", "QLD", "VIC", "TAS"))


origin_dates <- combined_counts |> 
  filter(pathogen == x_pathogen, target == x_target,
         is_complete) |> 
  group_by(location) |> 
  filter(notification_date == max(notification_date)) |> 
  arrange(location) |> 
  select(location, origin_date = notification_date) |> 
  ungroup()

n_days_ahead <- as.numeric(ymd("2027-01-01") - max(origin_dates$origin_date)) - 1

training_counts <- process_counts_hierarchical(combined_counts, x_pathogen, x_target, n_days_ahead, origin_dates) |> 
  mutate(year = match(year(date), unique(year(date))))



year_ranges <- training_counts |> 
  ungroup() |> filter(location == location[[1]]) |> 
  mutate(ix = row_number()) |>
  group_by(year) |> 
  summarise(ix_year_start = min(ix), ix_year_end = max(ix))

data_list <- create_stan_data_hierarchical(training_counts, n_days_ahead, x_pathogen)
data_list$N_years <- max(training_counts$year)
data_list$year_start <- year_ranges$ix_year_start
data_list$year_end <- year_ranges$ix_year_end
data_list$sigma_log10_total_inf <- 0.1

library(cmdstanr)
model <- cmdstan_model("src/gp_renewal_hierarchical.stan")



fit <- model$sample(
  data = data_list,
  seed = 1,
  chains = 4,
  parallel_chains = 4,
  init = 0.01,
  refresh = 50,
  iter_warmup = 1000, iter_sampling = 1000
)

print("Saving fit object and training data")

fs::dir_create(fs::path(results_dir, round_id))

fit_file <- str_c(results_dir, "/", round_id, "/", x_pathogen, "_", x_target, ".rds")
fit$save_object(fit_file)
arrow::write_parquet(training_counts, str_c(results_dir, "/", round_id, "/", x_pathogen, "_", x_target, "_training_counts.parquet"))

# print("Writing trajectories...")
# write_forecast_file(round_id, x_pathogen, x_location, x_target, results_dir = results_dir)

print("Calculating performance metrics...")

calculate_performance_metrics_hierarchical(round_id, x_pathogen, x_target, results_dir = results_dir)

print("Making plots...")

make_fit_plots_hierarchical(round_id, x_pathogen, x_target, results_dir = results_dir)

print("Complete")
