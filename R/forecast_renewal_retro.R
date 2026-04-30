
source("R/functions.R")

args <- commandArgs(trailingOnly = TRUE)

results_dir <- "saved_results/results_normalised"

ix_job <- as.integer(args[1])

jobs <- make_retrospective_jobs() |> 
  filter(round_id >= "2025-05-08", round_id <= "2025-10-09") |> 
  anti_join(
    ignore_jobs <- tribble(
      ~location, ~pathogen, ~target,
      "NZ", "RSV", "hospital_incidence",
      "NZ", "flu", "hospital_incidence",
    ),
    by = join_by(location, pathogen, target)
  )

print(str_c("Generated ", nrow(jobs), " jobs"))

if(ix_job > nrow(jobs)) {
  print("Exiting without completing a job (ix_job > nrow(jobs))")
  q('no')
}

job <- jobs[ix_job, ]

x_location <- job$location[[1]]
x_pathogen <- job$pathogen[[1]]
x_target <- job$target[[1]]
round_id <- job$round_id[[1]]

origin_date <- get_origin_date(combined_counts, x_pathogen, x_location, x_target)

print(str_c("Performing job ", ix_job, " (", x_location, "-", x_pathogen, "-", x_target, ")"))


year <- as.character(year(round_id))
combined_counts <- read_csv(str_c("../forecast-hub/processed-data/", year, "/", round_id, "/combined-count-smoothed-", round_id, ".csv"))

n_days_ahead <- as.numeric(ymd("2027-01-01") - origin_date)

training_counts <- process_counts(combined_counts, x_location, x_pathogen, x_target, n_days_ahead, origin_date)


data_list <- create_stan_data(training_counts, n_days_ahead, x_pathogen, x_location)

fit <- fit_model(data_list, model_file = "src/gp_renewal_norm.stan")

print("Saving fit object and training data")

fs::dir_create(fs::path(results_dir, round_id))

fit_file <- str_c(results_dir, "/", round_id, "/", x_pathogen, "_", x_location, ".rds")
fit$save_object(fit_file)
arrow::write_parquet(training_counts, str_c(results_dir, "/", round_id, "/", x_pathogen, "_", x_location, "_training_counts.parquet"))

print("Writing trajectories...")

write_forecast_file(round_id, x_pathogen, x_location, x_target, results_dir = results_dir)

print("Making plots...")

make_fit_plots(round_id, x_pathogen, x_location, x_target, results_dir = results_dir)

print("Calculating performance metrics...")

calculate_performance_metrics(round_id, x_pathogen, x_location, x_target, results_dir = results_dir)

print("Complete")


