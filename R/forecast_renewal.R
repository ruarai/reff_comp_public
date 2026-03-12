


source("R/functions.R")

args <- commandArgs(trailingOnly = TRUE)

ix_job <- as.integer(args[1])
round_id <- str_trim(args[2])

combined_counts <- read_csv(str_c("../forecast-hub/processed-data/2025/", round_id, "/combined-count-", round_id, ".csv"))

jobs <- distinct(combined_counts, location, pathogen)

print(str_c("Generated ", nrow(jobs), " jobs from combined_counts"))

if(ix_job > nrow(jobs)) {
  print("Exiting without completing a job (ix_job > nrow(jobs))")
  q('no')
}

job <- jobs[ix_job, ]

x_location <- job$location[[1]]
x_pathogen <- job$pathogen[[1]]

origin_date <- get_origin_date(combined_counts, x_pathogen, x_location)

print(str_c("Performing job ", ix_job, " (", x_location, "-", x_pathogen, ")"))

n_days_ahead <- 365

training_counts <- process_counts(combined_counts, x_location, x_pathogen, n_days_ahead) |> 
  # Save some extra data in case we need it later
  mutate(origin_date = origin_date)

data_list <- create_stan_data(training_counts, n_days_ahead, x_pathogen, x_location)

fit <- fit_model(data_list)

print("Saving fit object and training data")

fs::dir_create(fs::path("results", round_id))

fit_file <- str_c("results/", round_id, "/", x_pathogen, "_", x_location, ".rds")
fit$save_object(fit_file)
arrow::write_parquet(training_counts, str_c("results/", round_id, "/", x_pathogen, "_", x_location, "_training_counts.parquet"))

print("Writing trajectories...")

write_forecast_file(round_id, x_pathogen, x_location)

print("Making plots...")

make_fit_plots(round_id, x_pathogen, x_location)

print("Complete")



