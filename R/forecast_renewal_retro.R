
source("R/functions.R")

args <- commandArgs(trailingOnly = TRUE)
ix_job <- as.integer(args[1])


jobs <- make_retrospective_jobs() |> 
  anti_join(
    ignore_jobs <- tribble(
      ~location, ~pathogen, ~target,
      "NZ", "RSV", "hospital_incidence",
      "NZ", "flu", "hospital_incidence",
    ),
    by = join_by(location, pathogen, target)
  ) |> 
  expand_grid(model_version = c("baseline_early")) |> 
  
  filter(round_id < ymd("2025-04-10"))

print(str_c("Generated ", nrow(jobs), " jobs"))

if(ix_job > nrow(jobs)) {
  print("Exiting without completing a job (ix_job > nrow(jobs))")
  q('no')
}

job <- jobs[ix_job, ]

x_location <- job$location[[1]]
x_pathogen <- job$pathogen[[1]]
x_target <- job$target[[1]]
x_version <- job$model_version[[1]]
round_id <- job$round_id[[1]]

results_dir <- str_c("saved_results/", x_version)

print(str_c("Performing job ", ix_job, " (", x_location, "-", x_pathogen, "-", x_target, "), saving to ", results_dir))

year <- as.character(year(round_id))
combined_counts <- read_csv(str_c("../forecast-hub/processed-data/", year, "/", round_id, "/combined-count-smoothed-", round_id, ".csv"))


origin_date <- get_origin_date(combined_counts, x_pathogen, x_location, x_target)

n_days_ahead <- 365

training_counts <- process_counts(combined_counts, x_location, x_pathogen, x_target, n_days_ahead, origin_date)


data_list <- create_stan_data(training_counts, n_days_ahead, x_pathogen, x_location)

model_file <- "src/gp_renewal.stan"

fit_file <- str_c(results_dir, "/", round_id, "/", x_pathogen, "_", x_location, "_", x_target, ".rds")

fit <- fit_model(data_list, model_file = model_file)

print("Saving fit object and training data")

fs::dir_create(fs::path(results_dir, round_id))

fit$save_object(fit_file)
arrow::write_parquet(training_counts, str_c(results_dir, "/", round_id, "/", x_pathogen, "_", x_location, "_", x_target, "_training_counts.parquet"))


print("Writing trajectories...")

write_forecast_file(round_id, x_pathogen, x_location, x_target, results_dir = results_dir)

print("Calculating performance metrics...")

calculate_performance_metrics(round_id, x_pathogen, x_location, x_target, results_dir = results_dir)

print("Making plots...")

make_fit_plots(round_id, x_pathogen, x_location, x_target, results_dir = results_dir)

print("Complete")


