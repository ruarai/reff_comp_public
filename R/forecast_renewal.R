

library(tidyverse)

source("R/functions.R")

args <- commandArgs(trailingOnly = TRUE)

ix_job <- as.integer(args[1])
round_id <- str_trim(args[2])

combined_counts <- read_csv(str_c("data/combined-count-", round_id, ".csv"))
combined_counts_smoothed <- read_csv(str_c("data/combined-count-smoothed-", round_id, ".csv"))

all_counts <- combined_counts |> 
  left_join(combined_counts_smoothed)

jobs <- distinct(all_counts, location, pathogen)

print(str_c("Generated ", nrow(jobs), " jobs from combined_counts_smoothed"))

if(ix_job > nrow(jobs)) {
  print("Exiting without completing a job (ix_job > nrow(jobs))")
}

job <- jobs[ix_job, ]

x_location <- job$location[[1]]
x_pathogen <- job$pathogen[[1]]

print(str_c("Performing job ", ix_job, " (", x_location, "-", x_pathogen, ")"))

training_dates <- get_training_dates()

combined_counts_outliers <- all_counts |> 
  left_join(training_dates) |> 
  filter(notification_date >= cutoff_date) |> 
  filter(location == x_location, pathogen == x_pathogen) |> 
  mutate(outlier = abs(log(cases_trend_and_remainder + 0.5) - log(cases_trend + 0.5)) > 1.38) |> 
  
  filter(
    notification_date < max(notification_date)
  )

print(str_c("Dropping ", sum(combined_counts_outliers$outlier), " outliers"))

n_days_ahead <- 365

training_counts <- combined_counts_outliers |> 
  filter(!outlier) |> 
  make_training_counts(x_location, x_pathogen, days_ahead = n_days_ahead)


n_days_data <- nrow(training_counts) - n_days_ahead


library(cmdstanr)

gi <- get_gi(x_pathogen)

model <- cmdstan_model("src/gp_renewal.stan")

# TODO: how to pick these M values? see paper on HSGP
data_list <- list(
  x = training_counts$x, y = training_counts$y[1:n_days_data],
  N = nrow(training_counts), N_data = n_days_data,
  M_annual = 40, M_noise = 40, c_factor = 1.2,
  
  gi = gi, N_gi = length(gi)
)

fit <- model$sample(
  data = data_list,
  seed = 1,
  chains = 4,
  parallel_chains = 4,
  refresh = 50
)

print("Saving fit object and training data")

fs::dir_create(fs::path("results", round_id))

fit$save_object(str_c("results/", round_id, "/", x_pathogen, "_", x_location, ".rds"))
arrow::write_parquet(training_counts, str_c("results/", round_id, "/", x_pathogen, "_", x_location, "_training_counts.parquet"))

print("Complete")



