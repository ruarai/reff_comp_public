
source("R/functions.R")

source("R/presentation/plot_theme.R")


args <- commandArgs(trailingOnly = TRUE)

ix_job <- as.integer(args[1])

jobs <- make_retrospective_jobs() |> 
  filter(round_id >= "2025-05-08", round_id <= "2025-08-01")


if(ix_job > nrow(jobs)) {
  print("Exiting without completing a job (ix_job > nrow(jobs))")
  q('no')
}

print(str_c("Generated ", nrow(jobs), " jobs"))

calc_energy_distance <- function(x, y) {
  z <- sort(unique(c(x, y)))
  
  cdf_x <- ecdf(x)(z)
  cdf_y <- ecdf(y)(z)
  
  dz <- diff(z)
  
  integral <- sum((cdf_x[-length(z)] - cdf_y[-length(z)])^2 * dz)
  
  return(2 * integral)
}


job <- jobs[ix_job, ]

x_location <- job$location[[1]]
x_pathogen <- job$pathogen[[1]]
round_id <- job$round_id[[1]]

print(str_c("Performing job ", ix_job, " (", round_id, ": ", x_location, "-", x_pathogen, ")"))


fit_final <- read_rds(str_c("saved_results/results_baseline/2025-10-09/", x_pathogen, "_", x_location, ".rds"))
training_counts <- read_parquet(str_c("saved_results/results_baseline/2025-10-09/", x_pathogen, "_", x_location, "_training_counts.parquet"))

fit_forecast <- read_rds(str_c("saved_results/results_baseline/", round_id, "/", x_pathogen, "_", x_location, ".rds"))
origin_date <- read_parquet(str_c("saved_results/results_baseline/", round_id, "/", x_pathogen, "_", x_location, "_training_counts.parquet"))$origin_date[[1]]


traj_final <- fit_final |> 
  spread_draws(log_R[ix_t]) |> 
  mutate(date = training_counts$date[ix_t]) |> 
  group_by(date) |> 
  summarise(log_R_final = list(log_R))

traj_forecast <- fit_forecast |> 
  spread_draws(log_R[ix_t]) |> 
  mutate(date = training_counts$date[ix_t]) |> 
  group_by(date) |> 
  summarise(log_R_forecast = list(log_R))


results <- traj_forecast |> 
  inner_join(traj_final, by = "date") |> 
  mutate(
    energy_dist = map2_dbl(log_R_forecast, log_R_final, calc_energy_distance)
  ) |> 
  select(date, energy_dist) |> 
  mutate(location = x_location, pathogen = x_pathogen,
         round_id = round_id, origin_date = origin_date, .before = 1)

results |> write_parquet(
  str_c("saved_results/results_baseline/", round_id, "/",  x_pathogen, "_", x_location, "_energy_score.parquet")
)

print("Complete")


