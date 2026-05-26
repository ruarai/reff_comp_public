

write_forecast_file_retro <- function(round_id, x_pathogen, x_location, x_target) {
  
  x_target <- "case_incidence"
  
  round_results <- fs::path("saved_results", "results_baseline", round_id)
  
  if(!file.exists(str_c(round_results, "/", x_pathogen, "_", x_location, ".rds"))) {
    print("No fit file found")
    return(NULL)
  }
  
  fit <- read_rds(str_c(round_results, "/", x_pathogen, "_", x_location, ".rds"))
  training_counts <- read_parquet(str_c(round_results, "/", x_pathogen, "_", x_location, "_training_counts.parquet"))
  
  origin_date <- training_counts$origin_date[[1]]
  
  case_trajectories <- fit |> 
    spread_draws(cases[ix_t], latent_cases_log[ix_t], inf_delayed[ix_t]) |> 
    mutate(date = training_counts$date[ix_t],
           latent_cases = exp(latent_cases_log))
  
  case_forecasts <- case_trajectories |>
    ungroup() |>
    transmute(
      round_id = ymd(round_id),
      origin_date = origin_date,
      target = x_target,
      horizon = as.integer(date - origin_date),
      location = x_location,
      pathogen = x_pathogen,
      output_type = "sample",
      output_type_id = .draw,
      value = as.integer(replace_na(cases, 0))
    ) |>
    filter(horizon >= -7, horizon <= 7 * 4)

  
  output_file <- str_c(round_results, "/", x_pathogen, "_", x_location, "_forecast-trajectories.parquet")
  
  arrow::write_parquet(case_forecasts, output_file)
}



jobs <- make_retrospective_jobs() |> 
  filter(round_id >= "2025-05-08", round_id <= "2025-10-09")


pmap(jobs[153:nrow(jobs),], function(round_id, location, pathogen) {
  print(glue::glue("{round_id}-{location}-{pathogen}"))
  
  write_forecast_file_retro(round_id, pathogen, location)
})


traj_files <- fs::dir_ls("saved_results/results_baseline/", recurse = TRUE, glob = "*forecast-trajectories.parquet")

trajectories <- open_dataset(traj_files)

trajectories |> write_parquet(
  "data/uom-reff_comp_2025_retrospectives.parquet", 
  compression = "gzip", compression_level = 9
)
