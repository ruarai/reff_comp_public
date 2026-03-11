
source("R/functions.R")

round_id <- "2025-05-08"

forecast_trajectories <- fs::dir_ls(str_c("results_baseline/", round_id), glob = "*forecast-trajectories.parquet")

forecasts <- open_dataset(forecast_trajectories)

model_name <- "uom-reff-composure"

fs::dir_create(fs::path("forecast_outputs", model_name))

rounds <- forecasts |> distinct(round_id) |> collect() |> pull(round_id) |> as.character()

output_path <- fs::path("forecast_outputs", model_name)

forecasts |> 
  collect() |> 
  write_parquet(fs::path(output_path, str_c(x_round, "-", model_name,".parquet")))


fs::dir_copy(output_path, "../forecast-hub/model-output/")
