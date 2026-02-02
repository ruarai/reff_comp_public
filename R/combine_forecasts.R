
source("R/functions.R")

forecasts <- open_dataset("forecast_outputs_raw/")

model_name <- "baseline"

fs::dir_create(fs::path("forecast_outputs", model_name))

rounds <- forecasts |> distinct(round_id) |> collect() |> pull(round_id) |> as.character()

for(x_round in rounds) {
  print(x_round)
  forecasts |> 
    filter(round_id == x_round) |> 
    collect() |> 
    write_parquet(fs::path("forecast_outputs", "baseline", str_c(x_round, "-uom-reff-composure-", model_name,".parquet")))
}

