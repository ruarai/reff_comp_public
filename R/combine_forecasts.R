
source("R/functions.R")

forecasts <- open_dataset("forecast_outputs_raw/")

rounds <- forecasts |> distinct(round_id) |> collect() |> pull(round_id) |> as.character()

for(x_round in rounds) {
  forecasts |> 
    filter(round_id == x_round) |> 
    collect() |> 
    write_parquet(fs::path("forecast_outputs", str_c(x_round, "-uom-reff-composure-2.parquet")))
}

