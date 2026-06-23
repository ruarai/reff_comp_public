

# TODO
# rerun this
# run the snaive forecasts
# Plot outputs, in f_2_performance, check if they make sense


model_version <- unique(jobs$model_version)


for(m in model_version) {
  perf_old <- read_parquet(glue::glue("saved_results/{m}_performance.parquet"))
  
  perf_filt <- perf_old |> 
    filter(!(location == "VIC" & round_id %in% seq(ymd("2025-04-10"), ymd("2025-05-01"), by = "7 days")))
  
  
  perf_files_new <- fs::dir_ls(glue::glue("saved_results/{m}/"), recurse = TRUE, glob = "*VIC_case_incidence_performance.parquet")
  
  perf_new <- collect(open_dataset(perf_files_new)) |> 
    filter(round_id %in% seq(ymd("2025-04-10"), ymd("2025-05-01"), by = "7 days"))
  
  perf_filt |> 
    bind_rows(perf_new) |> 
    write_parquet(glue::glue("saved_results/{m}_performance.parquet"))
}