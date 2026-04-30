

results_name <- "results_no_holidays"


perf_files <- fs::dir_ls(str_c("saved_results/", results_name), glob = "*_performance.parquet", recurse = TRUE)

perf_dataset <- open_dataset(perf_files)

perf_dataset |> write_parquet(str_c("saved_results/", results_name, "_performance.parquet"))

