

upload_output_to_mediaflux <- function(target_dir) {
  library(processx)
  
  x <- run(
    "/apps/unimelb-mf-clients/current/bin/unimelb-mf-upload",
    c("--mf.config", "data/mflux.cfg",
      "--nb-workers", "4",
      "--aggregate-transmission",
      "--dest", "/projects/proj-5050_reffcomp_model_evaluation-1128.4.1617",
      target_dir),
    echo = TRUE
  )
  
  message(x$stdout)
}


download_dir_from_mediaflux <- function(target_dir) {
  library(processx)

  x <- run(
    "/apps/unimelb-mf-clients/current/bin/unimelb-mf-download",
    c("--mf.config", "data/mflux.cfg",
      "--nb-workers", "4",
      "--out", "saved_results/baseline_2",
      str_c("/projects/proj-5050_reffcomp_model_evaluation-1128.4.1617/baseline_2/2025-07-03")
    ),
    echo = TRUE
  )

  message(x$stdout)
}


consolidate_result <- function(results_name) {
  perf_files <- fs::dir_ls(str_c("saved_results/", results_name), glob = "*_performance.parquet", recurse = TRUE)
  
  perf_dataset <- open_dataset(perf_files)
  
  perf_dataset |> write_parquet(str_c("saved_results/", results_name, "_performance.parquet"))
  
  
  upload_output_to_mediaflux(str_c("saved_results/", results_name))
  
  fs::dir_delete(str_c("saved_results/", results_name))
}


consolidate_result("no_annual")
consolidate_result("no_holidays")
consolidate_result("no_holidays_no_annual")