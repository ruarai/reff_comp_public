

source("R/functions.R")
source("R/paper/common.R")

run_name <- "no_Rt_holidays"

rhat_files <- fs::dir_ls(glue::glue("saved_results/{run_name}/"), recurse = TRUE, glob = "*_convergence.parquet")
fit_files <- fs::dir_ls(glue::glue("saved_results/{run_name}/"), recurse = TRUE, glob = "*.rds")

length(rhat_files) == length(fit_files)

rhat_tbl <- map(rhat_files, read_parquet) |> 
  bind_rows(.id = "id")

max(rhat_tbl$rhat, na.rm = TRUE)


