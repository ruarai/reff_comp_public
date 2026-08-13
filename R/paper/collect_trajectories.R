
source("R/functions.R")

trajectory_files <- fs::dir_ls("saved_results/baseline/", recurse = TRUE, glob = "*-trajectories.parquet")



traj_quants <- map(
  trajectory_files,
  function(x) {
    open_dataset(x) |> 
      filter(target == "case_incidence") |> 
      group_by(location, pathogen, target, origin_date, round_id, horizon) |> 
      summarise(lower_95 = quantile(value, 0.025),
                upper_95 = quantile(value, 0.975),
                
                lower_80 = quantile(value, 0.1),
                upper_80 = quantile(value, 0.9),
                
                lower_50 = quantile(value, 0.25),
                upper_50 = quantile(value, 0.75)) |> 
      collect()
  },
  .progress = TRUE
)

traj_quants |> 
  bind_rows() |> 
  write_parquet("saved_results/baseline_quants.parquet")



