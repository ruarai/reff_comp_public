
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


download_round <- function(round_id, results_dir = "baseline") {
  library(processx)

  x <- run(
    "/apps/unimelb-mf-clients/current/bin/unimelb-mf-download",
    c("--mf.config", "data/mflux.cfg",
      "--nb-workers", "4",
      "--out", glue::glue("saved_results/{results_dir}"),
      glue::glue("/projects/proj-5050_reffcomp_model_evaluation-1128.4.1617/{results_dir}/{round_id}")
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

download_traj <- function() {
  library(processx)
  
  
  jobs <- make_retrospective_jobs() |> 
    anti_join(
      ignore_jobs <- tribble(
        ~location, ~pathogen, ~target,
        "NZ", "RSV", "hospital_incidence",
        "NZ", "flu", "hospital_incidence",
      ),
      by = join_by(location, pathogen, target)
    ) |> 
    expand_grid(model_version = "baseline")
  
  jobs |> 
    distinct(round_id) |> 
    pmap(function(round_id) {
      jobs |> 
        filter(round_id == !!round_id) |> 
        mutate(path = glue::glue("/projects/proj-5050_reffcomp_model_evaluation-1128.4.1617/baseline/{round_id}/{pathogen}_{location}_forecast-trajectories.parquet")) |> 
        pull(path) |> 
        str_c(collapse = "")
    })
  
  jobs |> 
    pmap(function(location, pathogen, target, round_id, ...) {
      
      fs::dir_create(glue::glue("saved_results/baseline/{round_id}/"))
      
      if(!fs::file_exists(glue::glue("saved_results/baseline/{round_id}/{pathogen}_{location}_forecast-trajectories.parquet"))) {
        
        x <- run(
          "/apps/unimelb-mf-clients/current/bin/unimelb-mf-download",
          c("--mf.config", "data/mflux.cfg",
            "--nb-workers", "4",
            "--out", glue::glue("saved_results/baseline/{round_id}/"),
            glue::glue("/projects/proj-5050_reffcomp_model_evaluation-1128.4.1617/baseline/{round_id}/{pathogen}_{location}_forecast-trajectories.parquet")
          ),
          echo = TRUE
        )
      }
      
    })
  
  
  
}