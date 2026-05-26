
source("R/functions.R")
fit_files <- fs::dir_ls("saved_results/baseline/2025-11-27/", recurse = TRUE, glob = "*.rds")

holiday_draws <- map(
  fit_files,
  function(fit_file) {
    
    fit <- read_rds(fit_file)
    
    round_id <- ymd(str_extract(fit_file, "\\d{4}-\\d{2}-\\d{2}"))
    job_name <- fs::path_file(fit_file) |> fs::path_ext_remove()
    pathogen <- str_split_fixed(job_name, "_", 3)[1]
    location <- str_split_fixed(job_name, "_", 3)[2]
    target <- str_split_fixed(job_name, "_", 3)[3]
    
    fit |> 
      
      spread_draws(school_holiday_Rt_effect_raw[ix_s],
                   school_holiday_asc_effect_raw[ix_s]) |> 
      mutate(location = !!location, pathogen = !!pathogen, target = !!target)
  },
  .progress = TRUE
) |> 
  bind_rows()

write_parquet(holiday_draws, "data/paper/holiday_draws.parquet")


gp_draws <- map(
  fit_files,
  function(fit_file) {
    
    fit <- read_rds(fit_file)
    
    round_id <- ymd(str_extract(fit_file, "\\d{4}-\\d{2}-\\d{2}"))
    job_name <- fs::path_file(fit_file) |> fs::path_ext_remove()
    pathogen <- str_split_fixed(job_name, "_", 3)[1]
    location <- str_split_fixed(job_name, "_", 3)[2]
    target <- str_split_fixed(job_name, "_", 3)[3]
    
    training_counts <- read_parquet(
      str_c("saved_results/baseline/", round_id, "/", job_name, "_training_counts.parquet")
    )
    
    fit |> 
      spread_draws(f_annual[ix_t], f_noise[ix_t], log_R[ix_t]) |> 
      filter(.draw %% 10 == 0) |> 
      mutate(date = training_counts$date[ix_t]) |> 
      mutate(round_id = round_id, location = location, pathogen = pathogen, target = target)
  },
  .progress = TRUE
) |> 
  bind_rows()

write_parquet(gp_draws, "data/paper/gp_draws.parquet")
