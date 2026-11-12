
source("R/functions.R")
source("R/paper/common.R")

fit_files <- fs::dir_ls("saved_results/baseline/2025-11-27/", recurse = TRUE, glob = "*.rds")


holiday_draws <- map(
  fit_files,
  function(fit_file) {
    
    round_id <- ymd(str_extract(fit_file, "\\d{4}-\\d{2}-\\d{2}"))
    job_name <- fs::path_file(fit_file) |> fs::path_ext_remove()
    pathogen <- str_split_fixed(job_name, "_", 3)[1]
    location <- str_split_fixed(job_name, "_", 3)[2]
    target <- str_split_fixed(job_name, "_", 3)[3]
    
    if(!(location %in% locations_paper)) {
      return(NULL)
    }
    
    fit <- read_rds(fit_file)
    
    
    fit |> 
      
      spread_draws(school_holiday_Rt_effect_raw[ix_s],
                   school_holiday_asc_effect_raw[ix_s],
                   scaling_annual, scaling_noise # Include for comparison
                   ) |> 
      mutate(location = !!location, pathogen = !!pathogen, target = !!target)
  },
  .progress = TRUE
) |> 
  bind_rows()

write_parquet(holiday_draws, "data/paper/holiday_draws.parquet")


public_holiday_draws <- map(
  fit_files,
  function(fit_file) {
    
    round_id <- ymd(str_extract(fit_file, "\\d{4}-\\d{2}-\\d{2}"))
    job_name <- fs::path_file(fit_file) |> fs::path_ext_remove()
    pathogen <- str_split_fixed(job_name, "_", 3)[1]
    location <- str_split_fixed(job_name, "_", 3)[2]
    target <- str_split_fixed(job_name, "_", 3)[3]
    
    if(!(location %in% locations_paper)) {
      return(NULL)
    }
    
    fit <- read_rds(fit_file)
    
    public_holidays <- get_public_holidays() |> 
      filter(location == !!location)
    
    fit |> 
      
      spread_draws(public_holiday_Rt_effect_raw[ix_p],
                   public_holiday_asc_effect_raw[ix_p]
      ) |> 
      mutate(location = !!location, pathogen = !!pathogen, target = !!target,
             public_holiday_name = unique(public_holidays$holiday_name)[ix_p])
  },
  .progress = TRUE
) |> 
  bind_rows()

write_parquet(public_holiday_draws, "data/paper/public_holiday_draws.parquet")


gp_draws <- map(
  fit_files,
  function(fit_file) {
    
    round_id <- ymd(str_extract(fit_file, "\\d{4}-\\d{2}-\\d{2}"))
    job_name <- fs::path_file(fit_file) |> fs::path_ext_remove()
    pathogen <- str_split_fixed(job_name, "_", 3)[1]
    location <- str_split_fixed(job_name, "_", 3)[2]
    target <- str_split_fixed(job_name, "_", 3)[3]
    
    if(!(location %in% locations_paper)) {
      return(NULL)
    }
    
    fit <- read_rds(fit_file)
    
    training_counts <- read_parquet(
      str_c("saved_results/baseline/", round_id, "/", job_name, "_training_counts.parquet")
    )
    
    fit |> 
      spread_draws(f_annual[ix_t], f_noise[ix_t], log_R[ix_t]) |> 
      filter(.draw %% 10 == 0) |>
      mutate(date = training_counts$date[ix_t]) |>
      # Figures 3 and 5 only show 2025
      filter(year(date) == 2025) |>
      mutate(location = location, pathogen = pathogen) |>
      select(.draw, date, location, pathogen, log_R, f_annual, f_noise)
  },
  .progress = TRUE
) |>
  bind_rows()

# Store draws as 32-bit floats to keep the file under GitHub's 100 MB limit
gp_draws |>
  mutate(.draw = as.integer(.draw)) |>
  arrow_table() |>
  (\(t) t$cast(schema(
    .draw = int32(), date = date32(), location = utf8(), pathogen = utf8(),
    log_R = float32(), f_annual = float32(), f_noise = float32()
  )))() |>
  write_parquet("data/paper/gp_draws.parquet")


delay_draws <- map(
  fit_files,
  function(fit_file) {
    
    round_id <- ymd(str_extract(fit_file, "\\d{4}-\\d{2}-\\d{2}"))
    job_name <- fs::path_file(fit_file) |> fs::path_ext_remove()
    pathogen <- str_split_fixed(job_name, "_", 3)[1]
    location <- str_split_fixed(job_name, "_", 3)[2]
    target <- str_split_fixed(job_name, "_", 3)[3]
    
    if(!(location %in% locations_paper)) {
      return(NULL)
    }
    
    fit <- read_rds(fit_file)
    
    fit |> 
      
      spread_draws(
        delay_mean, delay_sd,
        delay_dist[days]
      ) |> 
      mutate(location = !!location, pathogen = !!pathogen, target = !!target,
             days = days - 1)
  },
  .progress = TRUE
) |> 
  bind_rows()


write_parquet(delay_draws, "data/paper/delay_draws.parquet")
