

source("R/functions.R")
source("R/paper/common.R")

fit_files <- fs::dir_ls("saved_results/baseline/", recurse = TRUE, glob = "*.rds")

library(posterior)

library(future.callr)
library(furrr)
plan(callr, workers = 4)


rhat_tbl <- future_map(
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
    
    # Calculate rhat/ess only for the parameters that are actually free
    # (i.e. not the derived parameters)
    params_diag_match <- c(
      "scaling_annual", "scaling_noise", "scaling_weekly_asc", "ls_annual",
      "ls_noise", "log_R_0", "phi", "beta",
      "beta_offset_weekly", "r_0", "log_I_0",
      "school_holiday", "public_holiday", "delay_mean", "delay_sd"
    )
    
    params_all <- fit$metadata()$model_params
    params_diag <- params_all[str_starts(params_all, str_c(params_diag_match, collapse = "|"))]
    
    draws_summ <- fit$summary(
      variables = params_diag, 
      rhat = rhat, ess_bulk = ess_bulk, mean = mean
    ) |> 
      mutate(location = !!location, pathogen = !!pathogen, target = !!target,
             round_id = !!round_id)
    
    
    
    return(draws_summ)
  },
  .progress = TRUE
) |> 
  bind_rows()

max(rhat_tbl$rhat, na.rm = TRUE)

rhat_tbl |>
  filter(rhat > 1.01) |> 
  distinct(round_id, location, pathogen, target)

rhat_tbl |>
  write_csv("data/paper/draws_summ_all_rounds.csv")
