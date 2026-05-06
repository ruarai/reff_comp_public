
library(tidyverse)
library(tidybayes)
library(bayesplot)
library(patchwork)
library(arrow)


source("R/holidays.R")
source("R/functions_plotting.R")
source("R/functions_fitting.R")


get_training_dates <- function() {
  locations <- c("ACT", "NSW", "NT", "QLD", "TAS", "VIC", "WA", "SA", "NZ")
  training_dates <- bind_rows(
    expand_grid(location = locations, pathogen = "flu", cutoff_date = ymd("2023-01-01")),
    expand_grid(
      location = c("ACT", "NSW", "NT", "QLD", "TAS", "VIC", "WA"),
      pathogen = "SARSCOV2", cutoff_date = ymd("2023-03-01")
    ),
    expand_grid(
      location = c("NZ", "SA"),
      pathogen = "SARSCOV2", cutoff_date = ymd("2024-09-01")
    ),
    expand_grid(location = locations, pathogen = "RSV", cutoff_date = ymd("2023-01-01")),
  )
  
}

make_retrospective_jobs <- function() {
  hub_path <- "../forecast-hub/"
  
  rounds <- fs::dir_ls(fs::path(hub_path, "processed-data", "2025"))
  
  tibble(
    round_dir = rounds
  ) |> 
    mutate(round_id = str_extract(round_dir, "\\d{4}-\\d{2}-\\d{2}"),
           counts_file = fs::path(round_dir, str_c("combined-count-", round_id, ".csv"))) |> 
    
    filter(fs::file_exists(counts_file)) |> 
    
    rowwise() |> 
    mutate(counts = list(read_csv(counts_file, show_col_types = FALSE)), 
           jobs = list(distinct(counts, location, pathogen, target))) |> 
    select(-counts) |> 
    unnest(jobs) |> 
    
    select(round_id, location, pathogen, target)
}

calculate_performance_metrics <- function(
    round_id, x_pathogen, x_location, x_target,
    results_dir = "results"
) {
  
  eval_counts <- read_csv("data/combined-count-smoothed-2026-04-15.csv", show_col_types = FALSE) |> 
    filter(location == x_location, pathogen == x_pathogen, target == x_target) |> 
    process_counts(x_location, x_pathogen, x_target, 0, today()) |> 
    
    select(date, cases_true = y)
  
  round_results <- fs::path(results_dir, round_id)
  
  job_name <- str_c(x_pathogen, "_", x_location, "_", x_target)
  
  fit <- read_rds(str_c(round_results, "/", job_name, ".rds"))
  training_counts <- arrow::read_parquet(str_c(round_results, "/", job_name, "_training_counts.parquet"))
  
  perf_data <- fit |> 
    spread_draws(latent_cases_weekly_log[ix_t], cases[ix_t], phi) |> 
    mutate(date = training_counts$date[ix_t],
           
           # Replace NA values for RNG cases with zero (rng freaks out at high phi)
           cases = replace_na(cases, 0)) |> 
    left_join(eval_counts) |> 
    drop_na(cases_true) |> 
    
    mutate(log_likelihood = dnbinom(cases_true, mu = exp(latent_cases_weekly_log), size = phi, log = TRUE),
           is_overpredict = cases > cases_true) |>
    group_by(date) |> 
    mutate(log_cases = log(cases + 1),
           log_cases_true = log(cases_true + 1)) |> 
    summarise(
      cases_true = cases_true[1], log_cases_true = log_cases_true[1],
      elpd = matrixStats::logSumExp(log_likelihood) - log(n()),
      bias = sum(is_overpredict) / n(),
      crps = scoringRules::crps_sample(cases_true, cases),
      crps_log = scoringRules::crps_sample(log_cases_true, log_cases),
      med_ae = abs(median(cases) - cases_true),
      med_ae_log = abs(median(log_cases) - log_cases_true),
      .groups = "drop"
    ) |> 
    mutate(location = x_location, pathogen = x_pathogen, target = x_target,
           round_id = round_id, origin_date = training_counts$origin_date[[1]],
           .before = 1)
  
  perf_data |> 
    write_parquet(str_c(round_results, "/", job_name, "_performance.parquet"))
}

label_pathogens <- function(df) {
  df |> mutate(
    pathogen = factor(pathogen, c("RSV", "flu", "SARSCOV2"), c("RSV", "Influenza", "SARS-CoV-2"))
  ) 
}

combine_forecasts <- function(round_id, dir = "results") {
  forecast_trajectories <- fs::dir_ls(fs::path(dir, round_id), glob = "*forecast-trajectories.parquet")
  
  if(length(forecast_trajectories) == 0) {
    stop("No forecast trajectories found.")
  }
  
  forecasts <- open_dataset(forecast_trajectories)
  
  model_name <- "uom-reff_comp"
  output_path <- fs::path("forecast_outputs", model_name)
  
  fs::dir_create(output_path)
  
  forecasts_out <- forecasts |>
    filter(horizon <= 7 * 8,
           output_type_id <= 4000) |> 
    collect() |> 
    mutate(horizon = as.integer(horizon),
           value = as.integer(replace_na(value, 0)))
  
  write_parquet(forecasts_out, fs::path(output_path, str_c(round_id, "-", model_name, ".parquet")))
}

save_plots_to_pdf <- function(plot_list, file_path, width = 8, height = 6) {
  pdf(file = file_path, width = width, height = height)
  
  for (plot in plot_list) {
    print(plot)
  }
  
  dev.off()
}

check_and_log_git_status <- function(round_id) {
  library(gert)
  is_clean <- nrow(git_status()) == 0
  if (!is_clean) {
    print(git_status())
    stop("Working directory is not clean.")
  } else {
    fs::dir_create(fs::path("results", round_id))
    write_file(git_info()$commit, fs::path("results", round_id, "git_commit"))
  }
}


make_intervals <- function(
    tbl,
    vars,
    sample_var = rlang::sym(".draw"),
    date_var = rlang::sym("date"),
    group_vars = "",
    probs = c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9),
    na.rm = FALSE
) {
  tbl |> 
    ungroup() |> 
    select(!!date_var, !!sample_var, all_of(vars), any_of(group_vars)) |> 
    pivot_longer(
      cols = c(-!!sample_var, -!!date_var, -any_of(group_vars)), 
      names_to = "variable", 
      values_to = "value"
    ) |> 
    pivot_wider(
      names_from = !!sample_var, 
      names_prefix = "sample_", 
      values_from = value
    ) |> 
    make_intervals_internal(probs = probs, na.rm = na.rm)
}

make_intervals_internal <- function(
    tbl, 
    probs = c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), 
    prefix = "sample_", 
    na.rm = FALSE
) {
  data_matrix <- tbl |> 
    select(starts_with(prefix)) |> 
    as.matrix()
  
  id_tbl <- tbl |> 
    select(!starts_with(prefix))
  
  medians <- data_matrix |> 
    matrixStats::rowMedians(na.rm = na.rm)
  
  quant_probs <- c(rev(1 - probs) / 2, 0.5 + probs / 2)
  quant_labels <- as.character(probs * 100)
  
  quant_names <- c(
    str_c("lower_", rev(quant_labels)), 
    str_c("upper_", quant_labels)
  )
  
  quants <- data_matrix |> 
    matrixStats::rowQuantiles(probs = quant_probs, na.rm = na.rm)
  
  colnames(quants) <- quant_names
  
  bind_cols(id_tbl, median = medians, as_tibble(quants)) |> 
    pivot_longer(
      cols = all_of(quant_names),
      names_to = c(".value", "quant"),
      names_pattern = "(lower|upper)_(.*)"
    ) |> 
    mutate(
      quant = factor(quant, levels = quant_labels) |> fct_rev()
    )
}

add_forecast_noise <- function(case_trajectories, training_counts, horizon_range = c(-7, 7 * 8)) {
  fit_data <- training_counts |> 
    filter(date < origin_date) |> 
    filter(date > max(date) - days(365)) |> 
    mutate(remainder = y_trend_and_remainder - y_trend,
           log_abs_remainder = log(abs(remainder) + 1),
           count = y_trend)
  
  require(mgcv)
  
  gam_fit <- gam(
    log_abs_remainder ~ s(count),
    data = fit_data
  )
  
  smoothed_forecast <- case_trajectories |> 
    mutate(horizon = as.integer(date - origin_date)) |> 
    filter(horizon >= horizon_range[1], horizon <= horizon_range[2]) |>
    rename(count = inf_delayed)
  
  forecast_noisy <- smoothed_forecast |> 
    ungroup() |> 
    mutate(
      pred_noise_mean = predict(gam_fit, newdata = smoothed_forecast),
      pred_noise_mean = pmax(1e-5, pmin(6, pred_noise_mean)),
      pred_noise_mean = replace_na(pred_noise_mean, 1.0),
      
      pred_noise = VGAM::rlaplace(n(), location = 0, scale = exp(pred_noise_mean)),
      count = pmax(0, round(count) + round(pred_noise))
    )
  
  return(forecast_noisy)
}
