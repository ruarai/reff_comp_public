


calculate_performance_metrics_hierarchical <- function(
    round_id, x_pathogen, x_target,
    results_dir = "results"
) {
  
  eval_counts <- read_csv("data/combined-count-smoothed-2026-04-15.csv", show_col_types = FALSE) |> 
    filter(pathogen == x_pathogen, target == x_target) |> 
    
    select(date = notification_date, location, cases_true = value)
  
  round_results <- fs::path(results_dir, round_id)
  
  job_name <- str_c(x_pathogen, "_", x_target)
  
  fit <- read_rds(str_c(round_results, "/", job_name, ".rds"))
  training_counts <- arrow::read_parquet(str_c(round_results, "/", job_name, "_training_counts.parquet"))
  
  date_index <- unique(training_counts$date)
  location_index <- unique(training_counts$location)
  
  perf_data <- fit |> 
    spread_draws(cases[i, ix_t]) |> 
    mutate(date = date_index[ix_t],
           location = location_index[i],
           # Replace NA values for RNG cases with zero (rng freaks out at high phi)
           cases = replace_na(cases, 0)) |> 
    left_join(eval_counts) |> 
    drop_na(cases_true) |> 
    
    mutate(is_overpredict = cases > cases_true) |>
    group_by(date, location) |> 
    mutate(log_cases = log(cases + 1),
           log_cases_true = log(cases_true + 1)) |> 
    summarise(
      cases_true = cases_true[1], log_cases_true = log_cases_true[1],
      bias = sum(is_overpredict) / n(),
      crps = scoringRules::crps_sample(cases_true, cases),
      crps_log = scoringRules::crps_sample(log_cases_true, log_cases),
      med_ae = abs(median(cases) - cases_true),
      med_ae_log = abs(median(log_cases) - log_cases_true),
      .groups = "drop"
    ) |> 
    mutate(pathogen = x_pathogen, target = x_target,
           round_id = round_id, origin_date = training_counts$origin_date[[1]],
           .before = 1)
  
  perf_data |> 
    write_parquet(str_c(round_results, "/", job_name, "_performance.parquet"))
}


make_fit_plots_hierarchical <- function(round_id, x_pathogen, x_target, results_dir = "results") {
  round_results <- fs::path(results_dir, round_id)
  
  job_name <- str_c(x_pathogen, "_", x_target)
  
  fit <- read_rds(str_c(round_results, "/", job_name, ".rds"))
  training_counts <- read_parquet(str_c(round_results, "/", job_name, "_training_counts.parquet"))
  
  origin_date <- min(training_counts$origin_date)
  
  combined_counts <- read_csv("data/combined-count-smoothed-2026-04-22.csv", show_col_types = FALSE) |> 
    filter(pathogen == x_pathogen, location %in% training_counts$location)
  
  eval_counts <- combined_counts |> 
    filter(notification_date > origin_date)
  
  draws <- fit |>
    tidybayes::tidy_draws() |> 
    
    # NA cases values can cause issues with bayesplot, so
    # drop the columns early
    select(-starts_with("cases"))
  
  pars <- c(
    "scaling_annual", "scaling_noise", "ls_noise", "ls_annual"
  ) |> 
    intersect(fit$metadata()$model_params)
  
  p_hist <- mcmc_hist(draws |> mutate(chain = .chain), pars = pars) +
    theme_bw()
  
  p_trace <- mcmc_trace(draws |> mutate(chain = .chain), pars = pars) +
    theme_bw() +
    
    labs(caption = job_name)
  
  p_hist | p_trace
  
  ggsave(str_c(round_results, "/", job_name, "_plot_MCMC.png"), width = 14, height = 10)
  
  
  p_pairs <- mcmc_pairs(draws |> mutate(chain = .chain), pars = pars, off_diag_args = list(size = 0.3, alpha = 0.5))
  
  ggsave(str_c(round_results, "/", job_name, "_plot_MCMC_pairs.png"), p_pairs, width = 16, height = 16)
  
  school_holidays <- get_school_holidays() |> 
    filter(start_date <= origin_date + 365)
  
  public_holidays <- get_public_holidays() |> 
    filter(date <= origin_date + 365)
  
  date_index <- unique(training_counts$date)
  location_index <- unique(training_counts$location)
  
  p_ann <- fit |> 
    spread_draws(f_annual[i, ix_t]) |> 
    mutate(date = date_index[ix_t], location = location_index[i]) |> 
    ggplot() + 
    fastribbon::stat_fastribbon(aes(x = date, y = f_annual)) +
    
    geom_hline(yintercept = 0, linetype = "22") +
    
    facet_wrap(~location, nrow = 1, scales = "free_y") +
    coord_cartesian(ylim = c(-0.2, 0.2)) +
    
    xlab(NULL) + ylab(NULL) +
    
    scale_fill_brewer() + theme_bw() + theme(legend.position = "none")
  
  p_noise <- fit |> 
    spread_draws(f_noise[i, ix_t]) |> 
    mutate(date = date_index[ix_t], location = location_index[i]) |> 
    ggplot() + 
    fastribbon::stat_fastribbon(aes(x = date, y = f_noise)) +
    
    geom_hline(yintercept = 0, linetype = "22") +
    
    facet_wrap(~location, nrow = 1, scales = "free_y") +
    coord_cartesian(ylim = c(-0.2, 0.2)) +
    
    xlab(NULL) + ylab(NULL) +
    
    scale_fill_brewer() + theme_bw() + theme(legend.position = "none")
  
  
  p_Rt <- fit |> 
    spread_draws(log_R[i, ix_t]) |> 
    mutate(date = date_index[ix_t], location = location_index[i]) |> 
    ggplot() + 
    fastribbon::stat_fastribbon(aes(x = date, y = exp(log_R))) +
    
    geom_hline(yintercept = 1, linetype = "22") +
    
    facet_wrap(~location, nrow = 1, scales = "free_y") +
    coord_cartesian(ylim = c(0.7, 1.3)) +
    
    xlab(NULL) + ylab(NULL) +
    
    scale_fill_brewer() + theme_bw() + theme(legend.position = "none")
  
  p_inf <- fit |> 
    spread_draws(inf_delayed[i, ix_t]) |> 
    mutate(date = date_index[ix_t], location = location_index[i]) |> 
    ggplot() + 
    fastribbon::stat_fastribbon(aes(x = date, y = inf_delayed),
                                linewidth = 0.0) +
    
    geom_line(aes(x = date, y = y_trend),
              training_counts) +
    geom_line(aes(x = notification_date, y = value_trend),
              eval_counts, colour = "orangered") +
    
    facet_wrap(~location, nrow = 1, scales = "free_y") +
    
    xlab(NULL) + ylab(NULL) +
    
    scale_fill_brewer() + theme_bw() + theme(legend.position = "none")
  
  p_joint <- p_ann / p_noise / p_Rt / p_inf
  
  ggsave(str_c(round_results, "/", job_name, "_Rt_inf.png"), p_joint, width = 24, height = 12)
  
  library(pammtools)
  p_forecast <- fit |> 
    spread_draws(cases[i, ix_t]) |> 
    mutate(date = date_index[ix_t], location = location_index[i]) |> 
    filter(date > origin_date - days(7 * 12), date <= origin_date + days(7 * 4)) |> 
    ggplot() + 
    fastribbon::stat_fastribbon(aes(x = date, y = cases), geom = "stepribbon",
                                linewidth = 0.0) +
    
    geom_line(aes(x = date, y = y),
              training_counts |> 
                filter(date > origin_date - days(7 * 12), date <= origin_date + days(7 * 4))) +
    geom_line(aes(x = notification_date, y = value),
              eval_counts |> 
                filter(notification_date <= origin_date + days(7 * 4)), colour = "orangered") +
    
    facet_wrap(~location, nrow = 2, scales = "free_y") +
    
    xlab(NULL) + ylab(NULL) +
    
    scale_fill_brewer() + theme_bw() + theme(legend.position = "none")
  
  ggsave(str_c(round_results, "/", job_name, "_plot_forecast.png"), p_forecast, width = 24, height = 8)
  
  
  plot_data_true_total <- combined_counts |> 
    filter(notification_date >= ymd("2023-01-01")) |> 
    group_by(location, year = floor_date(notification_date, "year")) |> 
    summarise(total_inf = sum(value))
  
  total_inf_data <- fit |> 
    spread_draws(inf[i, ix_t]) |> 
    filter(.draw %% 10 == 0) |> 
    mutate(date = date_index[ix_t], location = location_index[i]) |> 
    group_by(location, year = floor_date(date, "year"), .draw) |> 
    summarise(total_inf = sum(inf))
  
  ggplot(total_inf_data |> filter(year < ymd("2026-01-01"))) + 
    geom_point(aes(x = year, y = total_inf),
               position = position_jitter(width = 50),
               size = 0.5, stroke = 0.5, alpha = 0.5) +
    
    geom_point(aes(x = year, y = total_inf),
               colour = "orangered",
               plot_data_true_total |> filter(year < ymd("2026-01-01"))) +
    facet_wrap(~location, scales = "free_y") +
    
    theme_bw()
  
  
  ggsave(str_c(round_results, "/", job_name, "_plot_total_inf.png"), width = 12, height = 8)
  
  
  fit |> 
    spread_draws(log10_total_inf_year_mean[ix_year], log10_total_inf_offset[i]) |> 
    mutate(year = unique(year(training_counts$date))[ix_year], location = location_index[i]) |> 
    mutate(latent_log10_total_inf_year = log10_total_inf_year_mean + log10_total_inf_offset) |> 
    
    ggplot() + 
    fastribbon::stat_fastribbon(aes(x = year, y = latent_log10_total_inf_year),
                                linewidth = 0.0) + 
    
    facet_wrap(~location) +
    
    scale_fill_brewer() + theme_bw() + theme(legend.position = "none")
  
  ggsave(str_c(round_results, "/", job_name, "_plot_latent_total_inf.png"), width = 12, height = 8)
}

