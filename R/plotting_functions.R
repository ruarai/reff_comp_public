
plot_reff_components <- function(trajectory_data, origin_date, holidays) {
  
  labels <- c(
    "f_annual" = "exp(f_annual)",
    "f_noise" = "exp(f_noise)",
    "log_R" = "Rt"
  )
  alpha_vals <- scales::rescale(rev(1/1.7^(1:8)), to = c(0.05, 0.99))
  
  trajectory_quants <- trajectory_data |> 
    make_intervals(c("f_annual", "f_noise", "log_R")) |> 
    mutate(variable = factor(variable, names(labels), labels = labels))
  
  trajectory_trace <- trajectory_data |> 
    ungroup() |> 
    filter(.draw %% 1000 == 0) |>
    select(date, .draw, f_annual, f_noise, log_R) |> 
    pivot_longer(-c(date, .draw)) |> 
    mutate(variable = factor(name, names(labels), labels = labels)) 
  
  ggplot(trajectory_quants) +
    
    geom_vline(aes(xintercept = date), tibble(date = seq(ymd("2023-01-01"), ymd("2027-01-01"), "1 year")),
               alpha = 0.3) +
    
    geom_rect(aes(xmin = start_date, xmax = end_date, ymin = 0, ymax = Inf),
              holidays, alpha = 0.1) +
    
    geom_ribbon(aes(x = date, ymin = exp(lower), ymax = exp(upper),
                    group = quant, alpha = quant, fill = variable)) +
    
    geom_line(aes(x = date, y = exp(value), group = .draw, colour = variable),
              alpha = 0.3,
              trajectory_trace) +
    
    geom_hline(yintercept = 1, alpha = 0.5) +
    
    facet_grid(rows = vars(variable)) +
    
    xlab("Date") + ylab("Value") +
    
    scale_alpha_manual(values = alpha_vals) +
    
    ggokabeito::scale_colour_okabe_ito(order = c(3, 5, 6)) +
    ggokabeito::scale_fill_okabe_ito(order = c(3, 5, 6)) +
    
    scale_y_log10() +
    
    geom_vline(xintercept = origin_date, linetype = "44") +
    
    theme_bw() +
    theme(legend.position = "none")
}

plot_latent_inf <- function(trajectory_data, training_counts, origin_date, holidays) {
  plot_quants <- make_intervals(trajectory_data, "inf")
  plot_quants_delayed <- make_intervals(trajectory_data, "inf_delayed")
  
  
  alpha_vals <- scales::rescale(rev(1/1.7^(1:8)), to = c(0.05, 0.99))
  
  
  ggplot() +
    
    geom_vline(aes(xintercept = date), tibble(date = seq(ymd("2023-01-01"), ymd("2027-01-01"), "1 year")),
               alpha = 0.4) +
    
    geom_rect(aes(xmin = start_date, xmax = end_date, ymin = -Inf, ymax = Inf),
              holidays, alpha = 0.1) +
    
    geom_ribbon(aes(x = date, ymin = lower, ymax = upper,
                    group = quant, alpha = quant),
                plot_quants,
                fill = ggokabeito::palette_okabe_ito(1)) +
    
    geom_ribbon(aes(x = date, ymin = lower, ymax = upper,
                    group = quant, alpha = quant),
                plot_quants_delayed,
                fill = ggokabeito::palette_okabe_ito(5)) +
    
    geom_line(aes(x = date, y = y_smooth),
              training_counts) +
    
    coord_cartesian(ylim = c(0, max(training_counts$y_smooth, na.rm = TRUE) * 1.5)) +
    
    xlab("Notification date") + ylab("Count") +
    
    scale_alpha_manual(values = alpha_vals) +
    
    theme_bw() +
    
    theme(strip.text = element_blank(),
          legend.position = "none")
}

plot_cases <- function(trajectory_data, training_counts, origin_date, holidays) {
  plot_quants <- make_intervals(trajectory_data, "cases")
  
  
  alpha_vals <- scales::rescale(rev(1/1.7^(1:8)), to = c(0.05, 0.99))
  
  
  ggplot(plot_quants) +
    
    geom_rect(aes(xmin = start_date, xmax = end_date, ymin = -Inf, ymax = Inf),
              holidays, alpha = 0.1) +
    
    geom_ribbon(aes(x = date, ymin = lower, ymax = upper,
                    group = quant, alpha = quant),
                fill = ggokabeito::palette_okabe_ito(5)) +
    
    coord_cartesian(ylim = c(0, max(training_counts$y_smooth, na.rm = TRUE) * 1.5)) +
    
    geom_vline(xintercept = origin_date, linetype = "44") +
    
    geom_line(aes(x = date, y = y),
              training_counts) +
    
    xlab("Notification date") + ylab("Count") +
    
    scale_alpha_manual(values = alpha_vals) +
    
    theme_bw() +
    
    theme(strip.text = element_blank(),
          legend.position = "none")
}

make_forecast_plot <- function(trajectory_data, training_counts, eval_counts, origin_date, holidays, x_location, x_pathogen) {
  
  forecast_plot_window <- c(origin_date - days(7 * 4), origin_date + days(7 * 12))
  
  p_forecast_latent <- plot_latent_inf(trajectory_data, training_counts, origin_date, holidays) +
    
    xlim(forecast_plot_window) +
    
    geom_line(aes(x = date, y = y_smooth), colour = "orangered",
              eval_counts) +
    
    geom_vline(xintercept = origin_date + days(4 * 7), linetype = "24")
  
  p_forecast_cases <- plot_cases(trajectory_data, training_counts, origin_date, holidays) +
    
    xlim(forecast_plot_window) +
    
    geom_line(aes(x = date, y = y), colour = "orangered",
              eval_counts) +
    
    geom_vline(xintercept = origin_date + days(4 * 7), linetype = "24") +
    
    labs(caption = str_c(x_location, " - ", x_pathogen))
  
  
  p_forecast_components <- plot_reff_components(trajectory_data, origin_date, holidays) +
    
    xlim(forecast_plot_window) +
    
    geom_vline(xintercept = origin_date + days(4 * 7), linetype = "24")
  
  p_forecast_components / p_forecast_latent / p_forecast_cases
}


make_fit_plots <- function(round_id, x_pathogen, x_location) {
  round_results <- fs::path("results", round_id)
  
  fit <- read_rds(str_c(round_results, "/", x_pathogen, "_", x_location, ".rds"))
  training_counts <- arrow::read_parquet(str_c(round_results, "/", x_pathogen, "_", x_location, "_training_counts.parquet"))
  
  
  origin_date <- training_counts |> 
    drop_na(y) |> 
    pull(date) |> 
    max()
  
  eval_counts <- read_csv("data/combined-count-2026-01-15.csv") |> 
    group_by(location, pathogen) |> 
    mutate(cases_trend_and_remainder = zoo::rollmean(cases, k = 7, fill = NA),
           outlier = FALSE) |> 
    filter(location == x_location, pathogen == x_pathogen) |> 
    make_training_counts(x_location, x_pathogen) |> 
    filter(date > origin_date)
  
  
  draws <- fit |>
    tidybayes::tidy_draws()

  pars <- c(
    "scaling_annual", "scaling_weekly", "scaling_noise", "ls_noise", "ls_annual", "phi",
    "delay_mean", "delay_sd"
  )

  p_hist <- mcmc_hist(draws |> mutate(chain = .chain), pars = pars) +
    theme_bw()

  p_trace <- mcmc_trace(draws |> mutate(chain = .chain), pars = pars) +
    theme_bw()

  p_hist / p_trace

  ggsave(str_c(round_results, "/", x_pathogen, "_", x_location, "_plot_MCMC.png"), width = 10, height = 10)

  p_pairs <- mcmc_pairs(draws |> mutate(chain = .chain), pars = pars, off_diag_args = list(size = 0.3, alpha = 0.5))

  ggsave(str_c(round_results, "/", x_pathogen, "_", x_location, "_plot_MCMC_pairs.png"), p_pairs, width = 16, height = 16)
  
  
  holidays <- get_school_holidays() |> 
    filter(location == x_location)
  
  trajectory_data <- fit |> 
    spread_draws(log_R[ix_x], inf[ix_x], inf_delayed[ix_x], cases[ix_x], f_annual[ix_x], f_noise[ix_x]) |> 
    mutate(date = training_counts$date[ix_x])
  
  p_reff <- plot_reff_components(trajectory_data, origin_date, holidays)
  
  p_latent_inf <- plot_latent_inf(trajectory_data, training_counts, origin_date, holidays) +
    
    geom_line(aes(x = date, y = y_smooth), colour = "orangered",
              eval_counts)
  
  p_reff / p_latent_inf
  
  ggsave(str_c(round_results, "/", x_pathogen, "_", x_location, "_plot_latent_inf.png"), width = 14, height = 8)
  
  
  make_forecast_plot(trajectory_data, training_counts, eval_counts, origin_date, holidays, x_location, x_pathogen)
  
  ggsave(str_c(round_results, "/", x_pathogen, "_", x_location, "_plot_forecast.png"), width = 10, height = 10)
  
}
