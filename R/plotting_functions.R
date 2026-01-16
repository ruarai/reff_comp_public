
plot_reff_components <- function(trajectory_data, origin_date, school_holidays, public_holidays) {
  
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
              school_holidays, alpha = 0.1) +
    
    geom_rug(aes(x = date), public_holidays) +
    
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
    
    coord_cartesian(ylim = c(1/2, 2)) +
    
    geom_vline(xintercept = origin_date, linetype = "44") +
    
    theme_bw() +
    theme(legend.position = "none")
}

plot_latent_inf <- function(trajectory_data, training_counts, origin_date, school_holidays, public_holidays) {
  plot_quants <- make_intervals(trajectory_data, "inf")
  plot_quants_delayed <- make_intervals(trajectory_data, "inf_delayed")
  
  
  alpha_vals <- scales::rescale(rev(1/1.7^(1:8)), to = c(0.05, 0.99))
  
  
  ggplot() +
    
    geom_vline(aes(xintercept = date), tibble(date = seq(ymd("2023-01-01"), ymd("2027-01-01"), "1 year")),
               alpha = 0.4) +
    
    geom_rect(aes(xmin = start_date, xmax = end_date, ymin = 0, ymax = Inf),
              school_holidays, alpha = 0.1) +
    
    geom_rug(aes(x = date), public_holidays) +
    
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

plot_cases <- function(trajectory_data, training_counts, origin_date, school_holidays, public_holidays) {
  plot_quants <- make_intervals(trajectory_data, "cases")
  
  
  alpha_vals <- scales::rescale(rev(1/1.7^(1:8)), to = c(0.05, 0.99))
  
  
  ggplot(plot_quants) +
    
    geom_rect(aes(xmin = start_date, xmax = end_date, ymin = 0, ymax = Inf),
              school_holidays, alpha = 0.1) +
    
    geom_rug(aes(x = date), public_holidays) +
    
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

make_forecast_plot <- function(trajectory_data, training_counts, eval_counts, origin_date, school_holidays, public_holidays, x_location, x_pathogen) {
  
  forecast_plot_window <- c(origin_date - days(7 * 4), origin_date + days(7 * 12))
  
  
  p_forecast_components <- plot_reff_components(trajectory_data, origin_date, school_holidays, public_holidays) +
    
    xlim(forecast_plot_window) +
    
    geom_vline(xintercept = origin_date + days(4 * 7), linetype = "24")
  
  p_forecast_latent <- plot_latent_inf(trajectory_data, training_counts, origin_date, school_holidays, public_holidays) +
    
    xlim(forecast_plot_window) +
    
    geom_line(aes(x = date, y = y_smooth), colour = "orangered",
              eval_counts) +
    
    geom_vline(xintercept = origin_date + days(4 * 7), linetype = "24")
  
  p_forecast_cases <- plot_cases(trajectory_data, training_counts, origin_date, school_holidays, public_holidays) +
    
    xlim(forecast_plot_window) +
    
    geom_line(aes(x = date, y = y), colour = "orangered",
              eval_counts) +
    
    geom_vline(xintercept = origin_date + days(4 * 7), linetype = "24") +
    
    labs(caption = str_c(x_location, " - ", x_pathogen))
  
  p_forecast_components / p_forecast_latent / p_forecast_cases
}

plot_holiday_effects <- function(fit, public_holidays) {
  school_holiday_draws <- fit |> 
    spread_draws(school_holiday_asc_effect[school_holiday_type],
                 school_holiday_Rt_effect[school_holiday_type]) |> 
    mutate(school_holiday_type = school_holiday_type - 1) |> 
    filter(school_holiday_type > 0) |> 
    mutate(holiday_name = case_when(
      school_holiday_type == 1 ~ "Summer",
      school_holiday_type == 2 ~ "Autumn",
      school_holiday_type == 3 ~ "Winter",
      school_holiday_type == 4 ~ "Spring"
    )) |> 
    ungroup() |> 
    mutate(holiday_name = factor(holiday_name, c("Summer", "Autumn", "Winter", "Spring")))
  
  
  public_holiday_labels <- public_holidays |> 
    distinct(holiday_name) |> 
    mutate(public_holiday_type = row_number())
  
  public_holiday_draws <- fit |> 
    spread_draws(public_holiday_asc_effect[public_holiday_type],
                 public_holiday_Rt_effect[public_holiday_type]) |> 
    mutate(public_holiday_type = public_holiday_type - 1) |> 
    filter(public_holiday_type > 0) |> 
    
    left_join(public_holiday_labels) |> 
    
    ungroup() |> 
    mutate(holiday_name = factor(holiday_name, public_holiday_labels$holiday_name))
  
  
  p_school_Rt <- ggplot(school_holiday_draws) +
    
    geom_vline(xintercept = 1, linetype = "22") +
    ggdist::stat_interval(aes(y = holiday_name, x = exp(school_holiday_Rt_effect)),
                          .width = c(0.5, 0.9, 0.95),
                          orientation = "horizontal") +
    
    geom_blank(
      aes(
        xmin = exp(quantile(public_holiday_draws$public_holiday_Rt_effect, 0.025)),
        xmax = exp(quantile(public_holiday_draws$public_holiday_Rt_effect, 0.975)),
      )
    ) +
    
    scale_colour_brewer() +
    
    xlab("Rt effect") + ylab(NULL) +
    scale_x_log10(breaks = scales::log_breaks(n = 7)) +
    
    theme_bw() +
    theme(legend.position = "none")
  
  
  p_school_asc <- ggplot(school_holiday_draws) +
    
    geom_vline(xintercept = 1, linetype = "22") +
    ggdist::stat_interval(aes(y = holiday_name, x = exp(school_holiday_asc_effect)),
                          .width = c(0.5, 0.9, 0.95),
                          orientation = "horizontal") +
    
    geom_blank(
      aes(
        xmin = exp(quantile(public_holiday_draws$public_holiday_asc_effect, 0.025)),
        xmax = exp(quantile(public_holiday_draws$public_holiday_asc_effect, 0.975)),
      )
    ) +
    
    scale_colour_brewer() +
    
    xlab("Ascertainment effect") + ylab(NULL) +
    scale_x_log10(breaks = scales::log_breaks(n = 7)) +
    
    theme_bw() +
    theme(legend.position = "none")
  
  
  p_public_Rt <- ggplot(public_holiday_draws) +
    
    geom_vline(xintercept = 1, linetype = "22") +
    ggdist::stat_interval(aes(y = holiday_name, x = exp(public_holiday_Rt_effect)),
                          .width = c(0.5, 0.9, 0.95),
                          orientation = "horizontal") +
    
    geom_blank(
      aes(
        xmin = exp(quantile(school_holiday_draws$school_holiday_Rt_effect, 0.025)),
        xmax = exp(quantile(school_holiday_draws$school_holiday_Rt_effect, 0.975)),
      )
    ) +
    
    scale_colour_brewer() +
    
    xlab("Rt effect") + ylab(NULL) +
    scale_x_log10(breaks = scales::log_breaks(n = 7)) +
    
    theme_bw() +
    theme(legend.position = "none")
  
  p_public_asc <- ggplot(public_holiday_draws) +
    
    geom_vline(xintercept = 1, linetype = "22") +
    ggdist::stat_interval(aes(y = holiday_name, x = exp(public_holiday_asc_effect)),
                          .width = c(0.5, 0.9, 0.95),
                          orientation = "horizontal") +
    
    geom_blank(
      aes(
        xmin = exp(quantile(school_holiday_draws$school_holiday_asc_effect, 0.025)),
        xmax = exp(quantile(school_holiday_draws$school_holiday_asc_effect, 0.975)),
      )
    ) +
    
    scale_colour_brewer() +
    
    xlab("Ascertainment effect") + ylab(NULL) +
    scale_x_log10(breaks = scales::log_breaks(n = 7)) +
    
    theme_bw() +
    theme(legend.position = "none")
  
  (p_school_Rt / p_public_Rt) + plot_layout(heights = c(1, 2)) |
    (p_school_asc / p_public_asc) + plot_layout(heights = c(1, 2))
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
    "scaling_annual", "scaling_weekly_Rt", "scaling_weekly_asc", "scaling_noise", "ls_noise", "ls_annual", "phi",
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
  
  mcmc_pairs(draws |> mutate(chain = .chain), pars = c("weekend_effect_Rt", "obs_offset_weekly[1]", "obs_offset_weekly[4]"), off_diag_args = list(size = 0.3, alpha = 0.5))
  
  school_holidays <- get_school_holidays() |> 
    filter(location == x_location)
  
  public_holidays <- get_public_holidays() |> 
    filter(location == x_location)
  
  trajectory_data <- fit |> 
    spread_draws(log_R[ix_x], inf[ix_x], inf_delayed[ix_x], cases[ix_x], f_annual[ix_x], f_noise[ix_x]) |> 
    mutate(date = training_counts$date[ix_x])
  
  p_reff <- plot_reff_components(trajectory_data, origin_date, school_holidays, public_holidays) +
    xlim(c(ymd(min(trajectory_data$date)), ymd(origin_date + days(365))))
  
  p_latent_inf <- plot_latent_inf(trajectory_data, training_counts, origin_date, school_holidays, public_holidays) +
    
    geom_line(aes(x = date, y = y_smooth), colour = "orangered",
              eval_counts) +
    xlim(c(ymd(min(trajectory_data$date)), ymd(origin_date + days(365))))
  
  p_reff / p_latent_inf
  
  ggsave(str_c(round_results, "/", x_pathogen, "_", x_location, "_plot_latent_inf.png"), width = 14, height = 8)
  
  
  make_forecast_plot(
    trajectory_data, training_counts, eval_counts, origin_date,
    school_holidays, public_holidays, x_location, x_pathogen
  )
  
  ggsave(str_c(round_results, "/", x_pathogen, "_", x_location, "_plot_forecast.png"), width = 10, height = 10)
  
  plot_holiday_effects(fit, public_holidays)
  
  ggsave(str_c(round_results, "/", x_pathogen, "_", x_location, "_holiday_effects.png"), width = 12, height = 8)
  
}
