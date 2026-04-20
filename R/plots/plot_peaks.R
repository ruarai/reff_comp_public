
library(tidyverse)
library(tidybayes)
library(bayesplot)
library(patchwork)

source("R/functions.R")


plot_peaks <- function(x_location, x_pathogen, round_id, y_upper) {
  round_results <- fs::path("saved_results", "results_baseline", round_id)
  
  fit <- read_rds(str_c(round_results, "/", x_pathogen, "_", x_location, ".rds"))
  training_counts <- arrow::read_parquet(str_c(round_results, "/", x_pathogen, "_", x_location, "_training_counts.parquet"))
  
  
  peak_date_range <- c(floor_date(ymd(round_id), "year"), floor_date(ymd(round_id), "year") + years(1))
  plot_date_range <- c(floor_date(ymd(round_id), "year") - years(1), floor_date(ymd(round_id), "year") + years(1))
  
  
  origin_date <- training_counts |> 
    drop_na(y) |> 
    pull(date) |> 
    max()
  
  eval_counts <- read_csv("data/combined-count-smoothed-2026-04-15.csv", show_col_types = FALSE) |> 
    group_by(location, pathogen) |> 
    filter(location == x_location, pathogen == x_pathogen, target == "case_incidence") |> 
    rename(date = notification_date)
  
  eval_peak <- eval_counts |> 
    filter(date > ymd("2025-01-01"), date < ymd("2026-01-01")) |> 
    summarise(date = date[value_trend == max(value_trend)][1],
              peak_case_incidence = max(value_trend),
              .groups = "drop")
    
  
  
  school_holidays <- get_school_holidays() |> 
    filter(location == x_location)
  
  public_holidays <- get_public_holidays() |> 
    filter(location == x_location)
  
  trajectory_data <- fit |> 
    spread_draws(latent_cases_log[ix_t], inf_delayed[ix_t]) |> 
    mutate(date = training_counts$date[ix_t],
           latent_cases = exp(latent_cases_log))
  
  
  
  peak_data <- trajectory_data |> 
    filter(date > ymd("2025-01-01"), date < ymd("2026-01-01")) |> 
    group_by(.draw) |> 
    summarise(date = date[latent_cases == max(latent_cases)][1],
              peak_case_incidence = max(latent_cases),
              .groups = "drop")
  
  upper_cases <- quantile(peak_data$peak_case_incidence, 0.9)
  
  plot_data <- peak_data |> 
    
    filter(peak_case_incidence < upper_cases)
  
  p_scatter <- ggplot() +
    geom_density2d_filled(aes(x = date, y = peak_case_incidence),
                          n = 512, bins = 6,
                          plot_data) +
    
    geom_point(aes(x = date, y = peak_case_incidence),
               peak_data |> slice_sample(n = 500),
               size = 0.1, stroke = 0.5, alpha = 0.4) +
    
    geom_line(aes(x = date, y = value_trend_and_remainder),
              alpha = 0.5,
              eval_counts) +
    
    geom_vline(xintercept = origin_date, linetype = "44") +
    
    geom_point(aes(x = date, y = peak_case_incidence),
               eval_peak, size = 2, colour = "orangered") +
    
    scale_y_log10() +
    scale_x_date(date_breaks = "month") +
    
    xlab(NULL) + ylab("Case incidence") +
    
    coord_cartesian(xlim = plot_date_range, ylim = c(100, y_upper)) +
    
    # scale_fill_brewer(palette = "Blues") +
    scale_fill_manual(values = c("transparent", RColorBrewer::brewer.pal(5, "Blues"))) +
    
    theme_bw() +
    
    ggtitle(glue::glue("{x_location}-{x_pathogen}")) +
    
    theme(legend.position = "none",
          axis.text.x = element_blank())
  
  p_scatter
  
  p_date <- ggplot(peak_data) +
    geom_histogram(aes(x = date, y = after_stat(ndensity)), binwidth = 7) +
    # stat_ecdf(aes(x = date)) +
    
    scale_x_date(date_breaks = "month", labels = scales::label_date_short()) +
    
    geom_vline(xintercept = origin_date, linetype = "44") +
    
    geom_vline(aes(xintercept = date),
               eval_peak, linetype = "44", colour = "orangered") +
    
    xlab(NULL) + ylab(NULL) +
    
    coord_cartesian(xlim = plot_date_range) +
    
    theme_bw()
  
  p_cases <- ggplot() +
    geom_histogram(aes(y = peak_case_incidence, x = after_stat(ndensity)), 
                   breaks = seq(0, y_upper, length.out = 30), peak_data) +
    # stat_ecdf(aes(y = peak_case_incidence), peak_data) +
    
    scale_y_log10() +
    geom_hline(aes(yintercept = peak_case_incidence),
               eval_peak, linetype = "44", colour = "orangered") +
    
    coord_cartesian(ylim = c(100, y_upper)) +
    
    xlab(NULL) + ylab(NULL) +
    
    theme_bw()
  
  p <- wrap_plots(
    p_scatter, p_cases,
    p_date, plot_spacer(),
    ncol = 2,
    widths = c(0.75, 0.25), heights = c(0.75, 0.25)
  )
  
  p
}


make_peaks_gif <- function(x_location, x_pathogen, y_upper) {
  rounds <- ymd("2025-05-08") + weeks(0:14)
  
  plots <- map(
    rounds,
    function(round_id) {
      plot_peaks(x_location, x_pathogen, round_id, y_upper = y_upper)
    }
  )
  
  
  out_file <- fs::path("plots", "peak_animations", str_c(x_location, "_", x_pathogen, ".gif"))
  
  plots_to_gif(plots, out_file, 12, 7)
}

plot_peaks("WA", "RSV", "2025-08-14", y_upper = 250)

make_peaks_gif("WA", "RSV", 250)
make_peaks_gif("WA", "flu", 500)
make_peaks_gif("WA", "SARSCOV2", 250)

plot_peaks("VIC", "flu", "2025-08-14", y_upper = 1500)

make_peaks_gif("VIC", "flu", 1500)
make_peaks_gif("VIC", "RSV", 500)
make_peaks_gif("VIC", "SARSCOV2", 500)

make_peaks_gif("NSW", "flu", 2500)
make_peaks_gif("NSW", "RSV", 1000)
make_peaks_gif("NSW", "SARSCOV2", 1000)


plot_peaks("WA", "RSV", ymd("2025-05-08"), y_upper = 500)


