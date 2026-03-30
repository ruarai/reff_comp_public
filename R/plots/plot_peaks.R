`
library(tidyverse)
library(tidybayes)
library(bayesplot)
library(patchwork)

source("R/functions.R")


plot_peaks <- function(x_location, x_pathogen, round_id, y_upper) {
  round_results <- fs::path("results", round_id)
  
  fit <- read_rds(str_c(round_results, "/", x_pathogen, "_", x_location, ".rds"))
  training_counts <- arrow::read_parquet(str_c(round_results, "/", x_pathogen, "_", x_location, "_training_counts.parquet"))
  
  
  origin_date <- training_counts |> 
    drop_na(y) |> 
    pull(date) |> 
    max()
  
  eval_counts <- read_csv("data/combined-count-smoothed-2026-01-15.csv") |> 
    group_by(location, pathogen) |> 
    mutate(outlier = FALSE, cases = NA,
           cases_trend_and_remainder = cases_trend) |> 
    filter(location == x_location, pathogen == x_pathogen) |> 
    make_training_counts(x_location, x_pathogen)
  
  
  
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
              latent_cases = max(latent_cases))
  
  
  eval_peak <- eval_counts |> 
    filter(date > ymd("2025-01-01"), date < ymd("2026-01-01")) |> 
    summarise(date = date[y_smooth == max(y_smooth)][1],
              y_smooth = max(y_smooth))
  
  p_peak <- ggplot() +
    
    geom_vline(xintercept = origin_date, linetype = "44") +
    
    geom_point(aes(x = date, y = latent_cases),
               size = 0.5, stroke = 0.5, alpha = 0.5,
               peak_data) +
    
    geom_line(aes(x = date, y = y_smooth),
              colour = "orangered", alpha = 0.7,
              eval_counts) +
    
    geom_point(aes(x = date, y = y_smooth),
               colour = "orangered", size = 3,
               eval_peak) +
    
    xlab("Notification date") + ylab("Count") +
    
    coord_cartesian(ylim = c(0, y_upper),
                    xlim = c(ymd("2024-03-01"), ymd("2026-01-01"))) +
    
    theme_bw() +
    
    theme(strip.text = element_blank(),
          legend.position = "none") +
    
    ggtitle(str_c(x_location, " - ", x_pathogen))
  
  
  
  p_date <- ggplot() +
    geom_histogram(aes(x = date),
                   binwidth = 7,
                   peak_data) +
    
    geom_vline(xintercept = eval_peak$date,
               linewidth = 1.0,
               colour = "orangered") +
    
    geom_vline(xintercept = origin_date, linetype = "44")  +
    
    coord_cartesian(xlim = c(ymd("2024-03-01"), ymd("2026-01-01"))) +
    
    ylab(NULL) + xlab(NULL) +
    
    theme_bw()
  
  p_cases <- ggplot() +
    geom_histogram(aes(y = latent_cases),
                   binwidth = 10,
                   peak_data) +
    
    geom_hline(yintercept = eval_peak$y_smooth,
               linewidth = 1.0,
               colour = "orangered") +
    
    coord_cartesian(ylim = c(0, y_upper)) +
    
    ylab(NULL) + xlab(NULL) +
    
    theme_bw()
  
  
  (p_cases + p_peak + plot_spacer() + p_date) +
    plot_layout(ncol = 2, heights = c(1, 0.2), widths = c(0.2, 1))
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

make_peaks_gif("WA", "RSV", 500)
make_peaks_gif("WA", "flu", 500)
make_peaks_gif("WA", "SARSCOV2", 250)

make_peaks_gif("VIC", "flu", 1000)
make_peaks_gif("VIC", "RSV", 500)
make_peaks_gif("VIC", "SARSCOV2", 500)

make_peaks_gif("NSW", "flu", 2500)
make_peaks_gif("NSW", "RSV", 1000)
make_peaks_gif("NSW", "SARSCOV2", 1000)


plot_peaks("WA", "RSV", ymd("2025-05-08"), y_upper = 500)


