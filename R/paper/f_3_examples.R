source("R/functions.R")
source("R/paper/common.R")

library(pammtools)


eval_counts <- read_csv("data/combined-count-smoothed-2026-04-15.csv", show_col_types = FALSE)

round_id <- "2025-07-10"
x_pathogen <- "SARSCOV2"
x_location <- "VIC"

plot_forecast <- function(x_pathogen, x_location, round_id) {
  x_target <- "case_incidence"
  
  job_name <- str_c(x_pathogen, "_", x_location, "_", x_target)
  
  fit <- read_rds(str_c("saved_results/baseline/", round_id, "/", job_name, ".rds"))
  training_counts <- read_parquet(str_c("saved_results/baseline/", round_id, "/", job_name, "_training_counts.parquet"))
  
  origin_date <- training_counts$origin_date[[1]] - days(1)
  
  trajectory_data <- fit |> 
    spread_draws(c(log_R, inf, inf_delayed, cases, f_annual,
                   f_noise, latent_cases_log)[ix_t],
                 log_R_0) |> 
    mutate(date = training_counts$date[ix_t],
           latent_cases = exp(latent_cases_log),
           holiday_effect = log_R - log_R_0 - f_annual - f_noise)
  
  plot_data <- trajectory_data |> 
    filter(date > origin_date - days(7 * 8),
           date <= origin_date + days(7 * 4))
  
  plot_data_forecast <- trajectory_data |> 
    filter(date > origin_date,
           date <= origin_date + days(7 * 4))
  
  plot_data_eval <- eval_counts |> 
    filter(location == x_location, pathogen == x_pathogen) |> 
    filter(notification_date > origin_date,
           notification_date <= origin_date + days(7 * 4))
  
  plot_data_training <- training_counts |> 
    filter(date >= origin_date - days(7 * 8)) |> 
    drop_na(y)
  
  
  school_holidays <- get_school_holidays() |> 
    filter(location == x_location,
           end_date > origin_date - days(7 * 8),
           start_date <= origin_date + days(7 * 4))
  
  public_holidays <- get_public_holidays() |> 
    filter(location == x_location,
           date > origin_date - days(7 * 8),
           date <= origin_date + days(7 * 4))
  
  pathogen_label <- c("RSV" = "RSV", "flu" = "Influenza", "SARSCOV2" = "SARS-CoV-2")[x_pathogen]
  
  p_cases <- ggplot(plot_data_forecast) +
    
    geom_rect(aes(xmin = date - 0.5, xmax = date + 0.5, ymin = -Inf, ymax = Inf),
              alpha = 0.2, fill = "#e31a1c",
              public_holidays) +
    geom_rect(aes(xmin = start_date - 0.5, xmax = end_date + 0.5, ymin = -Inf, ymax = Inf),
              school_holidays, alpha = 0.1) +
    fastribbon::stat_fastribbon(aes(x = date, y = cases, alpha = after_stat(level)),
                                fill = colour_fc,
                                linewidth = 0.0) +
    
    scale_alpha_manual(values = c(0.2, 0.5, 1)) +
    
    geom_vline(xintercept = origin_date,
               linetype = "44") +
    
    geom_line(aes(x = date, y = y), plot_data_training) +
    geom_point(aes(x = date, y = y), plot_data_training) +
    
    geom_line(aes(x = notification_date, y = value),
              plot_data_eval) +
    geom_point(aes(x = notification_date, y = value),
               plot_data_eval) +
    
    scale_fill_brewer() +
    xlab(NULL) + ylab("Count") +
    
    coord_cartesian(xlim = c(origin_date - days(7 * 8), origin_date + days(7 * 4))) +
    
    plot_theme_paper +
    
    theme(legend.position = "none")
  
  p_cases
}


plots <- list(
  plot_forecast("RSV", "QLD", "2025-04-04") +
    ggtitle(NULL, "<b>i</b>. RSV – QLD – Round 2025-04-04"),
  # Good result (closely matches)
  
  plot_forecast("flu", "VIC", "2025-06-05") +
    ggtitle(NULL, "<b>ii</b>. Influenza – VIC – Round 2025-06-05"),
  
  plot_forecast("SARSCOV2", "ACT", "2025-04-24") +
    ggtitle(NULL, "<b>iii</b>. SARS-CoV-2 – ACT – Round 2025-04-24"),
  # Bad result (increase not predicted, although holidays perform ok)
  
  plot_forecast("flu", "NSW", "2025-07-03") +
    ggtitle(NULL, "<b>iv</b>. Influenza – NSW – Round 2025-07-03")
  # Bad result (f_ann pushes down too early)
)


p_examples <- (
  wrap_plots(plots, ncol = 2) /
  plot_spacer() +
    plot_layout(heights = c(1, 0.1)))+ 
  plot_annotation("<b>A</b>. Example high- and low-skill forecasts",
                  theme = plot_theme_annotation)




ggsave(
  "paper_figures/f_3_examples.png",
  p_examples,
  width = 10, height = 6,
  bg = "white"
)

library(magick)

image_composite(
  image = image_read("paper_figures/f_3_examples.png"),
  composite_image = image_read("paper_figures/f_3_examples_overlay.png"),
  operator = "over"
) |> 
  image_write(path = "paper_figures/f_3_examples.png", format = "png")






