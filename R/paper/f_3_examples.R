source("R/functions.R")
source("R/paper/common.R")

library(pammtools)


eval_counts <- read_csv("data/combined-count-smoothed_final_anonymised.csv", show_col_types = FALSE)


gp_draws <- read_parquet("data/paper/gp_draws.parquet")

round_id <- "2025-04-04"
x_pathogen <- "RSV"
x_location <- "QLD"
title <- NULL

plot_forecast <- function(x_pathogen, x_location, round_id, title) {
  # Draws are read from data/paper/plot_data/ (or extracted from the raw fit if that file is absent)
  trajectory_data <- load_fit_trajectories("f_3", "baseline", round_id, x_pathogen, x_location)

  origin_date <- trajectory_data$origin_date[[1]]

  # Observed counts (before and after the forecast date) are from the anonymised final dataset
  training_counts <- read_final_counts(x_pathogen, x_location) |>
    filter(date <= origin_date)

  plot_data <- trajectory_data |> 
    filter(date > origin_date - days(7 * 8),
           date <= origin_date + days(7 * 4))
  
  plot_data_backcast <- trajectory_data |> 
    filter(date > origin_date - days(7 * 8),
           date <= origin_date)
  
  plot_data_forecast <- trajectory_data |> 
    filter(date > origin_date,
           date <= origin_date + days(7 * 4))
  
  plot_data_eval <- eval_counts |> 
    filter(location == x_location, pathogen == x_pathogen, target == "case_incidence") |>
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
  
  plots_common <- list(
    geom_vline(xintercept = origin_date,
               linetype = "44"),
    xlab(NULL),
    plot_theme_paper,
    theme(legend.position = "none", axis.text.x = element_blank(),
          panel.grid.major.x = element_gridline),
    scale_alpha_manual(values = c(0.2, 0.5, 1))
  )
  
  p_cases <- ggplot(plot_data_forecast) +
    
    geom_rect(aes(xmin = date - 0.5, xmax = date + 0.5, ymin = -Inf, ymax = Inf),
              alpha = 0.2, fill = "#e31a1c",
              public_holidays) +
    geom_rect(aes(xmin = start_date - 0.5, xmax = end_date + 0.5, ymin = -Inf, ymax = Inf),
              school_holidays, alpha = 0.1) +
    fastribbon::stat_fastribbon(aes(x = date, y = cases, alpha = after_stat(level)),
                                fill = colour_fc,
                                linewidth = 0.0) +
    
    fastribbon::stat_fastribbon(aes(x = date, y = cases, alpha = after_stat(level)),
                                plot_data_backcast, fill = "#3182BD",
                                linewidth = 0.0) +
    
    geom_line(aes(x = date, y = y), plot_data_training) +
    geom_point(aes(x = date, y = y), plot_data_training, size = 1.3) +
    geom_point(aes(x = date, y = y), plot_data_training, size = 0.5, stroke = 0.1, colour = "white") +
    
    geom_line(aes(x = notification_date, y = value),
              plot_data_eval) +
    geom_point(aes(x = notification_date, y = value),
               plot_data_eval, size = 1.3) +
    geom_point(aes(x = notification_date, y = value),
               plot_data_eval, size = 0.5, stroke = 0.1, colour = "white") +
    
    ylab("Count") +
    
    coord_cartesian(xlim = c(origin_date - days(7 * 8), origin_date + days(7 * 4))) +
    
    ggtitle(title, "Case incidence") +
    
    plots_common
  
  gp_eval <- gp_draws |> 
    filter(location == x_location, pathogen == x_pathogen,
           date > origin_date - days(7 * 8),
           date < origin_date + days(7 * 4))
  
  p_Rt <- ggplot(plot_data_forecast) +
    
    geom_rect(aes(xmin = date - 0.5, xmax = date + 0.5, ymin = -Inf, ymax = Inf),
              alpha = 0.2, fill = "#e31a1c",
              public_holidays) +
    geom_rect(aes(xmin = start_date - 0.5, xmax = end_date + 0.5, ymin = -Inf, ymax = Inf),
              school_holidays, alpha = 0.1) +
    fastribbon::stat_fastribbon(aes(x = date - 0.5, y = exp(log_R), alpha = after_stat(level)),
                                fill = colour_fc, geom = "stepribbon",
                                linewidth = 0.0) +
    fastribbon::stat_fastribbon(aes(x = date - 0.5, y = exp(log_R), alpha = after_stat(level)),
                                plot_data_backcast, fill = "#3182BD", geom = "stepribbon",
                                linewidth = 0.0) +
    fastribbon::stat_fastribbon(aes(x = date, y = exp(log_R)),
                                gp_eval, .width = 0.5, geom = "linerange",
                                linewidth = 0.7, colour = "black", alpha = 0.8) +
    
    geom_hline(yintercept = 1, linetype = "11", linewidth = 0.7, alpha = 0.5) +
    
    scale_y_continuous(breaks = c(0.75, 1, 1.25)) +
    
    ylab("<i>R</i>(<i>t</i>)") +
    
    coord_cartesian(ylim = c(0.6, 1.4),
                    xlim = c(origin_date - days(7 * 8), origin_date + days(7 * 4))) +
    
    ggtitle(NULL, "Effective reproduction number <i>R</i>(<i>t</i>)") +
    
    plots_common
  
  p_ann <- ggplot(plot_data_forecast) +

    fastribbon::stat_fastribbon(aes(x = date - 0.5, y = exp(f_annual), alpha = after_stat(level)),
                                fill = colour_fc, geom = "stepribbon",
                                linewidth = 0.0) +
    fastribbon::stat_fastribbon(aes(x = date - 0.5, y = exp(f_annual), alpha = after_stat(level)),
                                plot_data_backcast, fill = "#3182BD", geom = "stepribbon",
                                linewidth = 0.0) +
    fastribbon::stat_fastribbon(aes(x = date, y = exp(f_annual)),
                                gp_eval, .width = 0.5, geom = "linerange",
                                linewidth = 0.7, colour = "black", alpha = 0.8) +

    geom_hline(yintercept = 1, linetype = "11", linewidth = 0.7, alpha = 0.5) +
    
    scale_y_log10(breaks = c(0.8, 1, 1.25), labels = c("-20%", "0%", "+25%")) +

    ylab("Effect") +

    coord_cartesian(ylim = c(0.75, 1.3),
                    xlim = c(origin_date - days(7 * 8), origin_date + days(7 * 4))) +

    ggtitle(NULL, "Effect of seasonal GP <i>f</i><sub>seasonal</sub>(<i>t</i>)") +
    
    plots_common
  
  p_noise <- ggplot(plot_data_forecast) +
    
    fastribbon::stat_fastribbon(aes(x = date - 0.5, y = exp(f_noise), alpha = after_stat(level)),
                                fill = colour_fc, geom = "stepribbon",
                                linewidth = 0.0) +
    fastribbon::stat_fastribbon(aes(x = date - 0.5, y = exp(f_noise), alpha = after_stat(level)),
                                plot_data_backcast, fill = "#3182BD", geom = "stepribbon",
                                linewidth = 0.0) +
    fastribbon::stat_fastribbon(aes(x = date, y = exp(f_noise)),
                                gp_eval, .width = 0.5, geom = "linerange",
                                linewidth = 0.7, colour = "black", alpha = 0.8) +
    
    geom_hline(yintercept = 1, linetype = "11", linewidth = 0.7, alpha = 0.5) +
    
    scale_y_log10(breaks = c(0.8, 1, 1.25), labels = c("-20%", "0%", "+25%")) +
    
    scale_x_date(labels = scales::label_date_short()) +
    
    ylab("Effect") +
    
    coord_cartesian(ylim = c(0.75, 1.3),
                    xlim = c(origin_date - days(7 * 8), origin_date + days(7 * 4))) +
    
    ggtitle(NULL, "Effect of noise GP <i>f</i><sub>noise</sub>(<i>t</i>)") +
    
    plots_common +
    
    theme(axis.text.x = element_text(colour = "black"))
  
  p_cases / p_Rt / p_ann / p_noise +
    plot_layout(heights = c(1.5, 1, 1, 1))
}


plots <- list(
  plot_forecast("RSV", "QLD", "2025-04-04", "<b>i</b>. RSV – Queensland – Round 2025-04-04<br>(exemplar good forecasting performance)"),
  # Good result (closely matches)
  
  # plot_forecast("flu", "VIC", "2025-06-05", "<b>ii</b>. Influenza – VIC – Round 2025-06-05"),
  
  # plot_forecast("SARSCOV2", "ACT", "2025-04-24", "<b>ii</b>. SARS-CoV-2 – ACT – Round 2025-04-24")
  # Bad result (increase not predicted, although holidays perform ok)
  
  # plot_forecast("flu", "NSW", "2025-07-03", "<b>ii</b>. Influenza – NSW – Round 2025-07-03")
  # Bad result (f_ann pushes down too early)
  
  # plot_forecast("RSV", "NSW", "2025-03-14", "<b>i</b>. RSV – NSW – Round 2025-03-14"),
  
  # plot_forecast("RSV", "flu", "2025-08-07", "<b>i</b>. RSV – Influenza – Round 2025-08-07")
  plot_forecast("SARSCOV2", "VIC", "2025-05-08", "<b>ii</b>. SARS-CoV-2 – Victoria – Round 2025-05-08<br>(exemplar poor forecasting performance)")
)


p_examples <- (
  wrap_plots(plots, ncol = 2) /
  plot_spacer() +
    plot_layout(heights = c(1, 0.1)))



ggsave(
  "paper_figures/f_3_examples.png",
  p_examples,
  width = 10, height = 8,
  bg = "white"
)

library(magick)

image_composite(
  image = image_read("paper_figures/f_3_examples.png"),
  composite_image = image_read("paper_figures/f_3_examples_overlay.png"),
  operator = "over"
) |> 
  image_write(path = "paper_figures/f_3_examples.png", format = "png")






