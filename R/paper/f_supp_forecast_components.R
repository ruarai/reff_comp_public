source("R/functions.R")
source("R/paper/common.R")

library(pammtools)

eval_counts <- read_csv("data/combined-count-smoothed_final_anonymised.csv", show_col_types = FALSE)

plot_forecast <- function(x_pathogen, x_location, round_id, title, model_version, y_lim) {
  # Draws are read from data/paper/plot_data/ (or extracted from the raw fit if that file is absent)
  trajectory_data <- load_fit_trajectories(
    "f_supp_components", model_version, round_id, x_pathogen, x_location,
    vars = c("log_R", "cases")
  )

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
    
    coord_cartesian(xlim = c(origin_date - days(7 * 8), origin_date + days(7 * 4)),
                    ylim = c(0, y_lim)) +
    
    ggtitle(title, "Case incidence") +
    
    plots_common
  
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
    
    geom_hline(yintercept = 1, linetype = "11", linewidth = 0.7, alpha = 0.5) +
    
    scale_y_continuous(breaks = c(0.75, 1, 1.25)) +
    
    ylab("<i>R</i>(<i>t</i>)") +
    
    coord_cartesian(ylim = c(0.6, 1.4),
                    xlim = c(origin_date - days(7 * 8), origin_date + days(7 * 4))) +
    
    ggtitle(NULL, "Effective reproduction number <i>R</i>(<i>t</i>)") +
    
    plots_common
  
  p_cases / p_Rt +
    plot_layout(heights = c(1.5, 1))
}


p_baseline_top <- plot_forecast("RSV", "QLD", "2025-04-04", "With holiday and seasonality effects", "baseline", y_lim = 400)
p_no_holidays_top <- plot_forecast("RSV", "QLD", "2025-04-04", "Without holiday or seasonality effects", "no_holidays_no_annual", y_lim = 400)

p_top <- (p_baseline_top | p_no_holidays_top) + 
  plot_annotation("<b>A</b>. Forecast of RSV case incidence in QLD (round 2025-04-04)",
                  theme = plot_theme_annotation) 


p_baseline_bottom <- plot_forecast("flu", "VIC", "2025-06-26", "With holiday and seasonality effects", "baseline", y_lim = 2500)
p_no_holidays_bottom <- plot_forecast("flu", "VIC", "2025-06-26", "Without holiday or seasonality effects", "no_holidays_no_annual", y_lim = 2500)

p_bottom <- (p_baseline_bottom | p_no_holidays_bottom) + 
  plot_annotation("<b>B</b>. Forecast of influenza case incidence in VIC (round 2025-06-26)",
                  theme = plot_theme_annotation) 

p_full <- wrap_elements(full = p_top) / wrap_elements(full = p_bottom)


ggsave(
  "paper_figures/supp/f_forecast_components.png",
  p_full,
  width = 10, height = 10,
  bg = "white"
)


