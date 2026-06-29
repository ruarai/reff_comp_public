
source("R/functions.R")
source("R/paper/common.R")

baseline_quants <- read_parquet("saved_results/baseline_quants.parquet")
combined_counts <- read_csv("data/combined-count-smoothed_paper.csv")


x_location <- "NSW"
x_pathogen <- "RSV"

n_weeks_ahead <- 4

plot_forecasts <- function(x_location, x_pathogen) {
  
  loc_quants <- baseline_quants |> 
    filter(location == x_location, pathogen == x_pathogen,
           horizon > 0, horizon <= n_weeks_ahead * 7,
           origin_date > ymd("2025-01-01")) |> 
    mutate(date = origin_date + horizon, .before = 4) |> 
    ungroup() |> 
    
    arrange(round_id) |> 
    mutate(ix_round = match(round_id, unique(round_id)),
           ix_row = (ix_round - 1) %% n_weeks_ahead)
  
  
  origin_dates <- loc_quants |> 
    distinct(origin_date, ix_row)
  
  labels <- loc_quants |> 
    distinct(origin_date, ix_row, ix_round)
  
  plot_counts <- combined_counts |> 
    filter(location == x_location, pathogen == x_pathogen,
           notification_date > min(loc_quants$date),
           notification_date < max(loc_quants$date))
  
  school_holidays <- get_school_holidays() |> 
    filter(location == x_location)
  
  public_holidays <- get_public_holidays() |> 
    filter(location == x_location)
  
  y_lim <- max(50, max(plot_counts$value_trend) * 3.5)
  
  
  ggplot() +
    geom_rect(aes(xmin = date - 0.5, xmax = date + 0.5, ymin = -Inf, ymax = Inf),
              alpha = 0.2, fill = "#e31a1c",
              public_holidays) +
    
    geom_rect(aes(xmin = start_date - 0.5, xmax = end_date + 0.5, ymin = -Inf, ymax = Inf),
              school_holidays, alpha = 0.1) +
    geom_ribbon(aes(x = date, ymin = lower_95, ymax = upper_95, group = round_id),
                loc_quants, fill = "#D6E6F2") +
    geom_ribbon(aes(x = date, ymin = lower_80, ymax = upper_80, group = round_id),
                loc_quants, fill = "#84B4D8") +
    geom_ribbon(aes(x = date, ymin = lower_50, ymax = upper_50, group = round_id),
                loc_quants, fill = "#3182BD") +
    
    geom_line(aes(x = notification_date, y = value),
              plot_counts) +
    
    geom_vline(aes(xintercept = origin_date + 1), origin_dates, linetype = "44") +
    
    geom_text(aes(x = origin_date + days(3), y = y_lim * 0.9, label = ix_round),
              labels, hjust = 0, size = 3) +
    
    coord_cartesian(ylim = c(0, y_lim),
                    xlim = c(min(loc_quants$date), max(loc_quants$date))) +
    
    ylab("Value (trend and remainder)") + xlab(NULL) +
    
    facet_wrap(~ix_row, ncol = 1) +
    
    plot_theme_paper +
    
    theme(strip.text.x = element_blank())
}

plot_forecasts("NT", "flu")

sort(locations_paper) |> 
  map(function(location) {
    
    p_RSV <- plot_forecasts(location, "RSV") +
      ggtitle(NULL, "<b>A</b>. RSV")
    p_flu <- plot_forecasts(location, "flu") +
      ggtitle(NULL, "<b>B</b>. Influenza")
    p_SARSCOV2 <- plot_forecasts(location, "SARSCOV2") +
      ggtitle(NULL, "<b>C</b>. SARS-CoV-2")
    
    p <- p_RSV /
      p_flu / 
      p_SARSCOV2
    
    ggsave(
      glue::glue("paper_figures/supp/f_forecasts_{location}.png"),
      p + 
        plot_annotation(glue::glue("All case incidence forecasts ({location})"),
                        theme = plot_theme_annotation),
      width = 10, height = 10
    )
    
  })

