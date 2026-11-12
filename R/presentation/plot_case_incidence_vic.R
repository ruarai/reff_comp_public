source("R/functions.R")
source("R/paper/common.R")

combined_counts <- read_csv("data/combined-count-smoothed_final_anonymised.csv")

combined_counts <- read_csv("../forecast-hub/processed-data/2026/2026-08-19/combined-count-smoothed-2026-08-19.csv")

training_dates <- get_training_dates()

plot_data <- combined_counts |> 
  label_pathogen() |> 
  left_join(training_dates) |> 
  # filter(notification_date >= cutoff_date, notification_date < ymd("2026-01-01")) |> 
  filter(location == "VIC", pathogen == "RSV")



plot_pathogen <- function(x_pathogen) {
  
  y_lims <- plot_data |> 
    group_by(pathogen, location) |> 
    summarise(ylim_trend = max(value_trend) * 1.5) |> 
    filter(pathogen == x_pathogen)
  
  
  scales <- map(1:nrow(y_lims),
                function(i) {
                  scale_y_continuous(limits = c(0, y_lims$ylim_trend[[i]]))
                }) |> 
    `names<-`(y_lims$location)
  
  
  ggplot() +
    
    geom_line(aes(x = notification_date, y = value), alpha = 0.2,
              plot_data |> filter(pathogen == x_pathogen)) +
    
    geom_line(aes(x = notification_date, y = value_trend),
              plot_data |> filter(pathogen == x_pathogen)) +
    
    
    scale_x_date(breaks = scales::breaks_width("3 months", offset = -31),
                 labels = scales::label_date_short(format = c("%Y", "%b"))) +
    
    xlab("Notification date") + ylab("Case incidence") +
    
    coord_cartesian(xlim = c(ymd("2023-01-01"), ymd("2026-06-01"))) +
    
    ggh4x::facetted_pos_scales(y = scales) +
    
    plot_theme_paper
}


plot_pathogen("RSV") +
  ggtitle("RSV case incidence — Victoria, Australia") +
  
  coord_cartesian(ylim = c(0, 500), xlim = c(ymd("2023-01-01"), ymd("2026-01-01"))) +
  
  geom_vline(aes(xintercept = year_start),
             linetype = "44",
             tibble(year_start = seq(ymd("2023-01-01"), ymd("2026-01-01"), "year")))
