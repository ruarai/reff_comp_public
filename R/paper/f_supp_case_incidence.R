source("R/functions.R")
source("R/paper/common.R")

combined_counts <- read_csv("data/combined-count-smoothed_final_anonymised.csv")

training_dates <- get_training_dates()

plot_data <- combined_counts |> 
  label_pathogen() |> 
  left_join(training_dates) |> 
  filter(notification_date >= cutoff_date, notification_date < ymd("2026-01-01"))


school_holidays <- get_school_holidays() |> 
  filter(location %in% locations_paper,
         start_date <= ymd("2026-01-01"))

public_holidays <- get_public_holidays() |> 
  filter(location %in% locations_paper,
         date <= ymd("2026-01-01"))

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
    
    geom_rect(aes(xmin = date - 0.5, xmax = date + 0.5, ymin = -Inf, ymax = Inf),
              alpha = 0.2, fill = "#e31a1c",
              public_holidays) +
    geom_rect(aes(xmin = start_date - 0.5, xmax = end_date + 0.5, ymin = -Inf, ymax = Inf),
              school_holidays, alpha = 0.1) +
    
    geom_line(aes(x = notification_date, y = value), alpha = 0.2,
              plot_data |> filter(pathogen == x_pathogen)) +
    
    geom_line(aes(x = notification_date, y = value_trend_and_remainder),
              plot_data |> filter(pathogen == x_pathogen)) +
    
    facet_wrap(~location * pathogen, ncol = 1, 
               labeller = \(x) label_value(x, multi_line = FALSE),
               scales = "free_y") +
    
    scale_x_date(breaks = scales::breaks_width("3 months", offset = -31),
                 labels = scales::label_date_short(format = c("%Y", "%b"))) +
    
    xlab("Notification date") + ylab("Case incidence") +
    
    coord_cartesian(xlim = c(ymd("2023-01-01"), ymd("2026-01-01"))) +
    
    ggh4x::facetted_pos_scales(y = scales) +
    
    plot_theme_paper
}


p_RSV <- plot_pathogen("RSV")
p_flu <- plot_pathogen("flu")
p_SARSCOV2 <- plot_pathogen("SARSCOV2")


ggsave(
  "paper_figures/supp/f_incidence_RSV.png", p_RSV,
  width = 10, height = 8,
  bg = "white"
)
ggsave(
  "paper_figures/supp/f_incidence_flu.png", p_flu,
  width = 10, height = 8,
  bg = "white"
)
ggsave(
  "paper_figures/supp/f_incidence_SARSCOV2.png", p_SARSCOV2,
  width = 10, height = 8,
  bg = "white"
)


