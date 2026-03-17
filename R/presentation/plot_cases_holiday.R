

library(tidyverse)
source("R/functions.R")
source("R/presentation/plot_theme.R")


data <- read_csv("../forecast-hub/processed-data/2025/2025-09-25/combined-count-2025-09-25.csv")

x_location <- "NSW"
x_pathogen <- "RSV"

date_range <- c(ymd("2025-03-01"), ymd("2025-07-01"))

school_hols <- get_school_holidays() |> 
  filter(location == x_location,
         end_date >= date_range[1],
         start_date <= date_range[2])


data |> 
  filter(location == x_location, pathogen == x_pathogen,
         notification_date >= date_range[1],
         notification_date < date_range[2]) |> 
  
  ggplot() +
  
  geom_rect(aes(xmin = start_date, xmax = end_date, ymin = -Inf, ymax = Inf),
            
            alpha = 0.3, fill = ggokabeito::palette_okabe_ito(5),
            
            school_hols) +
  
  geom_line(aes(x = notification_date, y = cases)) +
  
  xlab("Date") + ylab("Case incidence") +
  
  scale_x_date(labels = scales::label_date_short()) +
  
  coord_cartesian(ylim = c(0, NA)) +
  
  ggtitle("RSV case incidence — New South Wales") +
  
  plot_theme_pres

