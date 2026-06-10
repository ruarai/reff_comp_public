public_holidays <- get_public_holidays() |> 
  filter(date >= ymd("2027-01-01"))

plot_data_public <- public_holidays |> 
  
  arrange(date) |> 
  mutate(holiday_name = factor(holiday_name, unique(holiday_name)))

plot_data_public_labels <- plot_data_public |> 
  group_by(location, holiday_name) |> 
  slice(1) |> 
  
  group_by(holiday_name, date) |> 
  slice(1) |> 
  mutate(holiday_label = str_wrap(holiday_name, width = 10),
         location = if_else(holiday_name == "Boxing Day", "NSW", location),
         location = if_else(holiday_name == "King's Birthday" & date < ymd("2024-10-01"),
                            "NSW", location))


plot_data_school <- get_school_holidays()


ggplot() +
  
  geom_tile(aes(x = date, y = location, fill = holiday_name),
            plot_data_public) +
  
  geom_text(aes(x = date, y = location, label = holiday_label),
            size = 3,
            plot_data_public_labels) +
  
  geom_linerange(aes(xmin = start_date, xmax = end_date, y = location),
                 plot_data_school, linewidth = 10.0,
                 alpha = 0.1) +
  
  scale_fill_discrete(name = "Public holiday name") +
  
  xlab(NULL) + ylab("Region") +
  
  coord_cartesian(xlim = c(ymd("2027-01-01"), ymd("2028-01-01"))) +
  
  plot_theme_paper +
  
  theme(legend.position = "bottom")
