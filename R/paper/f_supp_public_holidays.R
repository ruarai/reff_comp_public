
source("R/functions.R")
source("R/paper/common.R")


public_holidays <- get_public_holidays()



combined_counts <- read_csv("data/combined-count-smoothed-2026-04-15.csv")


public_holidays |> 
  expand_grid(pathogen = c("RSV", "SARSCOV2", "flu")) |> 
  filter(location == "NSW") |> 
  
  group_by(pathogen) |> 
  
  filter(!(holiday_name == "Easter" & lag(holiday_name) == "Easter"),
         holiday_name != "Boxing Day",
         !str_detect(holiday_name, "Observed")) |> 
  
  expand_grid(days_around = -14:14) |> 
  mutate(date_offset = date + days_around) |> 
  
  left_join(combined_counts |> select(date_offset = notification_date, location, pathogen, value, value_trend, value_trend_and_remainder)) |> 
  
  group_by(location, pathogen, date, holiday_name) |> 
  mutate(y = value_trend_and_remainder / mean(value_trend_and_remainder)) |> 
  
  label_pathogen() |> 
  
  ggplot() +
  
  geom_line(aes(x = days_around, y = y, group = interaction(date, pathogen),
            colour = factor(year(date)))) +
  scale_colour_brewer(type = "qual", palette = 2, name = "Year") +
  
  ylab("Relative case counts") +
  
  geom_vline(xintercept = 0, linetype = "44") +
  
  xlab("Days around public holiday") +
  
  coord_cartesian(ylim = c(0, 2)) +
  
  facet_grid(rows = vars(holiday_name), cols = vars(pathogen_label)) +
  
  plot_theme_paper +
  
  theme(legend.position = "bottom") +
  
  ggtitle(NULL, "Relative case incidence around public holidays — New South Wales, Australia")


public_holidays |> 
  expand_grid(pathogen = c("RSV", "SARSCOV2", "flu")) |> 
  
  group_by(pathogen, location) |> 
  
  filter(!(holiday_name == "Easter" & lag(holiday_name) == "Easter"),
         !str_detect(holiday_name, "Observed"),
         holiday_name %in% c("Labour Day", "Easter", "Christmas Day", "ANZAC Day",
                             "Australia Day", "New Year's Day"),
         location != "NZ") |> 
  
  expand_grid(days_around = -14:14) |> 
  mutate(date_offset = date + days_around) |> 
  
  left_join(combined_counts |> select(date_offset = notification_date, location, pathogen, value, value_trend, value_trend_and_remainder)) |> 
  
  group_by(location, pathogen, date, holiday_name) |> 
  mutate(y = value_trend_and_remainder / mean(value_trend_and_remainder)) |> 
  group_by(holiday_name, days_around, pathogen) |> 
  summarise(y = median(y, na.rm = TRUE)) |> 
  
  ggplot() +
  
  geom_line(aes(x = days_around, y = y)) +
  
  geom_vline(xintercept = 0, linetype = "44") +
  
  xlab("Days around public holiday") +
  
  coord_cartesian(ylim = c(0, 2)) +
  
  facet_grid(rows = vars(holiday_name), cols = vars(pathogen)) +
  
  plot_theme_paper +
  
  theme(legend.position = "bottom")

plot_data_public <- public_holidays |> 
  filter(location %in% locations_paper) |> 
  filter(year(date) == 2024) |> 
  
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


plot_data_school <- get_school_holidays() |> 
  filter(location %in% locations_paper) |> 
  filter(year(end_date) == 2024 | year(start_date) == 2024)


ggplot() +
  
  geom_tile(aes(x = date, y = location, fill = holiday_name),
            plot_data_public) +
  
  geom_text(aes(x = date, y = location, label = holiday_label),
            size = 3,
            plot_data_public_labels) +
  
  geom_linerange(aes(xmin = start_date, xmax = end_date, y = location),
                 plot_data_school, linewidth = 20.0,
                 alpha = 0.1) +
  
  scale_fill_discrete(name = "Public holiday name") +
  
  xlab(NULL) + ylab("Region") +
  
  coord_cartesian(xlim = c(ymd("2024-01-01"), ymd("2025-01-01"))) +
  
  plot_theme_paper +
  
  theme(legend.position = "bottom")


