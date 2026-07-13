


tbl_all <- read_csv("data/location_age_hours_all.csv", n_max = 5600)

count_worked_less <- read_csv(
  "data/count_worked_less.csv", skip = 11, n_max = 5371 - 11,
  col_names = c(
    "date", "location", "age_group", "method",
    "count", "count_rse", "count_annotation"
  )
) |>
  fill(date, location, age_group) |> 
  select(-X8)

count_all <- read_csv(
  "data/count_all.csv", skip = 8, n_max = 5371 - 11,
  col_names = c(
    "date", "location", "age_group", "method",
    "count_all", "count_all_rse", "count_all_annotation"
  )
) |>
  fill(date, location, age_group) |> 
  select(-X8)

plot_data <- count_worked_less |>
  left_join(count_all) |> 
  
  
  filter(date != "Total") |> 
  filter(location %in% c("Queensland", "New South Wales", "Victoria", "Western Australia")) |> 
  
  mutate(date = my(date, quiet = FALSE))


plot_data |> 

  ggplot() +
  
  geom_line(aes(x = date, y = count / count_all, colour = age_group)) +
  
  scale_x_date(date_breaks = "1 year",
               labels = scales::label_date_short()) +
  
  facet_wrap(~location, scales = "free_y") +
  
  coord_cartesian(ylim = c(0, NA)) +
  
  theme_bw()



plot_data_summed <- plot_data |> 
  
  group_by(date, location) |> 
  summarise(count = sum(count),
            count_all = sum(count_all))

fit_data <- plot_data_summed |> 
  ungroup() |> 
  mutate(t = as.numeric(date - min(date)),
         doy = yday(date),
         
         count = round(count), count_all = round(count_all),
         
         count_other = count_all - count)


library(mgcv)

fit <- gam(
  cbind(count, count_other) ~ s(t, k = 20) + s(doy, bs = "cc") + location,
  data = fit_data, family = "binomial"
)

gratia::draw(fit, rug = FALSE)

plot_data_fitted <- fit_data %>%
  mutate(pred = predict(fit, newdata = .))


ggplot() +
  
  geom_line(aes(x = date, y = count / count_all),
            plot_data_summed) +
  
  geom_line(aes(x = date, y = plogis(pred)),
            colour = "red", alpha = 0.5,
            plot_data_fitted) +
  
  scale_x_date(date_breaks = "1 year",
               labels = scales::label_date_short()) +
  
  # facet_grid(cols = vars(location)) +
  facet_wrap(~location) +
  
  coord_cartesian(ylim = c(0, NA)) +
  
  theme_bw()


 |> 
  
  ggplot() +
  
  geom_line(aes(x = date, y = count / count_all)) +
  
  scale_x_date(date_breaks = "1 year",
               labels = scales::label_date_short()) +
  
  facet_wrap(~location, scales = "free_y") +
  
  coord_cartesian(ylim = c(0, NA)) +
  
  theme_bw()
