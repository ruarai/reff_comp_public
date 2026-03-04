

source("R/functions.R")
source("R/presentation/plot_theme.R")

school_hols <- get_school_holidays() |> 
  filter(location == "VIC")

school_hols_long <- get_school_holidays_long() |> 
  filter(location == "VIC")


public_holidays <- get_public_holidays() |> 
  filter(location == "VIC")

date_range <- c(ymd("2024-12-01"), ymd("2025-06-15"))

school_holiday_effect <- tibble(
  date = seq(date_range[1], date_range[2], "days"),
  location = "VIC",
  
  summer_effect = -0.1,
  winter_effect = -0.2
) |> 
  left_join(school_hols_long) |> 
  mutate(school_holiday_type = replace_na(school_holiday_type, 0) + 1) |> 
  
  mutate(effect = case_when(
    school_holiday_type == 1 ~ 0,
    school_holiday_type == 2 ~ summer_effect,
    school_holiday_type == 3 ~ winter_effect
  ))
  
  
annotations <- list(geom_segment(
  data = data.frame(
    x = ymd(c("2025-02-02", "2025-04-26")),
    y = c(-0.0102, -0.0109),
    xend = ymd(c("2025-02-02", "2025-04-26")),
    yend = c(-0.0924, -0.1952)
  ),
  mapping = aes(x = x, y = y, xend = xend, yend = yend),
  arrow = arrow(30L, unit(0.1, "inches"), "last", "closed"),
  linewidth = 1.0,
  inherit.aes = FALSE
),
  ggtext::geom_richtext(
    data = data.frame(
      x = as.Date(c("2025-05-25", "2025-02-21")),
      y = c(-0.1, -0.05),
      label = c("<i>ω</i><sub>between-term</sub>", "<i>ω</i><sub>summer</sub>")
    ),
    label.colour = "white",
    size = 8,
    mapping = aes(x = x, y = y, label = label),
    size = 3.88,
    fontface = "plain",
    inherit.aes = FALSE
  ))

ggplot(school_holiday_effect) +
  
  geom_rect(aes(xmin = start_date, xmax = end_date + 1, ymin = -Inf, ymax = Inf, fill = factor(school_holiday_type)),
            
            alpha = 0.2,
            
            school_hols) +
  
  ggokabeito::scale_fill_okabe_ito(order = c(3, 5)) +
  
  geom_step(aes(x = date, y = effect),
            linewidth = 0.75) +
  
  coord_cartesian(xlim = date_range, ylim = c(-0.25, 0.1)) +
  
  scale_x_date(breaks = scales::date_breaks(),
               labels = scales::label_date_short(format = c("", "%b", "%d", "%H:%M"))) +
  
  xlab("Date") + ylab("School holiday effect") +
  
  plot_theme_pres +
  
  annotations +
  
  theme(legend.position = "none")



public_holiday_effects <- tibble(
  holiday_name = c("New Year's Day", "Australia Day", "Labour Day", "Easter", 
              "ANZAC Day", "King's Birthday", "Friday before the AFL Grand Final", 
              "Melbourne Cup", "Christmas Day", "Boxing Day"),
  effect = c(0.06, -0.015, -0.027, -0.026, -0.07, -0.029, 0.03, -0.096, 
             -0.02, -0.021)
)

public_holiday_effect <- tibble(
  date = seq(date_range[1], date_range[2], by = "days")
) |> 
  left_join(public_holidays) |> 
  left_join(public_holiday_effects) |> 
  mutate(effect = replace_na(effect, 0))

p_public <- ggplot(public_holiday_effect) +
  geom_step(aes(x = date, y = effect),
            linewidth = 0.75) +
  
  coord_cartesian(xlim = date_range, ylim = c(-0.08, 0.08)) +
  
  scale_x_date(breaks = scales::date_breaks(),
               labels = scales::label_date_short(format = c("", "%b", "%d", "%H:%M"))) +
  
  xlab(NULL) + ylab("<i>ω</i><sub>P</sub>(t)") +
  ggtitle("Public holiday effect <i>ω</i><sub>P</sub>(t)") +
  
  plot_theme_pres +
  
  theme(legend.position = "none",
        axis.title.y = element_markdown(),
        plot.margin = margin(t = 0.1, b = 0.1, l = 0.1, r = 0.1),
        plot.title = element_markdown())

p_school <- ggplot(school_holiday_effect) +
  
  geom_step(aes(x = date, y = effect),
            linewidth = 0.75) +
  
  coord_cartesian(xlim = date_range, ylim = c(-0.25, 0.1)) +
  
  scale_x_date(breaks = scales::date_breaks(),
               labels = scales::label_date_short(format = c("", "%b", "%d", "%H:%M"))) +
  
  xlab(NULL) + ylab("<i>ω</i><sub>S</sub>(t)") +
  ggtitle("School holiday effect <i>ω</i><sub>S</sub>(t)") +
  
  plot_theme_pres +
  
  theme(legend.position = "none",
        axis.title.y = element_markdown(),
        plot.margin = margin(t = 0.1, b = 0.1, l = 0.1, r = 0.1),
        plot.title = element_markdown())

n_days <- as.numeric(date_range[2] - date_range[1])

effect_noise <- tibble(t = 0:n_days) |>
  mutate(y_noise = simulate_matern(t, 0.05, 14, seed = 3),
         y_noise = y_noise + simulate_matern(t, 0.2, 120, seed = 4),
         date = date_range[1] + t)
  
p_noise <- ggplot(effect_noise) +
  
  geom_line(aes(x = date, y = y_noise),
            linewidth = 0.75) +
  
  scale_x_date(breaks = scales::date_breaks(),
               labels = scales::label_date_short(format = c("", "%b", "%d", "%H:%M"))) +
  
  scale_y_continuous(breaks = c(-0.2, 0, 0.2)) +
  
  xlab(NULL) + ylab("<i>f</i><sub>noise</sub>(t) + <i>f</i><sub>ann</sub>(t)") +
  ggtitle("Noise terms <i>f</i><sub>noise</sub>(t) + <i>f</i><sub>ann</sub>(t)") +
  
  plot_theme_pres +
  
  theme(legend.position = "none",
        axis.title.y = element_markdown(),
        plot.margin = margin(t = 0.1, b = 0.1, l = 0.1, r = 0.1),
        plot.title = element_markdown())


p_public / p_school / p_noise
