


source("R/functions.R")
source("R/presentation/plot_theme.R")

library(pammtools)


round_id <- "2025-07-03"
# round_id <- "2025-10-09"
x_pathogen <- "flu"
x_location <- "VIC"
training_counts <- arrow::read_parquet(str_c("saved_results/results_baseline/", round_id, "/", x_pathogen, "_", x_location, "_training_counts.parquet"))
eval_counts <- arrow::read_parquet(str_c("saved_results/results_baseline/2025-10-09/", x_pathogen, "_", x_location, "_training_counts.parquet"))

origin_date <- training_counts$origin_date[[1]]

school_holidays <- get_school_holidays() |> 
  filter(location == x_location)

public_holidays <- get_public_holidays() |> 
  filter(location == x_location)

results <- read_rds(str_c("saved_results/results_baseline/", round_id, "/", x_pathogen, "_", x_location, ".rds"))

date_range <- c(ymd("2025-05-15"), ymd("2025-09-10"))
date_range <- c(ymd("2023-12-01", ymd("2025-12-01")))

p_case_inc <- ggplot() +
  
  geom_rect(aes(xmin = start_date, xmax = end_date + 1, ymin = -Inf, ymax = Inf),
            school_holidays, alpha = 0.1) +
  
  geom_vline(xintercept = training_counts$origin_date[[1]], linetype = "44", linewidth = 1.0) +
  
  geom_line(aes(x = date + 0.5, y = y),
            training_counts) +
  
  scale_fill_brewer(palette = "Greens") +
  
  coord_cartesian(xlim = date_range, ylim = c(0, 1200)) +
  
  scale_x_date(breaks = scales::date_breaks(width = "3 months"),
               labels = scales::label_date_short()) +
  
  xlab(NULL) + ylab("Incidence") +
  
  plot_theme_pres +
  
  theme(legend.position = "none",
        panel.grid.major = element_blank()) +
  
  ggtitle(NULL, "Case incidence — Victoria, influenza")

p_case_inc

trajectory_data <- results |> 
  spread_draws(log_R[ix_t], inf[ix_t], inf_delayed[ix_t],
               cases[ix_t], f_annual[ix_t], f_noise[ix_t],
               latent_cases_log[ix_t], log_R_0) |> 
  mutate(date = training_counts$date[ix_t],
         latent_cases = exp(latent_cases_log)) |> 
  
  mutate(holiday_effect = log_R - log_R_0 - f_annual - f_noise)

hide_forecast <- FALSE


p_hol <- ggplot(trajectory_data |> filter(!(hide_forecast & date > origin_date), date >= date_range[1], date <= date_range[2])) +
  
  geom_rect(aes(xmin = date, xmax = date + 1, ymin = -Inf, ymax = Inf),
            public_holidays, alpha = 0.1)  +
  
  geom_rect(aes(xmin = start_date, xmax = end_date + 1, ymin = -Inf, ymax = Inf),
            school_holidays, alpha = 0.1)  +
  
  geom_rect(aes(xmin = date, xmax = date + 1, ymin = -Inf, ymax = Inf),
            public_holidays, alpha = 0.1)  +
  
  fastribbon::stat_fastribbon(aes(x = date, y = log_R - f_noise - f_annual - log_R_0), geom = "stepribbon") +
  
  geom_vline(xintercept = training_counts$origin_date[[1]], linetype = "44", linewidth = 1.0) +
  geom_hline(yintercept = 0, linetype = "22", linewidth = 0.7) +
  
  scale_fill_brewer() +
  
  coord_cartesian(xlim = date_range, ylim = c(-0.25, 0.25)) +
  scale_y_continuous(breaks = c(-0.2, 0, 0.2)) +
  
  scale_x_date(breaks = scales::date_breaks(width = "3 months"),
               labels = scales::label_date_short()) +
  
  xlab(NULL) + ylab("Effect") +
  
  plot_theme_pres +
  
  theme(legend.position = "none") +
  
  ggtitle(NULL, "Holiday effects <i>ω</i><sub>S</sub>(<i>t</i>) + <i>ω</i><sub>P</sub>(<i>t</i>)")  +
  
  plot_theme_pres +
  
  theme(legend.position = "none",
        plot.subtitle = ggtext::element_markdown(),
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.margin = margin(t = 0.25, b = 0.25, l = 0.1, r = 0.1, unit = "cm"))


p_noise <- ggplot(trajectory_data |> filter(date >= date_range[1], date <= date_range[2])
                  |> filter(!(hide_forecast & date > origin_date))) +
  
  fastribbon::stat_fastribbon(aes(x = date, y = f_noise), geom = "stepribbon") +
  
  geom_vline(xintercept = training_counts$origin_date[[1]], linetype = "44", linewidth = 1.0) +
  geom_hline(yintercept = 0, linetype = "22", linewidth = 0.7) +
  geom_line(aes(x = date + 0.5, y = y),
            training_counts) +
  
  geom_point(aes(x = date + 0.5, y = y),
             training_counts) +
  
  scale_fill_brewer(palette = "Blues") +
  
  coord_cartesian(xlim = date_range, ylim = c(-0.15, 0.15)) +
  scale_y_continuous(breaks = c(-0.1, 0, 0.1)) +
  
  scale_x_date(breaks = scales::date_breaks(width = "3 months"),
               labels = scales::label_date_short()) +
  
  xlab(NULL) + ylab("Effect") +
  
  plot_theme_pres +
  
  theme(legend.position = "none",
        plot.subtitle = ggtext::element_markdown(),
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.margin = margin(t = 0.25, b = 0.25, l = 0.1, r = 0.1, unit = "cm")) +
  
  ggtitle(NULL, "Noise <i>f</i><sub>noise</sub>(<i>t</i>)")

p_ann <- ggplot(trajectory_data |> filter(date >= date_range[1], date <= date_range[2]) |> 
                  filter(!(hide_forecast & date > origin_date))) +
  
  fastribbon::stat_fastribbon(aes(x = date, y = f_annual), geom = "stepribbon") +
  
  geom_vline(xintercept = training_counts$origin_date[[1]], linetype = "44", linewidth = 1.0) +
  geom_hline(yintercept = 0, linetype = "22", linewidth = 0.7) +
  
  geom_line(aes(x = date + 0.5, y = y),
            training_counts) +
  
  geom_point(aes(x = date + 0.5, y = y),
             training_counts) +
  
  scale_fill_brewer(palette = "Blues") +
  
  coord_cartesian(xlim = date_range, ylim = c(-0.15, 0.15)) +
  scale_y_continuous(breaks = c(-0.1, 0, 0.1)) +
  
  scale_x_date(breaks = scales::date_breaks(width = "3 months"),
               labels = scales::label_date_short()) +
  
  xlab(NULL) + ylab("Effect") +
  
  plot_theme_pres +
  
  theme(legend.position = "none",
        plot.subtitle = ggtext::element_markdown(),
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.margin = margin(t = 0.25, b = 0.25, l = 0.1, r = 0.1, unit = "cm")) +
  
  ggtitle(NULL, "Annual periodic variation <i>f</i><sub>ann</sub>(<i>t</i>)")


p_Rt_exp_forecast <- ggplot(trajectory_data |> filter(date >= date_range[1], date <= date_range[2]) |> 
                              filter(!(hide_forecast & date > origin_date))) +
  
  fastribbon::stat_fastribbon(aes(x = date, y = exp(log_R)), geom = "stepribbon") +
  
  geom_rect(aes(xmin = start_date, xmax = end_date + 1, ymin = -Inf, ymax = Inf),
            school_holidays, alpha = 0.1) +
  
  geom_vline(xintercept = training_counts$origin_date[[1]], linetype = "44", linewidth = 1.0) +
  geom_hline(yintercept = 1, linetype = "22", linewidth = 0.7) +
  
  geom_line(aes(x = date + 0.5, y = y),
            training_counts) +
  
  geom_point(aes(x = date + 0.5, y = y),
             training_counts) +
  
  scale_fill_brewer(palette = "Blues") +
  
  coord_cartesian(xlim = date_range, ylim = c(0.6, 1.4)) +
  scale_y_continuous(breaks = c(0.75, 1, 1.25)) +
  
  scale_x_date(breaks = scales::date_breaks(width = "1 months"),
               labels = scales::label_date_short()) +
  
  xlab(NULL) + ylab("R(t)") +
  
  plot_theme_pres +
  
  theme(legend.position = "none",
        plot.subtitle = ggtext::element_markdown(),
        panel.grid.major.x = element_blank(),
        plot.margin = margin(t = 0.25, b = 0.25, l = 0.1, r = 0.1, unit = "cm")) +
  
  ggtitle(NULL, "Instantaneous reproduction number <i>R</i>(<i>t</i>)")


(p_case_inc + theme(axis.text.x = element_blank())) / p_Rt_exp_forecast

(p_Rt_exp_forecast + theme(axis.text.x = element_blank())) / p_hol / p_ann / 
  (p_noise + theme(axis.text.x = element_text()))

p_hol / p_ann / p_noise / p_Rt_exp_forecast

# (p_Rt_exp + ggtitle(NULL, "R(<i>t</i>) — forecasted")) / 
#   (p_Rt_exp_latest + ggtitle(NULL, "R(<i>t</i>) — end-of-season estimate"))

p_inf <- ggplot(trajectory_data) +
  
  geom_rect(aes(xmin = start_date, xmax = end_date + 1, ymin = -Inf, ymax = Inf),
            school_holidays, alpha = 0.1)  +
  
  geom_rect(aes(xmin = date, xmax = date + 1, ymin = -Inf, ymax = Inf),
            public_holidays, alpha = 0.1)  +
  
  fastribbon::stat_fastribbon(aes(x = date, y = inf), geom = "stepribbon") +
  
  geom_vline(xintercept = training_counts$origin_date[[1]], linetype = "44", linewidth = 1.0) +
  
  scale_fill_brewer(palette = "Oranges") +
  
  coord_cartesian(xlim = date_range, ylim = c(0, 1500)) +
  
  scale_x_date(breaks = scales::date_breaks(),
               labels = scales::label_date_short()) +
  
  xlab(NULL) + ylab("Incidence") +
  
  plot_theme_pres +
  
  theme(legend.position = "none") +
  
  ggtitle(NULL, "Latent infection incidence <i>z</i>(<i>t</i>)")  +
  
  plot_theme_pres +
  
  theme(legend.position = "none",
        plot.subtitle = ggtext::element_markdown(),
        axis.text.x = element_blank(),
        plot.margin = margin(t = 0.25, b = 0.25, l = 0.1, r = 0.1, unit = "cm"))

p_latent_cases <- ggplot(trajectory_data) +
  
  geom_rect(aes(xmin = start_date, xmax = end_date + 1, ymin = -Inf, ymax = Inf),
            school_holidays, alpha = 0.1)  +
  
  geom_rect(aes(xmin = date, xmax = date + 1, ymin = -Inf, ymax = Inf),
            public_holidays, alpha = 0.1)  +
  
  fastribbon::stat_fastribbon(aes(x = date, y = latent_cases), geom = "stepribbon") +
  
  geom_vline(xintercept = training_counts$origin_date[[1]], linetype = "44", linewidth = 1.0) +
  
  scale_fill_brewer(palette = "Greens") +
  
  coord_cartesian(xlim = date_range, ylim = c(0, 1500)) +
  
  scale_x_date(breaks = scales::date_breaks(),
               labels = scales::label_date_short()) +
  
  xlab(NULL) + ylab("Incidence") +
  
  ggtitle(NULL, "Latent case incidence <i>y</i>(<i>t</i>) (without day-of-week effect)")  +
  
  plot_theme_pres +
  
  theme(legend.position = "none",
        plot.subtitle = ggtext::element_markdown(),
        axis.text.x = element_blank(),
        plot.margin = margin(t = 0.25, b = 0.25, l = 0.1, r = 0.1, unit = "cm"))

p_cases <- ggplot(trajectory_data |> filter(date >= date_range[1], date <= date_range[2])) +
  
  fastribbon::stat_fastribbon(aes(x = date, y = cases), geom = "stepribbon") +
  
  # geom_rect(aes(xmin = date, xmax = date + 1, ymin = -Inf, ymax = Inf),
  #           public_holidays, alpha = 0.1) +
  # 
  # geom_rect(aes(xmin = start_date, xmax = end_date + 1, ymin = -Inf, ymax = Inf),
  #           school_holidays, alpha = 0.1) +
  
  geom_vline(xintercept = training_counts$origin_date[[1]], linetype = "44", linewidth = 1.0) +
  
  geom_line(aes(x = date + 0.5, y = y),
            training_counts) +
  
  geom_point(aes(x = date + 0.5, y = y),
             training_counts) +
  
  scale_fill_brewer(palette = "Greens") +
  
  coord_cartesian(xlim = date_range, ylim = c(0, 1500)) +
  
  scale_x_date(breaks = scales::date_breaks(),
               labels = scales::label_date_short()) +
  
  xlab(NULL) + ylab("Incidence") +
  
  plot_theme_pres +
  
  theme(legend.position = "none",
        plot.margin = margin(t = 0.25, b = 0.25, l = 0.1, r = 0.1, unit = "cm")) +
  
  ggtitle(NULL, "Forecast case incidence")

p_cases

(p_Rt_exp_forecast + theme(axis.text.x = element_blank())) / p_cases



p_cases +
  
  geom_line(aes(x = date + 0.5, y = y),
            eval_counts |> filter(date >= training_counts$origin_date[[1]],
                                  date <= date_range[2])) +
  
  geom_point(aes(x = date + 0.5, y = y), size = 1.5, stroke = 1.5,
             eval_counts |> filter(date >= training_counts$origin_date[[1]],
                                   date <= date_range[2])) +
  
  geom_point(aes(x = date + 0.5, y = y), size = 0.7, stroke = 0.7, colour = "white",
             eval_counts |> filter(date >= training_counts$origin_date[[1]],
                                   date <= date_range[2])) +
  
  coord_cartesian(xlim = c(date_range[1],
                           ymd("2025-08-01")),
                  ylim = c(0, 1700))
