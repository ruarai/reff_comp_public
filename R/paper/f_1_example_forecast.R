


source("R/functions.R")
source("R/paper/plot_theme.R")

library(pammtools)


round_id <- "2025-07-03"
x_pathogen <- "flu"
x_location <- "VIC"
training_counts <- arrow::read_parquet(str_c("saved_results/results_baseline/", round_id, "/", x_pathogen, "_", x_location, "_training_counts.parquet"))

origin_date <- training_counts$origin_date[[1]] - days(1)


eval_counts <- read_csv("data/combined-count-smoothed-2026-04-22.csv", show_col_types = FALSE) |> 
  filter(notification_date > origin_date, location == x_location, pathogen == x_pathogen)

school_holidays <- get_school_holidays() |> 
  filter(location == x_location)

public_holidays <- get_public_holidays() |> 
  filter(location == x_location)

results <- read_rds(str_c("saved_results/results_baseline/", round_id, "/", x_pathogen, "_", x_location, ".rds"))

date_range <- c(ymd("2024-05-01"), ymd("2025-10-01"))

year_starts <- tibble(year_start = ymd(c("2024-01-01", "2025-01-01", "2026-01-01")))

plot_common <- list(
  geom_rect(aes(xmin = start_date, xmax = end_date + 1, ymin = -Inf, ymax = Inf),
            school_holidays, alpha = 0.1),
  scale_x_date(breaks = scales::date_breaks(width = "3 months"),
               labels = scales::label_date_short()),
  plot_theme_paper
)
plot_common_above <- list(
  geom_vline(xintercept = origin_date, linetype = "44", linewidth = 0.5),
  geom_vline(aes(xintercept = year_start), year_starts, linewidth = 0.5, alpha = 0.5, linetype = "24")
)

p_case_inc <- ggplot() +
  
  plot_common +
  
  geom_line(aes(x = date + 0.5, y = y),
            training_counts) +
  
  plot_common_above +
  
  coord_cartesian(xlim = date_range, ylim = c(0, 1200)) +
  
  xlab(NULL) + ylab("Incidence") +
  
  theme(legend.position = "none",
        panel.grid.major = element_blank()) +
  
  ggtitle(NULL, "Influenza case incidence — Victoria, Australia")

p_case_inc

trajectory_data <- results |> 
  spread_draws(c(log_R, inf, inf_delayed, cases, f_annual,
                 f_noise, latent_cases_log)[ix_t],
               log_R_0) |> 
  mutate(date = training_counts$date[ix_t],
         latent_cases = exp(latent_cases_log),
         holiday_effect = log_R - log_R_0 - f_annual - f_noise)

p_hol <- ggplot(trajectory_data |> filter(date >= date_range[1], date <= date_range[2])) +
  
  plot_common +
  
  fastribbon::stat_fastribbon(aes(x = date, y = log_R - f_noise - f_annual - log_R_0), geom = "stepribbon") +
  
  plot_common_above +
  
  geom_hline(yintercept = 0, linetype = "22", linewidth = 0.7) +
  
  scale_fill_brewer() +
  
  coord_cartesian(xlim = date_range, ylim = c(-0.25, 0.25)) +
  scale_y_continuous(breaks = c(-0.2, 0, 0.2)) +
  
  xlab(NULL) + ylab("Effect") +
  
  theme(legend.position = "none") +
  
  ggtitle(NULL, "Holiday effects <i>ω</i><sub>S</sub>(<i>t</i>) + <i>ω</i><sub>P</sub>(<i>t</i>)")  +
  
  theme(legend.position = "none",
        plot.subtitle = ggtext::element_markdown(),
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.margin = margin(t = 0.25, b = 0.25, l = 0.1, r = 0.1, unit = "cm"))

p_hol

p_noise <- ggplot(trajectory_data |> filter(date >= date_range[1], date <= date_range[2])) +
  
  plot_common +
  
  fastribbon::stat_fastribbon(aes(x = date, y = f_noise), geom = "stepribbon") +
  
  plot_common_above +
  
  geom_hline(yintercept = 0, linetype = "22", linewidth = 0.7) +
  
  scale_fill_brewer(palette = "Blues") +
  
  coord_cartesian(xlim = date_range, ylim = c(-0.15, 0.15)) +
  scale_y_continuous(breaks = c(-0.1, 0, 0.1)) +
  
  xlab(NULL) + ylab("Effect") +
  
  theme(legend.position = "none",
        plot.subtitle = ggtext::element_markdown(),
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.margin = margin(t = 0.25, b = 0.25, l = 0.1, r = 0.1, unit = "cm")) +
  
  ggtitle(NULL, "Noise <i>f</i><sub>noise</sub>(<i>t</i>)")

p_ann <- ggplot(trajectory_data |> filter(date >= date_range[1], date <= date_range[2])) +
  
  plot_common +
  
  fastribbon::stat_fastribbon(aes(x = date, y = f_annual), geom = "stepribbon") +
  
  plot_common_above +
  
  geom_hline(yintercept = 0, linetype = "22", linewidth = 0.7) +
  
  scale_fill_brewer(palette = "Blues") +
  
  coord_cartesian(xlim = date_range, ylim = c(-0.15, 0.15)) +
  scale_y_continuous(breaks = c(-0.1, 0, 0.1)) +
  
  xlab(NULL) + ylab("Effect") +
  
  theme(legend.position = "none",
        plot.subtitle = ggtext::element_markdown(),
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.margin = margin(t = 0.25, b = 0.25, l = 0.1, r = 0.1, unit = "cm")) +
  
  ggtitle(NULL, "Annual periodic variation <i>f</i><sub>ann</sub>(<i>t</i>)")


p_Rt_exp <- ggplot(trajectory_data |> filter(date >= date_range[1], date <= date_range[2])) +
  
  plot_common +
  
  fastribbon::stat_fastribbon(aes(x = date, y = exp(log_R)), geom = "stepribbon") +
  
  plot_common_above +
  
  geom_hline(yintercept = 1, linetype = "22", linewidth = 0.7) +
  
  scale_fill_brewer(palette = "Blues") +
  
  coord_cartesian(xlim = date_range, ylim = c(0.7, 1.3)) +
  scale_y_continuous(breaks = c(0.75, 1, 1.25)) +
  
  xlab(NULL) + ylab("R(t)") +
  
  theme(legend.position = "none",
        plot.subtitle = ggtext::element_markdown(),
        panel.grid.major.x = element_blank(),
        plot.margin = margin(t = 0.25, b = 0.25, l = 0.1, r = 0.1, unit = "cm")) +
  
  ggtitle(NULL, "Instantaneous reproduction number <i>R</i>(<i>t</i>)")

p_cases_long <- ggplot(trajectory_data |> filter(date >= date_range[1], date <= date_range[2])) +
  
  plot_common +
  
  fastribbon::stat_fastribbon(aes(x = date - 0.5, y = cases), geom = "stepribbon",
                              linewidth = 0.0) +
  
  plot_common_above +
  
  geom_line(aes(x = date, y = y),
            linewidth = 0.5,
            training_counts |> filter(date >= date_range[1])) +
  
  geom_line(aes(x = notification_date, y = value),
            colour = "orangered",
            eval_counts |> filter(notification_date <= date_range[2])) +
  
  scale_fill_brewer(palette = "Blues") +
  
  coord_cartesian(xlim = date_range, ylim = c(0, 1500)) +
  
  scale_y_continuous(labels = scales::label_comma()) +
  
  xlab(NULL) + ylab("Incidence") +
  
  theme(legend.position = "none",
        plot.margin = margin(t = 0.25, b = 0.25, l = 0.1, r = 0.1, unit = "cm"))

p_cases_long


p_Rt <- ((p_cases_long + theme(axis.text.x = element_blank())) / p_Rt_exp) + 
  plot_annotation("<b>A</b>. Conditioned case incidence and inferred <i>R</i>(<i>t</i>)",
                  theme = theme(plot.title = element_markdown(family = "Helvetica"),
                                plot.margin = margin(t = 0.2, l = 0.2, b = 0.2, r = 0.2, unit = "cm")))


p_decomp <- (p_hol / p_ann / (p_noise + theme(axis.text.x = element_text(margin = margin(t = 0.3, unit = "cm"), colour = "black")))) + 
  plot_annotation("<b>B</b>. Decomposed log <i>R</i>(<i>t</i>) additive effects",
                  theme = theme(plot.title = element_markdown(family = "Helvetica"),
                                plot.margin = margin(t = 0.2, l = 0.2, b = 0.2, r = 0.2, unit = "cm")))


n_days_ahead <- 7 * 4

date_range_forecast <- c(origin_date - days(7 * 8), origin_date + days(n_days_ahead))

p_cases <- ggplot(trajectory_data |> filter(date >= date_range_forecast[1], date <= date_range_forecast + days(1))) +
  
  plot_common +
  
  fastribbon::stat_fastribbon(aes(x = date - 0.5, y = cases), geom = "stepribbon",
                              linewidth = 0.0) +
  
  plot_common_above +
  
  geom_line(aes(x = date, y = y),
            training_counts |> filter(date >= date_range_forecast[1])) +
  geom_point(aes(x = date, y = y), size = 0.9,
             training_counts |> filter(date >= date_range_forecast[1])) +
  
  geom_line(aes(x = notification_date, y = value),
            colour = "orangered",
            eval_counts |> filter(notification_date <= date_range_forecast[2])) +
  geom_point(aes(x = notification_date, y = value), size = 0.9,
             colour = "orangered",
             eval_counts |> filter(notification_date <= date_range_forecast[2])) +
  
  scale_fill_brewer(palette = "Blues") +
  
  coord_cartesian(xlim = date_range_forecast, ylim = c(0, 1500)) +
  
  scale_x_date(date_breaks = "2 weeks", labels = scales::label_date_short()) +
  
  scale_y_continuous(labels = scales::label_comma()) +
  
  xlab(NULL) + ylab("Incidence") +
  
  theme(legend.position = "none",
        plot.margin = margin(t = 0.25, b = 0.25, l = 0.1, r = 0.1, unit = "cm"))

p_cases


p_forecast <- (p_cases) + 
  plot_annotation("<b>C</b>. Forecast case incidence",
                  theme = theme(plot.title = element_markdown(family = "Helvetica"))) +
  
  theme(plot.margin = margin(t = 0.2, l = 0.2, b = 0.2, r = 0.2, unit = "cm"))

p_example <- (wrap_elements(full = p_Rt) | wrap_elements(full = p_decomp)) /
  wrap_elements(full = p_forecast)

ggsave(
  "paper_figures/f_1_example.png", p_example,
  width = 9, height = 7,
  bg = "white"
)


  