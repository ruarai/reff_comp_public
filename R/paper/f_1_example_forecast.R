


source("R/functions.R")
source("R/paper/common.R")

library(pammtools)


round_id <- "2025-07-03"
x_pathogen <- "flu"
x_location <- "VIC"
job_name <- str_c(round_id, "/", x_pathogen, "_", x_location, "_case_incidence")
training_counts <- arrow::read_parquet(str_c("saved_results/baseline/", job_name, "_training_counts.parquet"))

origin_date <- training_counts$origin_date[[1]] - days(1)


eval_counts <- read_csv("data/combined-count-smoothed-2026-04-22.csv", show_col_types = FALSE) |> 
  filter(notification_date > origin_date, location == x_location, pathogen == x_pathogen)

school_holidays <- get_school_holidays() |> 
  filter(location == x_location)

public_holidays <- get_public_holidays() |> 
  filter(location == x_location)

results <- read_rds(str_c("saved_results/baseline/", job_name, ".rds"))

date_range <- c(ymd("2024-05-01"), ymd("2025-10-01"))

colour_fc <- "#02818a"

plot_margin <- 0.1

plot_common <- list(
  
  geom_rect(aes(xmin = date - 0.5, xmax = date + 0.5, ymin = -Inf, ymax = Inf),
            alpha = 0.2, fill = "#e31a1c",
            public_holidays),
  geom_rect(aes(xmin = start_date - 0.5, xmax = end_date + 0.5, ymin = -Inf, ymax = Inf),
            school_holidays, alpha = 0.1),
  scale_x_date(breaks = scales::date_breaks(width = "3 months"),
               labels = scales::label_date_short()),
  plot_theme_paper,
  theme(plot.margin = margin(t = plot_margin, b = plot_margin, l = 0.1, r = 0.1, unit = "cm"))
)
plot_common_above <- list(
  geom_vline(xintercept = origin_date, linetype = "44", linewidth = 0.5),
  geom_vline(xintercept = origin_date + days(28), linetype = "44", linewidth = 0.5, alpha = 0.5)
)

trajectory_data <- results |> 
  spread_draws(c(log_R, inf, inf_delayed, cases, f_annual,
                 f_noise, latent_cases_log)[ix_t],
               log_R_0) |> 
  mutate(date = training_counts$date[ix_t],
         latent_cases = exp(latent_cases_log),
         holiday_effect = log_R - log_R_0 - f_annual - f_noise)


traj_data_backcast_long <- trajectory_data |> 
  filter(date >= date_range[1], date <= origin_date + 1)

traj_data_forecast_long <- trajectory_data |> 
  filter(date >= origin_date + 1, date <= date_range + 1)

p_hol <- ggplot() +
  
  plot_common +
  
  fastribbon::stat_fastribbon(aes(x = date - 0.5, y = holiday_effect, alpha = after_stat(level)),
                              fill = "#3182BD",
                              geom = "stepribbon",
                              traj_data_backcast_long,
                              linewidth = 0.0) +
  
  fastribbon::stat_fastribbon(aes(x = date - 0.5, y = holiday_effect, alpha = after_stat(level)),
                              geom = "stepribbon",
                              fill = colour_fc,
                              traj_data_forecast_long,
                              linewidth = 0.0) +
  
  scale_alpha_manual(values = c(0.2, 0.5, 1)) +
  
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
        panel.grid.major.x = element_blank())

p_hol

p_noise <- ggplot() +
  
  plot_common +
  
  fastribbon::stat_fastribbon(aes(x = date - 0.5, y = f_noise, alpha = after_stat(level)),
                              fill = "#3182BD",
                              traj_data_backcast_long,
                              linewidth = 0.0) +
  
  fastribbon::stat_fastribbon(aes(x = date - 0.5, y = f_noise, alpha = after_stat(level)),
                              fill = colour_fc,
                              traj_data_forecast_long,
                              linewidth = 0.0) +
  
  plot_common_above +
  
  geom_hline(yintercept = 0, linetype = "22", linewidth = 0.7) +
  
  scale_fill_brewer(palette = "Blues") +
  
  coord_cartesian(xlim = date_range, ylim = c(-0.15, 0.15)) +
  scale_y_continuous(breaks = c(-0.1, 0, 0.1)) +
  
  xlab(NULL) + ylab("Effect") +
  
  theme(legend.position = "none",
        plot.subtitle = ggtext::element_markdown(),
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank()) +
  
  ggtitle(NULL, "Noise <i>f</i><sub>noise</sub>(<i>t</i>)")

p_ann <- ggplot() +
  
  plot_common +
  
  fastribbon::stat_fastribbon(aes(x = date - 0.5, y = f_annual, alpha = after_stat(level)),
                              fill = "#3182BD",
                              traj_data_backcast_long,
                              linewidth = 0.0) +
  
  fastribbon::stat_fastribbon(aes(x = date - 0.5, y = f_annual, alpha = after_stat(level)),
                              fill = colour_fc,
                              traj_data_forecast_long,
                              linewidth = 0.0) +
  
  plot_common_above +
  
  geom_hline(yintercept = 0, linetype = "22", linewidth = 0.7) +
  
  scale_fill_brewer(palette = "Blues") +
  
  coord_cartesian(xlim = date_range, ylim = c(-0.15, 0.15)) +
  scale_y_continuous(breaks = c(-0.1, 0, 0.1)) +
  
  xlab(NULL) + ylab("Effect") +
  
  theme(legend.position = "none",
        plot.subtitle = ggtext::element_markdown(),
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank()) +
  
  ggtitle(NULL, "Annual periodic variation <i>f</i><sub>ann</sub>(<i>t</i>)")


p_Rt_exp <- ggplot() +
  
  plot_common +
  
  fastribbon::stat_fastribbon(aes(x = date - 0.5, y = exp(log_R), alpha = after_stat(level)),
                              geom = "stepribbon",
                              fill = "#3182BD",
                              traj_data_backcast_long,
                              linewidth = 0.0) +
  
  fastribbon::stat_fastribbon(aes(x = date - 0.5, y = exp(log_R), alpha = after_stat(level)),
                              geom = "stepribbon",
                              fill = colour_fc,
                              traj_data_forecast_long,
                              linewidth = 0.0) +
  
  plot_common_above +
  
  geom_hline(yintercept = 1, linetype = "22", linewidth = 0.7) +
  
  scale_fill_brewer(palette = "Blues") +
  
  coord_cartesian(xlim = date_range, ylim = c(0.7, 1.3)) +
  scale_y_continuous(breaks = c(0.75, 1, 1.25)) +
  
  xlab(NULL) + ylab("<i>R</i>(<i>t</i>)") +
  
  theme(legend.position = "none",
        plot.subtitle = ggtext::element_markdown(),
        panel.grid.major.x = element_blank()) +
  
  ggtitle(NULL, "Instantaneous reproduction number <i>R</i>(<i>t</i>)")

p_Rt_exp

p_cases_long <- ggplot(trajectory_data |> filter(date >= date_range[1], date <= date_range[2])) +
  
  plot_common +
  
  plot_common_above +
  
  geom_line(aes(x = date, y = y),
            linewidth = 0.4,
            training_counts |> filter(date >= date_range[1])) +
  
  scale_fill_brewer(palette = "Blues") +
  
  coord_cartesian(xlim = date_range, ylim = c(0, 1000)) +
  
  scale_y_continuous(labels = scales::label_comma()) +
  
  xlab(NULL) + ylab("Incidence") +
  
  theme(legend.position = "none") +
  
  ggtitle(NULL, "Influenza case incidence — Victoria")

p_cases_long


p_Rt <- ((p_cases_long + theme(axis.text.x = element_blank())) / p_Rt_exp) + 
  plot_annotation("<b>A</b>. Conditioned case incidence and inferred <i>R</i>(<i>t</i>)",
                  theme = plot_theme_annotation)

p_decomp <- (p_hol / p_ann / (p_noise + theme(axis.text.x = element_text(margin = margin(t = 0.3, unit = "cm"), colour = "black", vjust = 1)))) + 
  plot_annotation("<b>B</b>. Decomposed log <i>R</i>(<i>t</i>) additive effects",
                  theme = plot_theme_annotation)

n_days_ahead <- 7 * 4

date_range_forecast <- c(origin_date - days(7 * 8), origin_date + days(n_days_ahead))

traj_data_backcast <- trajectory_data |> 
  filter(date >= date_range_forecast[1], date <= origin_date)

traj_data_forecast <- trajectory_data |> 
  filter(date > origin_date, date <= date_range_forecast)

p_cases <- ggplot() +
  
  plot_common +
  
  fastribbon::stat_fastribbon(aes(x = date, y = cases, alpha = after_stat(level)),
                              fill = "#3182BD",
                              traj_data_backcast,
                              linewidth = 0.0) +
  
  fastribbon::stat_fastribbon(aes(x = date, y = cases, alpha = after_stat(level)),
                              fill = colour_fc,
                              traj_data_forecast,
                              linewidth = 0.0) +
  
  scale_alpha_manual(values = c(0.2, 0.5, 1)) +
  
  plot_common_above +
  
  geom_line(aes(x = date, y = y),
            training_counts |> filter(date >= date_range_forecast[1])) +
  geom_point(aes(x = date, y = y), size = 0.9,
             training_counts |> filter(date >= date_range_forecast[1])) +
  
  geom_line(aes(x = notification_date, y = value),
            eval_counts |> filter(notification_date <= date_range_forecast[2])) +
  geom_point(aes(x = notification_date, y = value), size = 0.9,
             eval_counts |> filter(notification_date <= date_range_forecast[2])) +
  
  scale_fill_brewer(palette = "Blues") +
  
  coord_cartesian(xlim = date_range_forecast, ylim = c(0, 1500)) +
  
  scale_x_date(breaks = scales::breaks_width("14 days", offset = -3),
                labels = scales::label_date_short()) +
  
  scale_y_continuous(labels = scales::label_comma()) +
  
  xlab(NULL) + ylab("Incidence") +
  
  theme(legend.position = "none")

p_cases


p_forecast <- (p_cases) + 
  plot_annotation("<b>C</b>. Forecast case incidence",
                  theme = plot_theme_annotation) +
  
  theme(plot.margin = margin(t = 0.2, l = 0.2, b = 0.2, r = 0.2, unit = "cm"))

p_example <- (wrap_elements(full = p_Rt) | wrap_elements(full = p_decomp)) /
  wrap_elements(full = p_forecast) /
  plot_spacer() +
  plot_layout(heights = c(1.2, 1, 0.1))

p_example

ggsave(
  "paper_figures/f_1_example.png", p_example,
  width = 9, height = 7.5,
  bg = "white"
)

library(magick)

image_composite(
  image = image_read("paper_figures/f_1_example.png"),
  composite_image = image_read("paper_figures/f_1_example_overlay.png"),
  operator = "over"
) |> 
  image_write(path = "paper_figures/f_1_example.png", format = "png")



  