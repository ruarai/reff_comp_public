

source("R/functions.R")
source("R/presentation/plot_theme.R")



results <- read_rds("results_baseline/2025-10-09/flu_VIC.rds")


p_Rt <- results |> 
  spread_draws(school_holiday_Rt_effect_raw[i]) |> 
  
  mutate(school_holiday_type = case_when(
    i == 1 ~ "Summer",
    i == 2 ~ "Between-term"
  )) |> 
  
  ggplot() +
  
  ggdist::stat_histinterval(aes(x = exp(school_holiday_Rt_effect_raw), y = school_holiday_type)) +
  
  geom_vline(xintercept = 1, linetype = "44", linewidth = 1) +
  
  scale_x_log10(breaks = c(0.8, 0.9, 1, 1.1)) +
  
  ylab(NULL) + xlab("Multiplicative effect on R(t)") +
  
  plot_theme_pres

p_asc <- results |> 
  spread_draws(school_holiday_asc_effect_raw[i]) |> 
  
  mutate(school_holiday_type = case_when(
    i == 1 ~ "Summer",
    i == 2 ~ "Between-term"
  )) |> 
  
  ggplot() +
  
  ggdist::stat_histinterval(aes(x = exp(school_holiday_asc_effect_raw), y = school_holiday_type)) +
  
  geom_vline(xintercept = 1, linetype = "44", linewidth = 1) +
  
  scale_x_log10(breaks = c(0.6, 0.7, 0.8, 0.9, 1, 1.1)) +
  
  ylab(NULL) + xlab("Multiplicative effect on reporting") +
  
  plot_theme_pres

p_Rt / p_asc

results |> 
  spread_draws(school_holiday_Rt_effect_raw[i], school_holiday_asc_effect_raw[i]) |> 
  
  mutate(school_holiday_type = case_when(
    i == 1 ~ "Summer",
    i == 2 ~ "Between-term"
  )) |> 
  
  ggplot() +
  
  geom_density2d_filled(aes(x = exp(school_holiday_Rt_effect_raw), y = exp(school_holiday_asc_effect_raw), group = i),
                        contour_var = "ndensity", ) +
  
  geom_vline(xintercept = 1, linetype = "44", linewidth = 1) +
  geom_hline(yintercept = 1, linetype = "44", linewidth = 1) +
  
  scale_fill_manual(
    values = c("white", colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))(10))
  ) +
  
  facet_wrap(~school_holiday_type, ncol = 1) +
  
  xlab("Multiplicative effect on R(t)") + ylab("Multiplicative effect on reporting") +
  
  plot_theme_pres +
  
  theme(legend.position = "none",
        panel.grid.major = element_blank())


public_holidays <- get_public_holidays() |> 
  filter(location == "VIC") |> 
  mutate(public_holiday_type = match(holiday_name, unique(holiday_name))) |> 
  distinct(public_holiday_type, holiday_name)

p_asc <- results |> 
  spread_draws(public_holiday_asc_effect_raw[public_holiday_type]) |> 
  
  left_join(public_holidays) |> 
  
  ggplot() +
  
  ggdist::stat_histinterval(aes(x = exp(public_holiday_asc_effect_raw), y = holiday_name)) +
  
  geom_vline(xintercept = 1, linetype = "44", linewidth = 1) +
  
  scale_x_log10(breaks = c(0.3, 0.5, 0.75, 1, 1.3)) +
  
  ylab(NULL) + xlab("Multiplicative effect on reporting") +
  
  plot_theme_pres

p_Rt <- results |> 
  spread_draws(public_holiday_Rt_effect_raw[public_holiday_type]) |> 
  
  left_join(public_holidays) |> 
  
  ggplot() +
  
  ggdist::stat_histinterval(aes(x = exp(public_holiday_Rt_effect_raw), y = holiday_name)) +
  
  geom_vline(xintercept = 1, linetype = "44", linewidth = 1) +
  
  scale_x_log10(breaks = scales::breaks_log(n = 6)) +
  
  ylab(NULL) + xlab("Multiplicative effect on R(t)") +
  
  plot_theme_pres

p_asc | p_Rt

round_results <- "results_baseline/2025-10-09"
x_pathogen <- "flu"
x_location <- "VIC"
training_counts <- arrow::read_parquet(str_c(round_results, "/", x_pathogen, "_", x_location, "_training_counts.parquet"))


trajectory_data <- results |> 
  spread_draws(log_R[ix_t], inf[ix_t], inf_delayed[ix_t],
               cases[ix_t], f_annual[ix_t], f_noise[ix_t],
               latent_cases_log[ix_t], log_R_0) |> 
  mutate(date = training_counts$date[ix_t],
         latent_cases = exp(latent_cases_log))


date_range <- c(ymd("2024-06-01"), ymd("2024-08-01"))
y_range <- c(-0.25, 0.25)

library(pammtools)

school_holidays <- get_school_holidays() |> 
  filter(location == x_location)

public_holidays <- get_public_holidays() |> 
  filter(location == "VIC")

p_Rt <- ggplot(trajectory_data) +
  
  geom_rect(aes(xmin = start_date, xmax = end_date + 1, ymin = -Inf, ymax = Inf),
            school_holidays, alpha = 0.1)  +
  
  geom_rect(aes(xmin = date, xmax = date + 1, ymin = -Inf, ymax = Inf),
            public_holidays, alpha = 0.1)  +
  
  fastribbon::stat_fastribbon(aes(x = date, y = log_R), geom = "stepribbon") +
  
  geom_hline(yintercept = 0, linetype = "44") +
  
  scale_fill_brewer() +
  
  coord_cartesian(xlim = date_range, ylim = y_range) +
  
  scale_x_date(breaks = scales::date_breaks(),
               labels = scales::label_date_short()) +
  
  xlab(NULL) + ylab("log R(t)") +
  
  plot_theme_pres +
  
  theme(legend.position = "none") +
  
  ggtitle(NULL, "Composed log R(<i>t</i>)")  +
  
  plot_theme_pres +
  
  theme(legend.position = "none",
        plot.subtitle = ggtext::element_markdown())

p_noise <- ggplot(trajectory_data) +
  
  fastribbon::stat_fastribbon(aes(x = date, y = f_noise), geom = "stepribbon") +
  
  geom_hline(yintercept = 0, linetype = "44") +
  
  scale_fill_brewer() +
  
  coord_cartesian(xlim = date_range, ylim = y_range) +
  
  scale_x_date(breaks = scales::date_breaks(),
               labels = scales::label_date_short()) +
  
  xlab(NULL) + ylab("Additive effect") +
  
  ggtitle(NULL, "Noise <i>f</i><sub>noise</sub>(<i>t</i>)")  +
  
  plot_theme_pres +
  
  theme(legend.position = "none",
        plot.subtitle = ggtext::element_markdown())

p_ann <- ggplot(trajectory_data) +
  
  fastribbon::stat_fastribbon(aes(x = date, y = f_annual), geom = "stepribbon") +
  
  geom_hline(yintercept = 0, linetype = "44") +
  
  scale_fill_brewer() +
  
  coord_cartesian(xlim = date_range, ylim = y_range) +
  
  scale_x_date(breaks = scales::date_breaks(),
               labels = scales::label_date_short()) +
  
  xlab(NULL) + ylab("Additive effect") +
  
  plot_theme_pres +
  
  theme(legend.position = "none") +
  
  ggtitle(NULL, "Annual variation <i>f</i><sub>ann</sub>(<i>t</i>)")  +
  
  plot_theme_pres +
  
  theme(legend.position = "none",
        plot.subtitle = ggtext::element_markdown())

p_holiday <- ggplot(trajectory_data) +
  
  geom_rect(aes(xmin = date, xmax = date + 1, ymin = -Inf, ymax = Inf),
            public_holidays, alpha = 0.1)  +
  
  geom_rect(aes(xmin = start_date, xmax = end_date + 1, ymin = -Inf, ymax = Inf),
            school_holidays, alpha = 0.1)  +
  
  geom_rect(aes(xmin = date, xmax = date + 1, ymin = -Inf, ymax = Inf),
            public_holidays, alpha = 0.1)  +
  
  fastribbon::stat_fastribbon(aes(x = date, y = log_R - f_noise - f_annual - log_R_0), geom = "stepribbon") +
  
  geom_hline(yintercept = 0, linetype = "44") +
  
  scale_fill_brewer() +
  
  coord_cartesian(xlim = date_range, ylim = y_range) +
  
  scale_x_date(breaks = scales::date_breaks(),
               labels = scales::label_date_short()) +
  
  xlab(NULL) + ylab("Additive effect") +
  
  plot_theme_pres +
  
  theme(legend.position = "none") +
  
  ggtitle(NULL, "Holiday effects")  +
  
  plot_theme_pres +
  
  theme(legend.position = "none",
        plot.subtitle = ggtext::element_markdown())

(p_holiday | p_noise | p_ann)

p_Rt


p_inf <- ggplot(trajectory_data) +
  
  geom_rect(aes(xmin = start_date, xmax = end_date + 1, ymin = -Inf, ymax = Inf),
            school_holidays, alpha = 0.1)  +
  
  geom_rect(aes(xmin = date, xmax = date + 1, ymin = -Inf, ymax = Inf),
            public_holidays, alpha = 0.1)  +
  
  fastribbon::stat_fastribbon(aes(x = date, y = inf), geom = "stepribbon") +
  
  scale_fill_brewer(palette = "Oranges") +
  
  coord_cartesian(xlim = date_range, ylim = c(0, 1100)) +
  
  scale_x_date(breaks = scales::date_breaks(),
               labels = scales::label_date_short()) +
  
  xlab(NULL) + ylab("Incidence") +
  
  plot_theme_pres +
  
  theme(legend.position = "none") +
  
  ggtitle(NULL, "Latent infection incidence <i>z</i>(<i>t</i>)")  +
  
  plot_theme_pres +
  
  theme(legend.position = "none",
        plot.subtitle = ggtext::element_markdown())

p_delayed_inf <- ggplot(trajectory_data) +
  
  geom_rect(aes(xmin = start_date, xmax = end_date + 1, ymin = -Inf, ymax = Inf),
            school_holidays, alpha = 0.1)  +
  
  geom_rect(aes(xmin = date, xmax = date + 1, ymin = -Inf, ymax = Inf),
            public_holidays, alpha = 0.1)  +
  
  fastribbon::stat_fastribbon(aes(x = date, y = inf_delayed), geom = "stepribbon") +
  
  scale_fill_brewer(palette = "Oranges") +
  
  coord_cartesian(xlim = date_range, ylim = c(0, 1100)) +
  
  scale_x_date(breaks = scales::date_breaks(),
               labels = scales::label_date_short()) +
  
  xlab(NULL) + ylab("Incidence") +
  
  plot_theme_pres +
  
  theme(legend.position = "none") +
  
  ggtitle(NULL, "Delayed infection incidence")

p_latent_case <- ggplot(trajectory_data) +
  
  geom_rect(aes(xmin = start_date, xmax = end_date + 1, ymin = -Inf, ymax = Inf),
            school_holidays, alpha = 0.1)  +
  
  geom_rect(aes(xmin = date, xmax = date + 1, ymin = -Inf, ymax = Inf),
            public_holidays, alpha = 0.1)  +
  
  fastribbon::stat_fastribbon(aes(x = date, y = latent_cases), geom = "stepribbon") +
  
  scale_fill_brewer(palette = "Greens") +
  
  coord_cartesian(xlim = date_range, ylim = c(0, 1100)) +
  
  scale_x_date(breaks = scales::date_breaks(),
               labels = scales::label_date_short()) +
  
  xlab(NULL) + ylab("Incidence") +
  
  ggtitle(NULL, "Latent case incidence <i>y</i>(<i>t</i>)")  +
  
  plot_theme_pres +
  
  theme(legend.position = "none",
        plot.subtitle = ggtext::element_markdown())


p_case <- ggplot(trajectory_data) +
  
  geom_rect(aes(xmin = start_date, xmax = end_date + 1, ymin = -Inf, ymax = Inf),
            school_holidays, alpha = 0.1) +
  
  geom_rect(aes(xmin = date, xmax = date + 1, ymin = -Inf, ymax = Inf),
            public_holidays, alpha = 0.1)  +
  
  fastribbon::stat_fastribbon(aes(x = date, y = cases), geom = "stepribbon") +
  
  geom_line(aes(x = date + 0.5, y = y),
            training_counts) +
  
  geom_point(aes(x = date + 0.5, y = y),
             training_counts) +
  
  scale_fill_brewer(palette = "Greens") +
  
  coord_cartesian(xlim = date_range, ylim = c(0, 1600)) +
  
  scale_x_date(breaks = scales::date_breaks(),
               labels = scales::label_date_short()) +
  
  xlab(NULL) + ylab("Incidence") +
  
  plot_theme_pres +
  
  theme(legend.position = "none") +
  
  ggtitle(NULL, "Reported case incidence")

patchwork::wrap_plots(p_Rt, p_inf, p_latent_case, p_case, ncol = 2)

