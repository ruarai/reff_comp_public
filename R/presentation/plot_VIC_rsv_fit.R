
source("R/functions.R")
source("R/presentation/plot_theme.R")

round_results <- "saved_results/baseline/2025-05-08/"

x_location <- "VIC"
x_pathogen <- "RSV"
x_target <- "case_incidence"

job_name <- str_c(x_pathogen, "_", x_location, "_", x_target)



fit <- read_rds("saved_results/baseline/2025-05-08/RSV_VIC_case_incidence.rds")
training_counts <- read_parquet(str_c(round_results, "/", job_name, "_training_counts.parquet"))


origin_date <- training_counts |> 
  drop_na(y) |> 
  pull(date) |> 
  max()

eval_counts <- read_csv("../forecast-hub/processed-data/2026/2026-08-19/combined-count-smoothed-2026-08-19.csv", show_col_types = FALSE) |> 
  filter(notification_date > origin_date, location == x_location, pathogen == x_pathogen)





school_holidays <- get_school_holidays() |> 
  filter(location == x_location, end_date >= ymd("2023-01-01"), start_date <= origin_date + 365)

public_holidays <- get_public_holidays() |> 
  filter(location == x_location, date >= ymd("2023-01-01"), date <= origin_date + 365)



trajectory_data <- fit |> 
  spread_draws(log_R[ix_t], inf[ix_t], inf_delayed[ix_t],
               cases[ix_t], f_annual[ix_t], f_noise[ix_t],
               latent_cases_log[ix_t]) |> 
  
  mutate(date = training_counts$date[ix_t],
         latent_cases = exp(latent_cases_log),
         
         cases = replace_na(cases, 0))


labels <- c(
  "f_annual" = "exp(f_annual)",
  "f_noise" = "exp(f_noise)",
  "log_R" = "Rt"
)
alpha_vals <- scales::rescale(rev(1/1.7^(1:8)), to = c(0.05, 0.99))

trajectory_quants <- trajectory_data |> 
  make_intervals(c("f_annual", "f_noise", "log_R")) |> 
  mutate(variable = factor(variable, names(labels), labels = labels))

xlimits <- c(ymd("2024-01-01"), ymd("2026-01-01"))

p_reff <- ggplot(trajectory_quants) +
  
  geom_vline(aes(xintercept = date), tibble(date = seq(ymd("2023-01-01"), ymd("2027-01-01"), "1 year")),
             alpha = 0.3) +
  
  geom_rect(aes(xmin = start_date, xmax = end_date, ymin = 0, ymax = Inf),
            school_holidays, alpha = 0.1) +
  
  pammtools::geom_stepribbon(aes(x = date - 0.5, ymin = exp(lower), ymax = exp(upper),
                                 group = quant, alpha = quant, fill = variable)) +
  
  geom_hline(yintercept = 1, alpha = 0.5) +
  
  facet_grid(rows = vars(variable)) +
  
  xlab(NULL) + ylab("Value") +
  
  scale_alpha_manual(values = alpha_vals) +
  
  ggokabeito::scale_colour_okabe_ito(order = c(3, 5, 6)) +
  ggokabeito::scale_fill_okabe_ito(order = c(3, 5, 6)) +
  
  scale_y_log10(breaks = c(0.66, 0.8, 1, 1.25, 1.5)) +
  # scale_y_log10(breaks = c(0.8, 1, 1.25), labels = c("-20%", "0%", "+25%")) +
  
  coord_cartesian(ylim = c(1/1.4, 1.4), xlim = xlimits) +
  
  scale_x_date(labels = scales::label_date_short()) +
  
  geom_vline(xintercept = origin_date, linetype = "44") +
  
  plot_theme_pres + 
  theme(legend.position = "none") +
  
  ggtitle("RSV case incidence — Victoria, Australia")


y_upper_latent <- trajectory_data |> 
  filter(date <= origin_date + days(26 * 7)) |> 
  group_by(date) |> summarise(q95 = quantile(latent_cases, 0.95)) |> summarise(q95 = max(q95)) |> pull(q95)

plot_quants <- make_intervals(trajectory_data, "inf")
plot_quants_delayed <- make_intervals(trajectory_data, "latent_cases")


alpha_vals <- scales::rescale(rev(1/1.7^(1:8)), to = c(0.05, 0.99))


trajectory_trace <- trajectory_data |> 
  ungroup() |> 
  filter(.draw %% 2000 == 0) |>
  select(date, .draw, inf)


p_latent_inf <- ggplot() +
  
  geom_vline(aes(xintercept = date), tibble(date = seq(ymd("2023-01-01"), ymd("2027-01-01"), "1 year")),
             alpha = 0.4) +
  
  geom_rect(aes(xmin = start_date, xmax = end_date, ymin = 0, ymax = Inf),
            school_holidays, alpha = 0.1) +
  
  # pammtools::geom_stepribbon(aes(x = date - 0.5, ymin = lower, ymax = upper,
  #                                group = quant, alpha = quant),
  #                            plot_quants,
  #                            fill = ggokabeito::palette_okabe_ito(1)) +
  
  pammtools::geom_stepribbon(aes(x = date - 0.5, ymin = lower, ymax = upper,
                                 group = quant, alpha = quant),
                             plot_quants_delayed,
                             fill = ggokabeito::palette_okabe_ito(5)) +
  
  geom_line(aes(x = date, y = y_trend_and_remainder),
            training_counts) +
  
  geom_vline(xintercept = origin_date, linetype = "44") +
  
  coord_cartesian(ylim = c(0, y_upper_latent * 1.25),
                  xlim = xlimits) +
  
  scale_x_date(labels = scales::label_date_short()) +
  
  geom_line(aes(x = notification_date, y = value_trend_and_remainder), eval_counts, colour = "orangered") +
  
  xlab("Notification date") + ylab("Count") +
  
  scale_alpha_manual(values = alpha_vals) +
  
  plot_theme_pres + 
  
  theme(strip.text = element_blank(),
        legend.position = "none")


p_reff / p_latent_inf









