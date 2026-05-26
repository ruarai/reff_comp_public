source("R/functions.R")
source("R/paper/common.R")

library(pammtools)

#==2025-04-24 NSW flu (grew too fast due to f_ann)==
#==2025-06-19 VIC RSV (confidently went down too soon due to f_ann?)==
#==2025-04-04 NSW RSV (pretty good, holidays + f_ann)==
#==2025-07-10 VIC SARSCOV2 (pretty good, informed by holidays + f_ann)==

eval_counts <- read_csv("data/combined-count-smoothed-2026-04-15.csv", show_col_types = FALSE)

round_id <- "2025-07-10"
x_pathogen <- "SARSCOV2"
x_location <- "VIC"
x_target <- "case_incidence"

job_name <- str_c(x_pathogen, "_", x_location, "_", x_target)

fit <- read_rds(str_c("saved_results/baseline/", round_id, "/", job_name, ".rds"))
training_counts <- read_parquet(str_c("saved_results/baseline/", round_id, "/", job_name, "_training_counts.parquet"))

origin_date <- training_counts$origin_date[[1]] - days(1)

trajectory_data <- fit |> 
  spread_draws(c(log_R, inf, inf_delayed, cases, f_annual,
                 f_noise, latent_cases_log)[ix_t],
               log_R_0) |> 
  mutate(date = training_counts$date[ix_t],
         latent_cases = exp(latent_cases_log),
         holiday_effect = log_R - log_R_0 - f_annual - f_noise)

plot_data <- trajectory_data |> 
  filter(date >= origin_date - days(7 * 8),
         date <= origin_date + days(7 * 4))

plot_data_eval <- eval_counts |> 
  filter(location == x_location, pathogen == x_pathogen) |> 
  filter(notification_date > origin_date,
         notification_date <= origin_date + days(7 * 4))

plot_data_training <- training_counts |> 
  filter(date >= origin_date - days(7 * 8)) |> 
  drop_na(y)


ggplot(plot_data) +
  fastribbon::stat_fastribbon(aes(x = date - 0.5, y = cases),
                              geom = "stepribbon") +
  
  geom_vline(xintercept = origin_date,
             linetype = "44") +
  
  geom_line(aes(x = date, y = y), plot_data_training) +
  geom_point(aes(x = date, y = y), plot_data_training) +
  
  geom_line(aes(x = notification_date, y = value),
            linetype = "41",
            plot_data_eval) +
  geom_point(aes(x = notification_date, y = value),
             plot_data_eval) +
  
  scale_fill_brewer() +
  xlab(NULL) + ylab("Count") +
  
  coord_cartesian(ylim = c(0, 400)) +
  
  plot_theme_paper +
  
  theme(legend.position = "none")





plot_data <- trajectory_data |> 
  filter(date >= origin_date,
         date <= origin_date + days(7 * 4))

plot_data_eval <- eval_counts |> 
  filter(location == x_location, pathogen == x_pathogen) |> 
  filter(notification_date >= origin_date,
         notification_date <= origin_date + days(7 * 4))

dates <- unique(plot_data$date)

plot_dates <- function(date_pair) {
  if(date_pair[1] >= date_pair[2]) {
    return(ggplot() + geom_blank() + theme_void())
  }
  
  plot_data_samples <- plot_data |> 
    filter(.draw %% 10 == 0) |> 
    filter(date %in% date_pair) |> 
    group_by(.draw) |> 
    mutate(cases_1 = latent_cases, cases_2 = lead(latent_cases)) |> 
    drop_na(cases_2)
  
  plot_data_obs <- plot_data_eval |> 
    filter(notification_date %in% date_pair) |> 
    mutate(cases_1 = value_trend, cases_2 = lead(value_trend))
  
  ggplot() +
    geom_point(aes(x = cases_1, y = cases_2),
               plot_data_samples,
               size = 0.5, stroke = 0.5, alpha = 0.25) +
    geom_point(aes(x = cases_1, y = cases_2),
               plot_data_obs,
               colour = "orangered") +
    xlab(NULL) + ylab(NULL) +
    theme_bw() +
    theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, unit = "cm"),
          panel.grid = element_blank())
}


plot_dates(dates[c(1, 28)])

plots <- expand_grid(
  date_1 = dates[c(1, 7, 14, 21, 29)],
  date_2 = dates[c(1, 7, 14, 21, 29)]
) |> 
  rowwise() |> 
  mutate(date_plot = list(plot_dates(c(date_1, date_2))))


wrap_plots(plots$date_plot)




