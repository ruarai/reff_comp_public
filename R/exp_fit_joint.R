


library(tidyverse)
library(tidybayes)
library(bayesplot)
library(patchwork)


source("R/functions.R")

training_dates <- get_training_dates()

round_id <- "2025-05-08"
combined_counts <- read_csv(str_c("../forecast-hub/processed-data/2025/", round_id, "/combined-count-", round_id, ".csv"))
# combined_counts_smoothed <- read_csv(str_c("../forecast-hub/processed-data/2025/", round_id, "/combined-count-smoothed-", round_id, ".csv"))

x_pathogen <- "RSV"

school_holidays <- get_school_holidays_long()
public_holidays <- get_public_holidays()


pathogen_counts <- combined_counts |> 
  filter(pathogen == "RSV", notification_date >= ymd("2023-01-01"),
         location %in% c("VIC", "NSW", "WA")) |> 
  select(date = notification_date, location, cases)


origin_dates <- pathogen_counts |> 
  group_by(location) |> 
  summarise(origin_date = max(date)) |> 
  mutate(end_date = max(origin_date) + days(365))

forecast_filler_dates <- origin_dates |> 
  rowwise() |> 
  mutate(date = list(seq(origin_date + days(1), end_date, by = "day"))) |> 
  unnest(date) |> 
  select(location, date)

training_counts <- bind_rows(
  pathogen_counts,
  forecast_filler_dates
) |> 
  
  left_join(school_holidays) |> 
  left_join(public_holidays) |> 
  mutate(school_holiday_type = replace_na(school_holiday_type, 0) + 1,
         
         public_holiday_type = match(holiday_name, unique(sort(holiday_name))),
         public_holiday_type = replace_na(public_holiday_type, 0) + 1,
         day_of_week = wday(date, week_start = 7),
         x = as.numeric(date - min(date)))

unique_locs <- sort(unique(training_counts$location))
unique_dates <- sort(unique(training_counts$date))

get_matrix <- function(data, value_col) {
  data |> 
    select(date, location, {{value_col}}) |> 
    pivot_wider(names_from = location, values_from = {{value_col}}) |> 
    select(all_of(unique_locs)) |> 
    as.matrix()
}

N_loc <- length(unique_locs)
N <- length(unique(training_counts$date))
N_public_holiday_types <- length(unique(training_counts$public_holiday_type)) - 1

# TODO improve this
N_data <- training_counts |> 
  group_by(location) |> 
  summarise(N_data = max(x[!is.na(cases)]) + 1) |> 
  arrange(location) |> 
  pull(N_data) |> 
  min()

gi <- get_gi("RSV")

data_list <- list(
  N = N, N_loc = N_loc,
  N_data = N_data,
  N_public_holiday_types = N_public_holiday_types,
  
  x = get_matrix(training_counts, x)[, 1],
  day_of_week = get_matrix(training_counts, day_of_week)[, 1],
  
  y = t(get_matrix(training_counts, cases)[1:N_data, ]),
  
  public_holiday_type = t(get_matrix(training_counts, public_holiday_type)),
  school_holiday_type = t(get_matrix(training_counts, school_holiday_type)),
  
  
  M_annual = 40,
  M_noise = calculate_M_matern(c = 1.25, length_scale = 10.0, x = unique(training_counts$x)),
  c_factor = 1.25,
  
  gi_reversed = rev(gi), N_gi = length(gi)
)

library(cmdstanr)
model <- cmdstan_model("src/gp_renewal.stan")



fit <- model$pathfinder(
  seed = 4, data = data_list,
  max_lbfgs_iters = 1000
)

fit <- model$sample(
  data = data_list,
  seed = 5,
  chains = 4,
  parallel_chains = 4,
  refresh = 25,
  iter_warmup = 250, iter_sampling = 250,
  init = 0.25
)


fit |> 
  spread_draws(phi[ix_loc]) |> 
  mutate(location = unique_locs[ix_loc]) |> 
  ggplot() +
  geom_point(aes(x = phi, y = location))


trajectory_data <- fit |> 
  spread_draws(log_R[ix_t, ix_loc], latent_cases_weekly[ix_t, ix_loc],
               f_annual[ix_t, ix_loc], f_noise[ix_t, ix_loc]) |> 
  mutate(date = unique_dates[ix_t],
         location = unique_locs[ix_loc])


trajectory_data |> 
  filter(.draw %% 50 == 0) |> 
  
  ggplot() +
  geom_line(aes(x = date, y = exp(log_R), group = .draw)) +
  
  facet_wrap(~location, ncol = 1) +
  
  theme_bw()


fit |> 
  spread_draws(school_holiday_Rt_effect_raw[school_holiday_type, ix_loc], school_holiday_asc_effect_raw[school_holiday_type, ix_loc]) |> 
  mutate(location = unique_locs[ix_loc]) |> 
  mutate(holiday_name = case_when(
    school_holiday_type == 1 ~ "Summer holidays",
    school_holiday_type == 2 ~ "Autumn holidays",
    school_holiday_type == 3 ~ "Winter holidays",
    school_holiday_type == 4 ~ "Spring holidays"
  )) |> 
  ungroup() |> 
  mutate(holiday_name = factor(holiday_name, c("Summer holidays", "Autumn holidays", "Winter holidays", "Spring holidays"))) |> 
  
  ggplot() +
  ggdist::stat_histinterval(aes(y = holiday_name, x = school_holiday_Rt_effect_raw)) +
  
  facet_wrap(~location) +
  
  theme_bw()



fit |> 
  tidy_draws() |> 
  mutate(chain = .chain) |> 
  mcmc_trace(pars = str_c("delay_mean[",1:3,"]"))











