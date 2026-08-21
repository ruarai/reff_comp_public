


source("R/functions.R")

fit_files <- fs::dir_ls("results/2026-08-19/", recurse = TRUE, glob = "*.rds")

library(future)
library(furrr)
library(future.callr)

plan(callr, workers = 4)


Rt_draws <- future_map(
  fit_files,
  function(fit_file) {
    
    round_id <- ymd(str_extract(fit_file, "\\d{4}-\\d{2}-\\d{2}"))
    job_name <- fs::path_file(fit_file) |> fs::path_ext_remove()
    pathogen <- str_split_fixed(job_name, "_", 3)[1]
    location <- str_split_fixed(job_name, "_", 3)[2]
    target <- str_split_fixed(job_name, "_", 3)[3]
    
    
    fit <- read_rds(fit_file)
    training_counts <- read_parquet(glue::glue("results/2026-08-19/{pathogen}_{location}_{target}_training_counts.parquet"))
    
    
    trajectory_data <- fit |> 
      spread_draws(c(log_R,f_annual,f_noise)[ix_t],
                   log_R_0) |> 
      mutate(date = training_counts$date[ix_t],
             holiday_effect = log_R - log_R_0 - f_annual - f_noise) |> 
      mutate(location = !!location, pathogen = !!pathogen, target = !!target) |> 
      group_by(location, pathogen, target, date) |> 
      summarise(R = exp(mean(log_R)),
                R_holidayless = exp(mean(log_R - holiday_effect)))
  },
  .progress = TRUE
) |> 
  bind_rows()

write_parquet(Rt_draws, "data/paper/Rt_draws.parquet")

Rt_draws <- read_parquet("data/paper/Rt_draws.parquet")

combined_counts <- read_csv("../forecast-hub/processed-data/2026/2026-08-19/combined-count-smoothed-2026-08-19.csv") |> 
  rename(date = notification_date)


get_pop_size <- function(location) { # 2025 ABS ERP
  case_when(
    location == "VIC" ~ 7.121 * 1e6,
    location == "NSW" ~ 8.641 * 1e6,
    location == "QLD" ~ 5.712 * 1e6,
    location == "SA" ~ 1.910 * 1e6,
    location == "WA" ~ 3.076 * 1e6,
    location == "TAS" ~ 0.579 * 1e6,
    location == "NT" ~ 0.267 * 1e6,
    location == "ACT" ~ 0.487 * 1e6
  )
}

plot_data <- combined_counts |> 
  expand_grid(lag = c(50, 100, 150, 200, 250, 300, 350, 400, 450, 500)) |> 
  group_by(location, pathogen, lag) |> 
  mutate(cumulative_cases = cumsum(value_trend) - lag(cumsum(value_trend), lag[1], default = 0)) |> 
  ungroup() |> 
  filter(pathogen == "RSV") |> 
  left_join(Rt_draws) |> 
  
  mutate(pop_size = get_pop_size(location)) |> 
  
  mutate(cumulative_frac = (cumulative_cases / pop_size),
         location = factor(location),
         log_R_holidayless = log(R_holidayless),
         day_of_year = yday(date)) |> 
  
  drop_na(log_R_holidayless)


ggplot(plot_data) +
  geom_path(aes(x = cumulative_frac, y = log_R_holidayless, colour = day_of_year)) +
  
  theme_bw() +
  
  scale_colour_viridis_c() +
  
  facet_grid(rows = vars(location), cols = vars(lag),
             scales = "free_x")

combined_counts |> 
  expand_grid(lag = seq(10, 500, by = 10)) |> 
  group_by(location, pathogen, lag) |> 
  mutate(cumulative_cases = cumsum(value_trend) - lag(cumsum(value_trend), lag[1], default = 0)) |> 
  ungroup() |> 
  left_join(Rt_draws) |> 
  
  mutate(pop_size = get_pop_size(location)) |> 
  
  mutate(cumulative_frac = (cumulative_cases / pop_size),
         location = factor(location),
         log_R_holidayless = log(R_holidayless),
         day_of_year = yday(date)) |> 
  
  drop_na(log_R_holidayless) |> 
  
  group_by(location, pathogen, lag) |> 
  summarise(r = cor(cumulative_frac, log_R_holidayless)) |> 
  
  ggplot() +
  geom_line(aes(x = lag, y = r, colour = pathogen)) +
  
  geom_vline(xintercept = 150, linetype = "44") +
  
  facet_wrap(~location) +
  
  theme_bw()




fit_data <- combined_counts |> 
  group_by(location, pathogen) |> 
  mutate(cumulative_cases = cumsum(value_trend) - lag(cumsum(value_trend), 70, default = 0)) |> 
  ungroup() |> 
  filter(pathogen == "SARSCOV2") |> 
  left_join(Rt_draws) |> 
  
  mutate(pop_size = get_pop_size(location)) |> 
  
  mutate(cumulative_frac = (cumulative_cases / pop_size),
         location = factor(location),
         log_R_holidayless = log(R_holidayless),
         day_of_year = yday(date)) |> 
  
  drop_na(log_R_holidayless)

ggplot(fit_data) +
  geom_point(aes(x = cumulative_frac, y = log_R_holidayless, colour = day_of_year)) +
  
  theme_bw() +
  
  scale_colour_viridis_c() +
  
  facet_wrap(~location)

library(mgcv)

fit <- gam(
  log_R_holidayless ~ s(day_of_year, bs = "cc", by = location) +
    cumulative_frac * location + location,
  
  data = fit_data
)

summary(fit)

plot(fit, all.terms = TRUE, pages = 1)

library(marginaleffects)

plot_predictions(fit, condition = c("cumulative_frac", "location")) +
  theme_bw()


lags <- seq(10, 500, by = 5)

var_explained <- map_dbl(
  lags,
  function(lag) {
    
    fit_data <- combined_counts |> 
      group_by(location, pathogen) |> 
      mutate(cumulative_cases = cumsum(value_trend) - lag(cumsum(value_trend), lag, default = 0)) |> 
      ungroup() |> 
      filter(pathogen == "SARSCOV2") |> 
      left_join(Rt_draws) |> 
      
      mutate(pop_size = get_pop_size(location)) |> 
      
      mutate(cumulative_frac = (cumulative_cases / pop_size),
             location = factor(location),
             log_R_holidayless = log(R_holidayless),
             day_of_year = yday(date)) |> 
      
      drop_na(log_R_holidayless)
    
    fit <- gam(
      log_R_holidayless ~ s(day_of_year, bs = "cc", by = location) +
        cumulative_frac * location + location,
      
      data = fit_data
    )
    
    summary(fit)$r.sq
  }
)

tibble(
  lag = lags, var_explained = var_explained
) |> 
  ggplot() +
  geom_line(aes(x = lag, y = var_explained)) +
  
  theme_bw()



