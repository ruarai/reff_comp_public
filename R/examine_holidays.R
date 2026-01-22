
x_location <- "WA"
x_pathogen <- "flu"


round_results <- fs::path("results", round_id)

fit <- read_rds(str_c(round_results, "/", x_pathogen, "_", x_location, ".rds"))
training_counts <- arrow::read_parquet(str_c(round_results, "/", x_pathogen, "_", x_location, "_training_counts.parquet"))


origin_date <- training_counts |> 
  drop_na(y) |> 
  pull(date) |> 
  max()

eval_counts <- read_csv("data/combined-count-2026-01-15.csv") |> 
  group_by(location, pathogen) |> 
  mutate(cases_trend_and_remainder = zoo::rollmean(cases, k = 7, fill = NA),
         outlier = FALSE) |> 
  filter(location == x_location, pathogen == x_pathogen) |> 
  make_training_counts(x_location, x_pathogen) |> 
  filter(date > origin_date)


school_holidays <- get_school_holidays() |> 
  filter(location == x_location)

public_holidays <- get_public_holidays() |> 
  filter(location == x_location)

trajectory_data <- fit |> 
  spread_draws(log_R[ix_t], inf[ix_t], inf_delayed[ix_t],
               cases[ix_t], f_annual[ix_t], f_noise[ix_t],
               latent_cases[ix_t]) |> 
  mutate(date = training_counts$date[ix_t])

holiday_period <- c(ymd("2024-12-01"), ymd("2025-02-07"))
holiday_period <- c(ymd("2023-12-01"), ymd("2024-02-07"))
y_upper <- 100

p_reff <- plot_reff_components(trajectory_data, origin_date, school_holidays, public_holidays) +
  xlim(holiday_period)

p_latent_inf <- plot_latent_inf(trajectory_data, training_counts, origin_date, school_holidays, public_holidays) +
  
  geom_line(aes(x = date, y = y_smooth), colour = "orangered",
            eval_counts) +
  xlim(holiday_period) +
  coord_cartesian(ylim = c(0, y_upper))

p_cases <- plot_cases(trajectory_data, training_counts, origin_date, school_holidays, public_holidays) +
  xlim(holiday_period) +
  coord_cartesian(ylim = c(0, y_upper))

p_reff / p_latent_inf / p_cases

