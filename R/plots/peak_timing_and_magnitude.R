

library(scoringRules)

results_path <- "saved_results/results_neither"

jobs <- make_retrospective_jobs() |> 
  filter(round_id >= "2025-05-08", round_id <= "2025-10-09")

peak_scores <- jobs |> 
  rename(x_location = location, x_pathogen = pathogen) |> 
  pmap(function(round_id, x_location, x_pathogen) {
    
    fit <- read_rds(str_c(results_path, "/", round_id, "/", x_pathogen, "_", x_location, ".rds"))
    training_counts <- read_parquet(
      str_c(results_path, "/", round_id, "/", x_pathogen, "_", x_location, "_training_counts.parquet")
    )
    
    origin_date <- training_counts$origin_date[[1]]
    
    eval_counts <- read_csv("data/combined-count-smoothed-2026-01-15.csv", show_col_types = FALSE)
    
    case_trajectories <- fit |> 
      spread_draws(latent_cases_log[ix_t]) |> 
      mutate(latent_cases = exp(latent_cases_log),
             date = training_counts$date[ix_t])
    
    
    case_peaks <- case_trajectories |> 
      filter(date >= ymd("2025-01-01"), date < ymd("2026-01-01")) |> 
      group_by(.draw) |> 
      summarise(
        peak_latent_cases = max(latent_cases),
        peak_date = first(date[latent_cases == max(latent_cases)])
      )
    
    true_peak <- eval_counts |> 
      filter(location == x_location, pathogen == x_pathogen) |> 
      filter(notification_date >= ymd("2025-01-01"), notification_date < ymd("2026-01-01")) |> 
      summarise(
        peak_latent_cases = max(cases_trend),
        peak_date = first(notification_date[cases_trend == max(cases_trend)])
      )
    
    ggplot() +
      geom_point(aes(x = peak_date, y = peak_latent_cases),
                 case_peaks) +
      
      geom_point(aes(x = peak_date, y = peak_latent_cases),
                 colour = "red",
                 true_peak) +
      
      coord_cartesian(ylim = c(0, NA))
    
    
    true_peak_vec <- c(as.numeric(true_peak$peak_date), true_peak$peak_latent_cases)
    
    ens_matrix <- case_peaks %>%
      select(peak_date, peak_latent_cases) %>%
      mutate(peak_date = as.numeric(peak_date)) %>%
      as.matrix() |> 
      t()
    
    # These don't really make sense given the different units
    # for peak timing vs peak scale
    energy_score <- es_sample(y = true_peak_vec, dat = ens_matrix)
    variogram_score <- vs_sample(y = true_peak_vec, dat = ens_matrix)
    
    crps_date <- crps_sample(y = true_peak_vec[1], dat = ens_matrix[1, ])
    crps_cases <- crps_sample(y = true_peak_vec[2], dat = ens_matrix[2, ])
    crps_cases_log <- crps_sample(y = log(true_peak_vec[2] + 1), dat = log(ens_matrix[2, ] + 1))
    
    tibble(
      round_id,
      location = x_location,
      pathogen = x_pathogen,
      energy_score,
      variogram_score,
      crps_date,
      crps_cases,
      crps_cases_log
    )
  },
  .progress = TRUE) |> 
  bind_rows()

# write_parquet(peak_scores, str_c(results_path, "/peak_scores.parquet"))

peak_scores <- read_parquet("saved_results/results_neither/peak_scores.parquet")


peak_scores |> 
  filter(location %in% c("WA", "NSW", "VIC", "QLD", "ACT"),
         pathogen != "SARSCOV2") |> 
  mutate(round_id = ymd(round_id)) |> 
  group_by(location, pathogen) |> 
  mutate(y_date = crps_date / min(crps_date),
         y_cases = crps_cases_log / min(crps_cases_log)) |> 
  ggplot() +
  geom_point(aes(x = round_id, y = crps_date),
             colour = ggokabeito::palette_okabe_ito(5)) +
  geom_line(aes(x = round_id, y = crps_date),
            colour = ggokabeito::palette_okabe_ito(5)) +
  # geom_point(aes(x = round_id, y = y_cases)) +
  # geom_line(aes(x = round_id, y = y_cases)) +
  
  facet_grid(rows = vars(location),
             cols = vars(pathogen),
             scales = "free_y") +
  
  theme_bw()


peak_scores |> 
  filter(location %in% c("WA", "NSW", "VIC", "QLD", "ACT")) |> 
  mutate(round_id = ymd(round_id)) |> 
  ggplot() +
  geom_point(aes(x = round_id, y = energy_score),
             colour = ggokabeito::palette_okabe_ito(5)) +
  geom_line(aes(x = round_id, y = energy_score),
            colour = ggokabeito::palette_okabe_ito(5)) +
  # geom_point(aes(x = round_id, y = y_cases)) +
  # geom_line(aes(x = round_id, y = y_cases)) +
  
  facet_grid(rows = vars(location),
             cols = vars(pathogen),
             scales = "free_y") +
  
  theme_bw()


plot_data <- bind_rows(plot_data, peak_scores |> 
  # filter(pathogen != "SARSCOV2") |> 
  mutate(round_id = ymd(round_id)) |> 
  filter(round_id == min(round_id)) |> 
  mutate(.model = "neither"))

ggplot(plot_data) +
  geom_point(aes(x = crps_cases_log, y = crps_date, colour = .model)) +
  geom_text(aes(x = crps_cases_log, y = crps_date, label = location),
            position = position_nudge(y = -5)) +
  facet_wrap(~pathogen, ncol = 1) +
  
  coord_cartesian(xlim = c(0, NA), ylim = c(0, NA)) +
  
  theme_bw()



