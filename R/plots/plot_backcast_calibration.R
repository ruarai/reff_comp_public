
source("R/functions.R")

fit_files <- fs::dir_ls("results/2026-04-22/", recurse = TRUE, glob = "*.rds")

backcast_draws <- map(
  fit_files,
  function(fit_file) {
    
    fit <- read_rds(fit_file)
    
    round_id <- ymd(str_extract(fit_file, "\\d{4}-\\d{2}-\\d{2}"))
    pathogen_location <- fs::path_file(fit_file) |> fs::path_ext_remove()
    pathogen <- str_split_fixed(pathogen_location, "_", 3)[1]
    location <- str_split_fixed(pathogen_location, "_", 3)[2]
    
    training_counts <- read_parquet(
      str_c("results/", round_id, "/", pathogen_location, "_training_counts.parquet")
    )
    
    fit |> 
      
      spread_draws(latent_cases_log[ix_t]) |> 
      filter(.draw %% 10 == 0) |> 
      mutate(date = training_counts$date[ix_t]) |> 
      mutate(round_id = round_id, location = location, pathogen = pathogen,
             latent_cases = exp(latent_cases_log))
  },
  .progress = TRUE
) |> 
  bind_rows()



forecast_quants <- backcast_draws |> 
  mutate(count = round(latent_cases)) |> 
  select(location, pathogen, notification_date = date, sample = .draw, count) |> 
  expand_grid(q = seq(0.0, 1.0, by = 0.1)) |> 
  group_by(location, pathogen, notification_date, q) |> 
  summarise(count_q = quantile(count, q[1]))


combined_counts <- read_csv("../forecast-hub/processed-data/2026/2026-04-22/combined-count-smoothed-2026-04-22.csv")


forecast_quants |> 
  filter(location != "NZ", location != "SA") |> 
  left_join(combined_counts, by = c("location", "pathogen", "notification_date")) |> 
  mutate(y = value_trend_and_remainder) |> 
  drop_na(y) |> 
  mutate(below_q = y <= count_q) |> 
  group_by(location, pathogen, q) |> 
  summarise(n_total = n(), n = sum(below_q)) |> 
  group_by(location, pathogen) |> 
  mutate(n = n / n_total) |> 
  
  ggplot() +
  geom_line(aes(x = q, y = n)) +
  
  geom_abline(intercept = 0, slope = 1, colour = ggokabeito::palette_okabe_ito(5)) +
  
  facet_grid(rows = vars(location), cols = vars(pathogen)) +
  
  coord_fixed(xlim = c(0,1), ylim = c(0,1)) +
  
  theme_bw() +
  
  ggtitle("Backcast calibration (uom-reffcomp)")










