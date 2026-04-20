
source("R/functions.R")

fit_files <- fs::dir_ls("saved_results/results_baseline/2025-10-09/", recurse = TRUE, glob = "*.rds")

fit_file <- fit_files[[1]]

gp_draws <- map(
  fit_files,
  function(fit_file) {
    
    fit <- read_rds(fit_file)
    
    round_id <- ymd(str_extract(fit_file, "\\d{4}-\\d{2}-\\d{2}"))
    pathogen_location <- fs::path_file(fit_file) |> fs::path_ext_remove()
    pathogen <- str_split_fixed(pathogen_location, "_", 2)[1]
    location <- str_split_fixed(pathogen_location, "_", 2)[2]
    
    training_counts <- read_parquet(
      str_c("saved_results/results_baseline/", round_id, "/", pathogen_location, "_training_counts.parquet")
    )
    
    fit |> 
      
      spread_draws(f_annual[ix_t], f_noise[ix_t]) |> 
      filter(.draw %% 10 == 0) |> 
      mutate(date = training_counts$date[ix_t]) |> 
      mutate(round_id = round_id, location = location, pathogen = pathogen)
  },
  .progress = TRUE
) |> 
  bind_rows()

ggplot(gp_draws |> filter(location != "NZ", year(date) == 2024, location != "SA")) +
  fastribbon::stat_fastribbon(aes(x = date, y = f_annual),
                              linewidth = 0.0) +
  
  geom_hline(yintercept = 0, linetype = "22") +
  
  facet_grid(rows = vars(location), cols = vars(pathogen)) +
  
  scale_y_continuous(breaks = c(-0.2, 0, 0.2)) +
  
  scale_x_date(labels = scales::label_date_short(format = c("", "%b"))) +
  
  scale_fill_brewer() +
  
  plot_theme_pres_2 + xlab(NULL) + ylab(NULL) +
  theme(legend.position = "none")



ggplot(gp_draws |> filter(location != "NZ", year(date) == 2024, location != "SA", pathogen == "RSV")) +
  fastribbon::stat_fastribbon(aes(x = date, y = f_annual),
                              linewidth = 0.0) +
  
  geom_hline(yintercept = 0, linetype = "22") +
  
  facet_grid(cols = vars(location)) +
  
  scale_x_date(labels = scales::label_date_short()) +
  
  scale_fill_brewer() +
  
  theme_bw() +
  theme(legend.position = "none")

aus_humidity <- read_csv("data/australia_humidity.csv")

aus_humidity |> 
  `colnames<-`(c("month", "NSW", "VIC", "QLD", "WA", "SA", "TAS", "ACT", "NT")) |> 
  pivot_longer(-month) |> 
  filter(name != "SA") |> 
  mutate(date = ym(str_c("2025-",month))) |> 
  
  ggplot() +
  
  geom_line(aes(x = date, y = value)) +
  
  facet_grid(cols = vars(name))

weekly_humidity <- daily_humidity |> 
  mutate(date = floor_date(date, "week")) |> 
  group_by(state, date) |> 
  summarise(rel_hum = mean(rel_hum, na.rm = TRUE))

gp_draws |> 
  group_by(location, pathogen, date) |> 
  summarise(f = median(f_annual + f_noise)) |> 
  
  left_join(daily_humidity |> rename(location = state)) |> 
  
  ggplot() +
  geom_point(aes(x = rel_hum, y = exp(f) )) +
  facet_grid(rows = vars(location), cols = vars(pathogen))

gp_draws |> 
  
  ggplot() +
  
  geom_line(aes(x = date, y = rel_hum),
            weekly_humidity |> rename(location = state)) +
  fastribbon::stat_fastribbon(aes(x = date, y = (f_annual + f_noise) * 100 + 30),
                              linewidth = 0.0) +
  
  scale_fill_brewer() +
  facet_grid(rows = vars(location), cols = vars(pathogen)) +
  
  theme_bw()
