
source("R/functions.R")

fit_files <- fs::dir_ls("results/2026-04-29/", recurse = TRUE, glob = "*.rds")

gp_draws <- map(
  fit_files,
  function(fit_file) {
    
    fit <- read_rds(fit_file)
    
    round_id <- ymd(str_extract(fit_file, "\\d{4}-\\d{2}-\\d{2}"))
    job_name <- fs::path_file(fit_file) |> fs::path_ext_remove()
    pathogen <- str_split_fixed(job_name, "_", 3)[1]
    location <- str_split_fixed(job_name, "_", 3)[2]
    target <- str_split_fixed(job_name, "_", 3)[3]
    
    training_counts <- read_parquet(
      str_c("results/", round_id, "/", job_name, "_training_counts.parquet")
    )
    
    fit |> 
      spread_draws(f_annual[ix_t], f_noise[ix_t], log_R[ix_t]) |> 
      filter(.draw %% 10 == 0) |> 
      mutate(date = training_counts$date[ix_t]) |> 
      mutate(round_id = round_id, location = location, pathogen = pathogen, target = target)
  },
  .progress = TRUE
) |> 
  bind_rows()

plot_data <- gp_draws |> filter(location != "NZ", location != "SA", year(date) == 2025)

ggplot(plot_data) +
  fastribbon::stat_fastribbon(aes(x = date, y = f_annual),
                              linewidth = 0.0) +
  
  geom_hline(yintercept = 0, linetype = "22") +
  
  facet_grid(rows = vars(location), cols = vars(pathogen)) +
  
  scale_y_continuous(breaks = c(-0.2, 0, 0.2)) +
  
  scale_x_date(labels = scales::label_date_short(format = c("", "%b"))) +
  
  scale_fill_brewer() +
  
  theme_bw() + xlab(NULL) + ylab(NULL) +
  theme(legend.position = "none")

plot_data |> 
  filter(location %in% c("WA", "VIC", "NSW")) |> 
  group_by(location, pathogen, date) |> 
  mutate(y = f_noise) |> 
  summarise(m = median(y),
            q25 = quantile(y, 0.25), q75 = quantile(y, 0.75)) |> 

  ggplot() +
  geom_line(aes(x = date, y = m, colour = location)) +
  geom_ribbon(aes(x = date, ymin = q25, ymax = q75, fill = location),
              alpha = 0.2) +
  
  facet_wrap(~pathogen) +
  
  geom_hline(yintercept = 0, linetype = "22") +
  
  scale_y_continuous(breaks = c(-0.2, 0, 0.2)) +
  
  coord_cartesian(ylim = c(-0.3, 0.3)) +
  
  scale_x_date(labels = scales::label_date_short(format = c("", "%b"))) +
  
  theme_bw() + xlab(NULL) + ylab(NULL) +
  theme(legend.position = "bottom")


