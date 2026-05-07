
source("R/functions.R")


round_id <- "2025-07-10"

x_pathogen <- "RSV"
x_target <- "case_incidence"

year <- as.character(year(round_id))
combined_counts <- read_csv(str_c("../forecast-hub/processed-data/", year, "/", round_id, "/combined-count-smoothed-", round_id, ".csv")) |> 
  filter(location %in% c("ACT", "NSW", "WA", "QLD", "VIC", "TAS"))

source("R/exp/process_data_hierarchical.R")

origin_dates <- combined_counts |> 
  filter(pathogen == x_pathogen, target == x_target,
         is_complete) |> 
  group_by(location) |> 
  filter(notification_date == max(notification_date)) |> 
  arrange(location) |> 
  select(location, origin_date = notification_date) |> 
  ungroup()

n_days_ahead <- as.numeric(ymd("2027-01-01") - max(origin_dates$origin_date)) - 1

training_counts <- process_counts_hierarchical(combined_counts, x_pathogen, x_target, n_days_ahead, origin_dates) |> 
  mutate(year = match(year(date), unique(year(date))))



year_ranges <- training_counts |> 
  ungroup() |> filter(location == location[[1]]) |> 
  mutate(ix = row_number()) |>
  group_by(year) |> 
  summarise(ix_year_start = min(ix), ix_year_end = max(ix))


ggplot(training_counts) +
  geom_tile(aes(x = date, y = location, fill = factor(year)))


data_list <- create_stan_data_hierarchical(training_counts, n_days_ahead, x_pathogen)
data_list$N_years <- max(training_counts$year)
data_list$year_start <- year_ranges$ix_year_start
data_list$year_end <- year_ranges$ix_year_end
data_list$sigma_log10_total_inf <- 0.1

library(cmdstanr)
model <- cmdstan_model("src/gp_renewal_hierarchical.stan")

fit <- model$sample(
  data = data_list,
  seed = 1,
  chains = 4,
  parallel_chains = 4,
  init = 0.01,
  refresh = 50,
  iter_warmup = 500, iter_sampling = 500
)

fit$save_object("saved_results/hierarchical_exp/fit_RSV_2.rds")

q('no')

fit <- read_rds("saved_results/hierarchical_exp/fit_RSV_2.rds")


draws <- tidy_draws(fit) |> mutate(chain = .chain)

date_index <- unique(training_counts$date)
location_index <- unique(training_counts$location)


trajectory_data <- fit |> 
  spread_draws(c(inf, inf_delayed, latent_cases_weekly_log,
                 cases)[i, ix_t]) |> 
  mutate(date = date_index[ix_t],
         location = location_index[i])

ggplot(trajectory_data |> filter(date >= ymd("2025-01-01"))) + 
  fastribbon::stat_fastribbon(aes(x = date, y = inf_delayed),
                              linewidth = 0.0) +
  
  
  geom_line(aes(x = notification_date, y = value_trend),
            combined_counts |> 
              filter(pathogen == "RSV",
                     location %in% location_index,
                     notification_date >= min(origin_dates$origin_date)),
            colour = "orangered") +
  
  geom_line(aes(x = date, y = y_trend),
            training_counts |> 
              filter(pathogen == "RSV")) +
  
  coord_cartesian(xlim = c(ymd("2025-01-01"), NA)) +
  
  facet_wrap(~location, scales = "free_y") +
  
  scale_fill_brewer() + theme_bw() + theme(legend.position = "none")



year_data <- fit |> 
  spread_draws(log10_total_inf_year[i, ix_year])

ggplot(year_data) + 
  fastribbon::stat_fastribbon(aes(x = ix_year, y = log10_total_inf_year),
                              linewidth = 0.0) +
  
  facet_wrap(~i, scales = "free_y") +
  
  scale_fill_brewer() + theme_bw() + theme(legend.position = "none")


fit |> 
  spread_draws(log10_total_inf_offset[i]) |> 
  filter(i > 1) |> 
  ggplot() +
  geom_histogram(aes(x = log10_total_inf_offset)) +
  
  facet_wrap(~i, scales = "free_y") +
  scale_fill_brewer() + theme_bw() + theme(legend.position = "none")

fit |> 
  spread_draws(log10_total_inf_year_mean[ix_year], log10_total_inf_offset[i]) |> 
  mutate(latent_log10_total_inf_year = log10_total_inf_year_mean + log10_total_inf_offset) |> 

  ggplot() + 
  fastribbon::stat_fastribbon(aes(x = ix_year, y = latent_log10_total_inf_year),
                              linewidth = 0.0) + 
  
  facet_wrap(~i) +
  
  scale_fill_brewer() + theme_bw() + theme(legend.position = "none")


combined_counts <- read_csv("data/combined-count-smoothed-2026-04-22.csv")

plot_data_true_total <- combined_counts |> 
  filter(location %in% location_index, pathogen == x_pathogen,
         notification_date >= ymd("2023-01-01")) |> 
  group_by(location, year = floor_date(notification_date, "year")) |> 
  summarise(total_inf = sum(value))

trajectory_data |> 
  filter(.draw %% 10 == 0) |> 
  group_by(location, year = floor_date(date, "year"), .draw) |> 
  summarise(total_inf = sum(inf)) |> 
  ggplot() + 
  geom_point(aes(x = year, y = total_inf),
             position = position_jitter(width = 50),
             size = 0.5, stroke = 0.5, alpha = 0.5) +
  
  geom_point(aes(x = year, y = total_inf),
             colour = "orangered",
             plot_data_true_total) +
  facet_wrap(~location, scales = "free_y") +
  
  theme_bw()
