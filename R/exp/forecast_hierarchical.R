


source("R/functions.R")


round_id <- "2025-06-12"

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

fit$save_object("saved_results/hierarchical_exp/fit_RSV.rds")

q('no')

fit <- read_rds("saved_results/hierarchical_exp/fit_RSV.rds")


draws <- tidy_draws(fit) |> mutate(chain = .chain)

index_ref <- training_counts |> 
  group_by(location) |> 
  mutate(ix_t = row_number()) |> 
  ungroup() |> 
  mutate(i = match(location, unique(location))) |> 
  select(i, ix_t, location, date)

trajectory_data <- fit |> 
  spread_draws(c(inf, inf_delayed, latent_cases_weekly_log,
                 cases)[i, ix_t]) |> 
  left_join(index_ref)

ggplot(trajectory_data) + 
  fastribbon::stat_fastribbon(aes(x = date, y = inf_delayed),
                              linewidth = 0.0) +
  
  facet_wrap(~location, scales = "free_y") +
  
  scale_fill_brewer() + theme_bw() + theme(legend.position = "none")



year_data <- fit |> 
  spread_draws(log10_total_inf_year[i, ix_year]) |> 
  filter(ix_year < 5)

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
  filter(ix_year < 5) |> 

  ggplot() + 
  fastribbon::stat_fastribbon(aes(x = ix_year, y = latent_log10_total_inf_year),
                              linewidth = 0.0) + 
  fastribbon::stat_fastribbon(aes(x = ix_year, y = log10_total_inf_year),
                              year_data, alpha = 0.5,
                              linewidth = 0.0) +
  
  facet_wrap(~i, scales = "free_y") +
  
  scale_fill_brewer() + theme_bw() + theme(legend.position = "none")
