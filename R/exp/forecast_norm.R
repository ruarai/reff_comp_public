


source("R/functions.R")


round_id <- "2026-04-29"

x_location <- "WA"
x_pathogen <- "flu"
x_target <- "case_incidence"

year <- as.character(year(round_id))
combined_counts <- read_csv(str_c("../forecast-hub/processed-data/", year, "/", round_id, "/combined-count-smoothed-", round_id, ".csv"))

origin_date <- get_origin_date(combined_counts, x_pathogen, x_location, x_target)


eval_counts <- read_csv("../forecast-hub/processed-data/2026/2026-04-22/combined-count-smoothed-2026-04-22.csv") |> 
  filter(notification_date > origin_date, location == x_location, pathogen == x_pathogen)

n_days_ahead <- as.numeric(ymd("2027-01-01") - origin_date)

training_counts <- process_counts(combined_counts, x_location, x_pathogen, x_target, n_days_ahead, origin_date)


data_list <- create_stan_data(training_counts, n_days_ahead, x_pathogen, x_location)



qlnorm(c(0.0005, 0.9995), 0, data_list$sigma_total_inf)

library(cmdstanr)
model <- cmdstan_model("src/gp_renewal_norm.stan")

fit <- model$optimize(
  data = data_list, seed = 5, iter = 1e4, init = 0.1,
  jacobian = TRUE
)
fit <- model$variational(
  data = data_list, seed = 3, iter = 1e4, init = 0.1
)

fit <- model$pathfinder(
  data = data_list, seed = 4, init = 0.1
)


fit <- model$sample(
  data = data_list,
  seed = 1,
  chains = 4,
  parallel_chains = 4,
  init = 0.01,
  refresh = 50,
  iter_warmup = 250, iter_sampling = 250
)

draws <- tidy_draws(fit) |> mutate(chain = .chain)

mcmc_hist(draws, c("log10_total_mean", "log10_total_inf_year[1]", "log10_total_inf_year[2]",
                   "log10_total_inf_year[3]", "log10_total_inf_year[4]"))

trajectory_data <- fit |> 
  spread_draws(log_R[ix_t], inf[ix_t], inf_delayed[ix_t],
               cases[ix_t], f_noise[ix_t],
               latent_cases_log[ix_t]
  ) |> 
  mutate(date = training_counts$date[ix_t],
         latent_cases = exp(latent_cases_log),
         cases = replace_na(cases, 0))


p_Rt <- ggplot(trajectory_data) + 
  fastribbon::stat_fastribbon(aes(x = date, y = exp(log_R))) +
  
  scale_fill_brewer() + theme_bw() + theme(legend.position = "none")


p_f_noise <- ggplot(trajectory_data) + 
  fastribbon::stat_fastribbon(aes(x = date, y = f_noise)) +
  
  scale_fill_brewer() + theme_bw() + theme(legend.position = "none")

p_cases <- ggplot(trajectory_data) + 
  fastribbon::stat_fastribbon(aes(x = date, y = inf_delayed))  + 
  fastribbon::stat_fastribbon(aes(x = date, y = latent_cases),
                              alpha = 0.2) +
  
  geom_line(aes(x = notification_date, y = value_trend), eval_counts,
            colour = "orangered")  +
  
  coord_cartesian(ylim = c(0, 400)) +
  
  scale_fill_brewer() + theme_bw() + theme(legend.position = "none")

p_Rt / p_cases


trajectory_data |> 
  filter((.draw - 1) %% 100 == 0) |>
  ggplot() +
  geom_line(aes(x = date, y = latent_cases, group = .draw)) +
  geom_line(aes(x = date, y = y_trend), training_counts,
            colour = "blue") +
  
  geom_line(aes(x = notification_date, y = value_trend), eval_counts,
            colour = "orangered") +
  
  theme_bw() + theme(legend.position = "none")


library(pammtools)
ggplot(trajectory_data) + 
  fastribbon::stat_fastribbon(aes(x = date, y = cases),
                              geom = "stepribbon") +
  
  geom_line(aes(x = notification_date, y = value), eval_counts,
            colour = "orange") +
  geom_line(aes(x = date, y = y), training_counts,
            colour = "black") +
  
  coord_cartesian(ylim = c(0, 400), xlim = c(origin_date - days(28), origin_date + days(7 * 8))) +
  
  scale_fill_brewer() + theme_bw() + theme(legend.position = "none")


