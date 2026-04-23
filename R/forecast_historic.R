
source("R/functions.R")

round_id <- "2026-04-15"
year <- as.character(year(round_id))

combined_counts <- read_csv(str_c("../forecast-hub/processed-data/", year, "/", round_id, "/combined-count-", round_id, ".csv"))

x_location <- "VIC"
x_pathogen <- "RSV"
x_target <- "case_incidence"

origin_date <- ymd("2020-01-01")

print(str_c("Performing job ", ix_job, " (", x_location, "-", x_pathogen, "-", x_target, ")"))

n_days_ahead <- 0

process_counts_historic <- function(
    combined_counts,
    x_location, x_pathogen, x_target
) {
  training_dates <- get_training_dates()
  
  n_days_ahead <- 0
  
  # Filter the counts for our given location and pathogen
  # and drop data prior to the training date
  reduced_counts <- combined_counts |> 
    filter(notification_date < ymd("2020-01-01")) |> 
    
    filter(location == x_location, pathogen == x_pathogen, target == x_target) |> 
    filter(is_complete) |> 
    
    filter(
      # Drop the last day of data (often issues)
      # If this improves, remove this line
      notification_date < max(notification_date)
    )
  
  training_counts <- reduced_counts |> 
    make_training_counts(x_location, x_pathogen, n_days_ahead, holidays_cutoff = ymd("2015-01-01"))
  
  return(training_counts)
}

training_counts <- process_counts_historic(combined_counts, x_location, x_pathogen, x_target) |> 
  # Save some extra data in case we need it later
  mutate(origin_date = origin_date)

data_list <- create_stan_data(training_counts, n_days_ahead = 0, x_pathogen, x_location)

fit <- fit_model(data_list, model_file = "src/gp_renewal.stan")

# fs::dir_create(fs::path("results", round_id))
# 
# fit_file <- str_c("results/", round_id, "/", x_pathogen, "_", x_location, "_", x_target, ".rds")
# fit$save_object(fit_file)




draws <- fit |>
  tidybayes::tidy_draws() |> 
  
  # NA cases values can cause issues with bayesplot, so
  # drop the columns early
  select(-starts_with("cases"))

pars <- c(
  "scaling_annual", "scaling_noise", "ls_noise", "ls_annual", "phi",
  "delay_mean", "delay_sd", "log_R_0"
)

p_hist <- mcmc_hist(draws |> mutate(chain = .chain), pars = pars) +
  theme_bw()

p_trace <- mcmc_trace(draws |> mutate(chain = .chain), pars = pars) +
  theme_bw()

p_hist | p_trace

p_pairs <- mcmc_pairs(draws |> mutate(chain = .chain), pars = pars, off_diag_args = list(size = 0.3, alpha = 0.5))

school_holidays <- get_school_holidays() |> 
  filter(location == x_location)

public_holidays <- get_public_holidays() |> 
  filter(location == x_location)

trajectory_data <- fit |> 
  spread_draws(log_R[ix_t], inf[ix_t], inf_delayed[ix_t],
               cases[ix_t], f_annual[ix_t], f_noise[ix_t],
               latent_cases_log[ix_t]) |> 
  mutate(date = training_counts$date[ix_t],
         latent_cases = exp(latent_cases_log),
         
         cases = replace_na(cases, 0))

p_reff <- plot_reff_components(trajectory_data, origin_date, school_holidays, public_holidays) +
  xlim(c(ymd(min(trajectory_data$date)), ymd(origin_date)))

p_latent_inf <- plot_latent_inf(trajectory_data, training_counts, origin_date, school_holidays, public_holidays) +
  
  xlim(c(ymd(min(trajectory_data$date)), ymd(origin_date))) +
  
  scale_y_log10()

p_reff / p_latent_inf






