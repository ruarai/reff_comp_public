
library(ggtext)

plot_theme_paper <- list(
  theme_minimal(),
  theme(panel.grid = element_blank(),
        axis.line = element_line(),
        axis.ticks = element_line(),
        text = element_text(family = "Helvetica", colour = "black"),
        line = element_line(linewidth = 0.7),
        axis.text.x = element_text(margin = margin(t = 0.3, unit = "cm"), colour = "black"),
        axis.text.y = element_text(margin = margin(r = 0.3, unit = "cm"), colour = "black"),
        plot.margin = margin(t = 0.2, l = 0.2, b = 0.2, r = 0.2, unit = "cm"),
        axis.ticks.length=unit(-0.1, "cm"),
        axis.title.y = element_markdown(),
        axis.title.x = element_markdown(),
        strip.text.x = element_text(size = 11, hjust = 0),
        plot.title = element_markdown(),
        plot.subtitle = element_markdown()
  )
)

plot_theme_annotation <- theme(plot.title = element_markdown(family = "Helvetica"),
                               plot.margin = margin(t = 0.2, l = 0.2, b = 0.2, r = 0.2, unit = "cm"))

element_gridline <- element_line(colour = "grey70", linetype = "28", linewidth = 0.3)
element_facet_background <- element_rect(fill = "white", colour = "grey50", linewidth = 0.5)

# colour_fc <- "#0f7ea3"
colour_fc <- "#9831bd"
colour_fit <- "#3182BD"

label_pathogen <- function(tbl) {
  tbl |> 
    mutate(pathogen_label = factor(
      pathogen,
      c("RSV", "flu", "SARSCOV2"),
      labels = c("RSV" = "RSV", "flu" = "Influenza", "SARSCOV2" = "SARS-CoV-2")
    ))
}

locations_paper <- c("NSW", "QLD", "VIC", "NT", "ACT")
# locations_paper <- c("NSW", "WA", "VIC", "NT", "ACT")

process_perf_data <- function(perf_data) {
  
  processed_perf <- perf_data |> 
    mutate(horizon = as.numeric(date - origin_date)) |> 
    
    # Drop some specific dates with known issues
    anti_join(
      tribble(
        ~location, ~pathogen, ~date,
        "WA", "flu", ymd("2025-08-16"),
        "WA", "flu", ymd("2025-08-23"),
        "WA", "flu", ymd("2025-08-30"),
        "WA", "RSV", ymd("2025-08-16"),
        "WA", "RSV", ymd("2025-08-23"),
        "WA", "RSV", ymd("2025-08-30")
      ),
      by = join_by(location, pathogen, date)
    ) |> 
    
    filter(date >= ymd("2025-01-01"), # Remove forecasts that actually were produced on 2024 data
           ymd(round_id) <= ymd("2025-11-27"),
           
           # Remove an early artificial VIC round for consistency with other locations
           !(round_id == "2025-04-10" & location == "VIC" & pathogen == "flu"),
           
           location %in% locations_paper)
  
  return(processed_perf)
}

read_perf <- function(output_dir, model_name = output_dir) {
  summ_perf_file <- str_c("saved_results/", output_dir, "_performance.parquet")
  
  if(fs::file_exists(summ_perf_file)) {
    perf_data <- read_parquet(summ_perf_file) |> 
      collect() |> 
      mutate(.model = model_name)
  } else {
    perf_data <- open_dataset(
      fs::dir_ls(str_c("saved_results/", output_dir), recurse = TRUE, glob = "*performance.parquet")
    ) |> 
      collect() |> 
      mutate(.model = model_name)
  }
  
  return(perf_data)
}

# Mean CRPS (log scale) per model, pathogen, jurisdiction, round and horizon week
# (0 = days 1-7, ..., 3 = days 22-28), written by summarise_performance.R. The figures use
# this instead of the per-day scores, which are not released.
read_perf_summary <- function(models) {
  read_parquet("saved_results/perf_summary.parquet") |>
    filter(.model %in% models)
}

# Row order of jurisdictions in the heatmaps
location_y <- tibble(location = c("NSW", "QLD", "ACT", "NT", "VIC")) |>
  mutate(y_offset = row_number())


# Plot data for figures that show individual forecasts (Figs 1, 3 and Supp. Fig 6) ----
#
# The public repository does not include the raw model fits in saved_results/. Instead,
# the posterior draws these figures need are saved to data/paper/plot_data/ the first time
# they are extracted from a raw fit, and are read from there afterwards. Delete a file in
# data/paper/plot_data/ to re-extract it.
#
# Only model output and dates are saved. Observed case counts always come from the
# anonymised final dataset via read_final_counts().

plot_data_dir <- "data/paper/plot_data"

load_fit_trajectories <- function(
    plot_name, model_version, round_id, x_pathogen, x_location,
    vars = c("log_R", "cases", "f_annual", "f_noise"),
    window = c(-7 * 8, 7 * 4), # days relative to the forecast origin date
    n_draws = 2000,
    x_target = "case_incidence"
) {
  job_name <- str_c(x_pathogen, "_", x_location, "_", x_target)
  plot_data_file <- fs::path(
    plot_data_dir,
    str_c(plot_name, "_", model_version, "_", round_id, "_", job_name, ".parquet")
  )
  
  if(fs::file_exists(plot_data_file)) {
    return(read_parquet(plot_data_file))
  }
  
  fit_dir <- str_c("saved_results/", model_version, "/", round_id, "/")
  fit <- read_rds(str_c(fit_dir, job_name, ".rds"))
  training_counts <- read_parquet(str_c(fit_dir, job_name, "_training_counts.parquet"))
  
  origin_date <- training_counts$origin_date[[1]] - days(1)
  
  vars_spec <- rlang::parse_expr(str_c("c(", str_c(vars, collapse = ", "), ")[ix_t]"))
  
  plot_data <- fit |> 
    spread_draws(!!vars_spec, log_R_0, ndraws = n_draws, seed = 1) |> 
    ungroup() |> 
    mutate(date = training_counts$date[ix_t], origin_date = origin_date) |> 
    filter(date > origin_date + days(window[1]), date <= origin_date + days(window[2])) |> 
    select(.draw, date, origin_date, all_of(vars), log_R_0)
  
  fs::dir_create(plot_data_dir)
  write_parquet(plot_data, plot_data_file)
  
  return(plot_data)
}

read_final_counts <- function(x_pathogen, x_location, x_target = "case_incidence") {
  read_csv("data/combined-count-smoothed_final_anonymised.csv", show_col_types = FALSE) |> 
    filter(location == x_location, pathogen == x_pathogen, target == x_target) |> 
    select(date = notification_date, y = value)
}
