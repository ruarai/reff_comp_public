

source("R/functions.R")

round_id <- max(list.dirs("results", full.names = FALSE))

forecast_trajectories <- fs::dir_ls(fs::path("results", round_id), glob = "*forecast-trajectories.parquet")

forecasts <- open_dataset(forecast_trajectories)

plot_theme <- list(
  xlab(NULL), ylab(NULL),
  theme_minimal(),
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "none",
    panel.background = element_rect(colour = "grey70", linewidth = 0.25)
  )
)


plot_model_ribbon <- function(forecasts, combined_counts) {
  model_submissions <- forecasts |> 
    filter(target == "case_incidence", horizon <= 28) |> 
    collect()
  
  
  origin_dates <- forecasts |> 
    distinct(origin_date, pathogen, location) |> 
    collect()
  
  observations <- combined_counts |> 
    left_join(origin_dates, by = join_by(location, pathogen)) |> 
    filter(notification_date > origin_date - 7 * 8)
  
  library(pammtools)
  library(ggh4x)
  library(ufs)
  
  ggplot() +
    
    fastribbon::stat_fastribbon(aes(x = origin_date + horizon, y = value),
                                model_submissions,
                                geom = "stepribbon",
                                linewidth = 0.0) +
    
    geom_line(aes(x = notification_date, y = value),
              observations) +
    
    facet_wrap(~location * pathogen, ncol = 3, scales = "free_y") +
    
    ggh4x::facet_grid2(rows = vars(location), cols = vars(pathogen),
                       scales = "free_y", independent = "y") +

    geom_vline(aes(xintercept = origin_date),
               origin_dates, linetype = "44") +
    
    scale_y_continuous(breaks = scales::breaks_extended(3)) +
    
    coord_cartesian(ylim = c(0, NA)) +
    
    plot_theme +
    
    scale_fill_brewer()
}

year <- as.character(year(round_id))

combined_counts <- read_csv(str_c("../forecast-hub/processed-data/", year, "/", round_id, "/combined-count-", round_id, ".csv"))

plot_model_ribbon(forecasts |> filter(location != "NZ"), combined_counts |> filter(location != "NZ"))
ribbon_plot_path <- fs::path("plots", str_c("uom-reff_comp_plots-ribbon-", round_id, ".png"))

ggsave(
  ribbon_plot_path,
  height = 12,
  width = 10
)

library(magick)
forecast_plot_files <- fs::dir_ls(fs::path("results", round_id), glob = "*plot_forecast.png")
latent_plot_files <- fs::dir_ls(fs::path("results", round_id), glob = "*plot_latent_inf.png")
diagnostic_plot_files <- fs::dir_ls(fs::path("results", round_id), glob = "*plot_MCMC.png")

image_write(
  c(ribbon_plot_path, forecast_plot_files, latent_plot_files) |> 
    image_read() |> 
    image_resize("50%"),
  path = fs::path("plots", str_c("uom-reff_comp_plots-", round_id, ".pdf")), format = "pdf"
)

image_write(
  c(diagnostic_plot_files) |> 
    image_read() |> 
    image_resize("50%"),
  path = fs::path("plots", str_c("uom-reff_comp_plots-diagnostic-", round_id, ".pdf")), format = "pdf"
)


