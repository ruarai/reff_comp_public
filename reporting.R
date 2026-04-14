

forecast_trajectories <- fs::dir_ls(fs::path("results", round_id), glob = "*forecast-trajectories.parquet")

forecasts <- open_dataset(forecast_trajectories)

peak_forecasts <- forecasts |> 
  filter(target %in% c("peak_case_incidence", "peak_day_of_year")) |> 
  collect() |> 
  mutate(.model = "reff-comp") |> 
  pivot_wider(names_from = "target")

peak_date_range <- c(floor_date(ymd(round_id), "year"), floor_date(ymd(round_id), "year") + years(1))

plots <- peak_forecasts |> 
  distinct(x_location = location, x_pathogen = pathogen, x_model = .model) |> 
  arrange(x_pathogen, x_location) |> 
  collect() |> 
  pmap(
    function(x_location, x_pathogen, x_model) {
      
      training_counts <- read_parquet(fs::path("results", round_id, str_c(x_pathogen, "_", x_location, "_training_counts.parquet")))
      
      origin_date <- training_counts$origin_date[[1]]
      
      plot_data_training <- training_counts |> 
        filter(year(date) == year(round_id))
      
      peak_data <- peak_forecasts |> 
        filter(location == x_location,
               pathogen == x_pathogen,
               .model == x_model) |> 
        
        drop_na(peak_case_incidence) |> 
        mutate(date = floor_date(round_id, "year") + peak_day_of_year - 1)
      
      upper_cases <- quantile(peak_data$peak_case_incidence, 0.9)
      
      plot_data <- peak_data |> 
        
        filter(peak_case_incidence < upper_cases)
      
      p_scatter <- ggplot() +
        geom_density2d_filled(aes(x = date, y = peak_case_incidence),
                              n = 512, bins = 6,
                              plot_data) +
        
        geom_point(aes(x = date, y = peak_case_incidence),
                   plot_data |> slice_sample(n = 500),
                   size = 0.1, stroke = 0.5, alpha = 0.4) +
        
        geom_line(aes(x = date, y = y),
                  plot_data_training) +
        
        geom_vline(xintercept = origin_date, linetype = "44") +
        
        # scale_y_log10() +
        scale_x_date(date_breaks = "month") +
        
        xlab(NULL) + ylab("Case incidence") +
        
        coord_cartesian(xlim = peak_date_range) +
        
        # scale_fill_brewer(palette = "Blues") +
        scale_fill_manual(values = c("transparent", RColorBrewer::brewer.pal(5, "Blues"))) +
        
        theme_bw() +
        
        ggtitle(glue::glue("{x_location}-{x_pathogen}-{x_model}")) +
        
        theme(legend.position = "none",
              axis.text.x = element_blank())
      
      p_date <- ggplot(peak_data) +
        geom_histogram(aes(x = date, y = after_stat(ndensity)), binwidth = 7) +
        stat_ecdf(aes(x = date)) +
        
        scale_x_date(date_breaks = "month", labels = scales::label_date_short()) +
        
        geom_vline(xintercept = origin_date, linetype = "44") +
        
        xlab(NULL) + ylab(NULL) +
        
        coord_cartesian(xlim = peak_date_range) +
        
        theme_bw()
      
      p_cases <- ggplot() +
        geom_histogram(aes(y = peak_case_incidence, x = after_stat(ndensity)), 
                       breaks = seq(0, upper_cases, length.out = 30), plot_data) +
        stat_ecdf(aes(y = peak_case_incidence), peak_data) +
        
        coord_cartesian(ylim = c(0, upper_cases)) +
        
        xlab(NULL) + ylab(NULL) +
        
        theme_bw()
      
      p <- wrap_plots(
        p_scatter, p_cases,
        p_date, plot_spacer(),
        ncol = 2,
        widths = c(0.75, 0.25), heights = c(0.75, 0.25)
      )
      
      return(p)
    }
  )

save_plots_to_pdf(plots, fs::path("plots", str_c(round_id, "_peaks.pdf")), width = 10, height = 7)
