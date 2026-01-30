
library(tidyverse)
library(arrow)

source("R/functions.R")

forecasts <- open_dataset("forecast_outputs_raw/")

alpha_vals <- scales::rescale(rev(1/1.7^(1:8)), to = c(0.05, 0.99))


location_pathogens <- forecasts |> distinct(location, pathogen) |> collect()


forecast_quants <- forecasts |> 
  filter(target == "case_incidence_smoothed", horizon > 0, horizon <= 28) |> 
  collect() |> 
  mutate(date = origin_date + horizon) |> 
  make_intervals(c("value"),
                 group_vars = c("round_id", "location", "pathogen"),
                 sample_var = rlang::sym("output_type_id"),
                 na.rm = TRUE)

eval_counts <- read_csv("data/combined-count-smoothed-2026-01-15.csv") |> 
  filter(notification_date >= ymd("2025-05-01"),
         notification_date <= ymd("2025-12-01"))

plots <- location_pathogens |> 
  arrange(pathogen, location) |> 
  rename(x_location = location, x_pathogen = pathogen) |> 
  pmap(
    function(x_location, x_pathogen) {
      eval_counts_plot <- eval_counts |> 
        filter(location == x_location, pathogen == x_pathogen)
      
      plot_data <- forecast_quants |> 
        filter(location == x_location, pathogen == x_pathogen) |> 
        arrange(round_id) |> 
        mutate(row = match(round_id, unique(round_id)) %% 4)
      
      ggplot(plot_data) +
        geom_ribbon(aes(x = date, ymin = lower, ymax = upper,
                        group = interaction(quant, round_id), alpha = quant),
                    fill = ggokabeito::palette_okabe_ito(5)) +
        
        
        geom_line(aes(x = notification_date, y = cases_trend_and_remainder),
                  alpha = 0.8, colour = "black",
                  eval_counts_plot) +
        
        
        geom_vline(aes(xintercept = date),
                   plot_data |> group_by(round_id) |> filter(date == min(date)),
                   linetype = "44") +
        
        
        xlab("Notification date") + ylab("Count") +
        
        facet_wrap(~row, ncol = 1) +
        
        scale_alpha_manual(values = alpha_vals) +
        
        # scale_y_log10() +
        scale_x_date(breaks = scales::date_breaks(),
                     labels = scales::label_date_short()) +
        
        theme_bw() +
        
        theme(strip.text = element_blank(),
              legend.position = "none") +
        
        ggtitle(str_c(x_location, " - ", x_pathogen))
    }
  )

plots[[1]]

pdf("plots/collected_forecasts_smoothed_2.pdf", width = 10, height = 6)

for(p in plots) {
  plot(p)
}

dev.off()


