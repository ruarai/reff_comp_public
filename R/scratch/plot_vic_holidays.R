source("R/functions.R")
source("R/paper/common.R")


holiday_draws <- read_parquet("data/paper/holiday_draws.parquet")
gp_draws <- read_parquet("data/paper/gp_draws.parquet")


plot_data <- holiday_draws |> 
  mutate(school_holiday_name = c("Summer holidays", "Between-term holidays")[ix_s]) |> 
  filter(location %in% "VIC") |> 
  
  label_pathogen()

plot_data_scaling <- holiday_draws |> 
  filter(location %in% "VIC") |> 
  group_by(location, pathogen) |> 
  summarise(scaling_noise = mean(scaling_noise),
            scaling_lower = mean(exp(scaling_noise)),
            scaling_upper = mean(exp(-scaling_noise))) |> 
  
  label_pathogen()


plot_effect <- function(plot_data, effect_name, holiday_name, pal = "Blues") {
  plot_data |> 
    filter(school_holiday_name == holiday_name) |>  
    rename(x = !!effect_name) |> 
    ggplot() +
    
    ggdist::stat_interval(aes(x = exp(x), y = pathogen_label),
                          linewidth = 4.0) +
    
    
    geom_vline(xintercept = 1, linetype = "44") +
    
    plot_theme_paper +
    
    # scale_colour_brewer(palette = pal) +
    
    scale_colour_manual(values = RColorBrewer::brewer.pal(4, pal)[2:4]) +
    
    # facet_grid(cols = vars(pathogen_label))+  
    
    theme(legend.position = "none",
          panel.grid.major = element_gridline,
          plot.margin = margin()) +
    
    xlab("Multiplicative effect") + ylab(NULL)
}

scale_Rt <- list(
  scale_x_log10(breaks = c(0.8, 1, 1.2)),
  coord_cartesian(xlim = exp(c(-0.3, 0.3)))
)

p_Rt_1 <- plot_effect(plot_data, "school_holiday_Rt_effect_raw", "Between-term holidays") +
  
  # geom_linerange(aes(y = location, xmin = scaling_lower, xmax = scaling_upper), plot_data_scaling,
  #            alpha = 0.8, linetype = "11") +
  scale_Rt

p_Rt_1

p_Rt <- (p_Rt_1) + 
  plot_annotation("Multiplicative effect of between-term school<br>holidays on <i>R</i>(<i>t</i>)",
                  theme = plot_theme_annotation)



p_asc_1 <- plot_effect(plot_data, "school_holiday_asc_effect_raw", "Between-term holidays", "Greens") +
  scale_Rt


p_asc <- (p_asc_1) + 
  plot_annotation("Multiplicative effect of between-term school<br>holidays on reporting",
                  theme = plot_theme_annotation)



p_school <- wrap_elements(full = p_Rt) /
  wrap_elements(full = p_asc)

p_school

ggsave(
  "est_effect_VIC.png",
  width = 5, height = 4.5,
  bg = "white"
)



ggsave(
  "paper_figures/f_5_estimated_effects.png",
  width = 10, height = 6.5,
  bg = "white"
)
