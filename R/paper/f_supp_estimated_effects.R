source("R/functions.R")
source("R/paper/common.R")


holiday_draws <- read_parquet("data/paper/holiday_draws.parquet")


plot_effect <- function(plot_data, effect_name, holiday_name, pal = "Blues") {
  plot_data |> 
    filter(school_holiday_name == holiday_name) |>  
    rename(x = !!effect_name) |> 
    ggplot() +
    
    ggdist::stat_interval(aes(x = exp(x), y = location),
                          linewidth = 4.0) +
    
    
    geom_vline(xintercept = 1, linetype = "44") +
    
    plot_theme_paper +
    
    # scale_colour_brewer(palette = pal) +
    
    scale_colour_manual(values = RColorBrewer::brewer.pal(4, pal)[2:4]) +
    
    facet_grid(cols = vars(pathogen_label))+  
    
    theme(legend.position = "none",
          panel.grid.major = element_gridline,
          plot.margin = margin()) +
    
    xlab("Multiplicative effect") + ylab(NULL)
}

scale_Rt <- list(
  scale_x_log10(breaks = c(0.6,  0.8, 1, 1.2, 1.67)),
  coord_cartesian(xlim = exp(c(-0.6, 0.6)))
)


plot_data <- holiday_draws |> 
  mutate(school_holiday_name = c("Summer holidays", "Between-term holidays")[ix_s]) |> 
  filter(location %in% locations_paper) |> 
  
  label_pathogen()

p_asc_2 <- plot_effect(plot_data, "school_holiday_asc_effect_raw", "Summer holidays", "Greens") +
  scale_Rt +
  theme(plot.margin = margin(t = 0.2, unit = "cm"))

p_Rt_2 <- plot_effect(plot_data, "school_holiday_Rt_effect_raw", "Summer holidays") +
  scale_Rt +
  theme(plot.margin = margin(t = 0.2, b = 0.2, unit = "cm"))

p_Rt <- (p_Rt_2) + 
  plot_annotation("<b>A</b>. Multiplicative effect of summer school holidays on <i>R</i>(<i>t</i>)",
                  theme = plot_theme_annotation)


p_asc <- (p_asc_2) + 
  plot_annotation("<b>B</b>. Multiplicative effect of summer school holidays on reporting",
                  theme = plot_theme_annotation)





p_school <- wrap_elements(full = p_Rt) /
  wrap_elements(full = p_asc)

p_school

ggsave(
  "paper_figures/f_supp_summer_holiday_estimates.png",
  width = 8, height = 5.5,
  bg = "white"
)



public_holiday_draws <- read_parquet("data/paper/public_holiday_draws.parquet")

plot_data_public <- public_holiday_draws |> 
  filter(location %in% locations_paper) |> 
  
  label_pathogen()

p_Rt_public <- plot_data_public |>
  
  ggplot() +
  
  ggdist::stat_interval(aes(x = exp(public_holiday_Rt_effect_raw), y = public_holiday_name),
                        linewidth = 2.0) +
  
  geom_vline(xintercept = 1, linetype = "44") +
  
  facet_grid(rows = vars(location), cols = vars(pathogen_label),
             scales = "free_y") +
  
  scale_colour_manual(values = RColorBrewer::brewer.pal(4, "Blues")[2:4]) +
  
  plot_theme_paper +
  
  scale_x_log10(breaks = c(0.5, 1, 2)) +
  coord_cartesian(xlim = exp(c(-1.1, 1.1))) +
  
  theme(panel.background = element_facet_background,
        panel.grid.major = element_gridline,
        legend.position = "none") +
  
  xlab("Multiplicative effect") + ylab(NULL) +
  
  ggtitle("<b>B</b>. Multiplicative effect of public holidays on <i>R</i>(<i>t</i>)")


p_asc_public <- plot_data_public |>
  
  ggplot() +
  
  ggdist::stat_interval(aes(x = exp(public_holiday_asc_effect_raw), y = public_holiday_name),
                        linewidth = 2.0) +
  
  geom_vline(xintercept = 1, linetype = "44") +
  
  facet_grid(rows = vars(location), cols = vars(pathogen_label),
             scales = "free_y") +
  
  scale_colour_manual(values = RColorBrewer::brewer.pal(4, "Greens")[2:4]) +
  
  plot_theme_paper +
  
  scale_x_log10(breaks = c(0.5, 1, 2)) +
  coord_cartesian(xlim = exp(c(-1.1, 1.1))) +
  
  theme(panel.background = element_facet_background,
        panel.grid.major = element_gridline,
        legend.position = "none") +
  
  xlab("Multiplicative effect") + ylab(NULL) +
  
  ggtitle("<b>B</b>. Multiplicative effect of public holidays on case reporting")


p_Rt_public / p_asc_public

ggsave(
  "paper_figures/f_supp_public_holiday_estimates.png",
  width = 10, height = 14,
  bg = "white"
)
