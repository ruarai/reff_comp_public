source("R/functions.R")
source("R/paper/common.R")


holiday_draws <- read_parquet("data/paper/holiday_draws.parquet")
gp_draws <- read_parquet("data/paper/gp_draws.parquet")


plot_data <- holiday_draws |> 
  mutate(school_holiday_name = c("Summer holidays", "Between-term holidays")[ix_s]) |> 
  filter(location %in% locations_paper) |> 
  
  label_pathogen()


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
  scale_x_log10(breaks = c(0.8, 1, 1.2)),
  coord_cartesian(xlim = exp(c(-0.3, 0.3)))
)

p_Rt_1 <- plot_effect(plot_data, "school_holiday_Rt_effect_raw", "Between-term holidays") +
  scale_Rt


p_Rt <- (p_Rt_1) + 
  plot_annotation("<b>A</b>. Multiplicative effect of between-term school<br>holidays on <i>R</i>(<i>t</i>)",
                  theme = plot_theme_annotation)



p_asc_1 <- plot_effect(plot_data, "school_holiday_asc_effect_raw", "Between-term holidays", "Greens") +
  scale_Rt


p_asc <- (p_asc_1) + 
  plot_annotation("<b>B</b>. Multiplicative effect of between-term school<br>holidays on reporting",
                  theme = plot_theme_annotation)

plot_data_gp <- gp_draws |> 
  filter(location %in% locations_paper,
         year(date) == 2025) |> 
  
  label_pathogen()

plot_data_scaling <- holiday_draws |> 
  filter(location %in% locations_paper) |> 
  group_by(location, pathogen) |> 
  summarise(scaling_noise = mean(scaling_noise),
            scaling_lower = mean(exp(scaling_noise)),
            scaling_upper = mean(exp(-scaling_noise))) |> 
  
  label_pathogen()


p_ann <- ggplot(plot_data_gp) +
  fastribbon::stat_fastribbon(aes(x = date, y = exp(f_annual)),
                              linewidth = 0.0) +
  
  geom_hline(aes(yintercept = scaling_lower), plot_data_scaling,
             alpha = 0.8, linetype = "11") +
  
  geom_hline(aes(yintercept = scaling_upper), plot_data_scaling,
             alpha = 0.8, linetype = "11") +
  
  geom_hline(yintercept = 1, linetype = "44") +
  
  facet_grid(rows = vars(fct_rev(location)), cols = vars(pathogen_label)) +
  
  scale_y_log10(breaks = c(0.8, 1, 1.25)) +
  
  scale_x_date(labels = scales::label_date_short(format = c("", "%b"))) +
  
  scale_fill_brewer() +
  
  plot_theme_paper +
  xlab(NULL) + ylab("Multiplicative effect") +
  theme(legend.position = "none",
        panel.grid.major.x = element_gridline,
        panel.background = element_facet_background,
        panel.spacing.x = unit(0.5, "cm")) + 
  plot_annotation("<b>C</b>. Multiplicative effect of periodic GP <i>f</i><sub>ann</sub>(<i>t</i>) on <i>R</i>(<i>t</i>)",
                  theme = theme(plot.title = element_markdown(family = "Helvetica"),
                                plot.margin = margin(t = 0.4, l = 0.2, b = 0.0, r = 0.2, unit = "cm")))

p_ann


p_school <- wrap_elements(full = p_Rt) /
  wrap_elements(full = p_asc)


wrap_elements(full = p_school) |
  wrap_elements(full = p_ann)


ggsave(
  "paper_figures/f_5_estimated_effects.png",
  width = 10, height = 6.5,
  bg = "white"
)

# 
# p_asc_2 <- plot_effect(plot_data, "school_holiday_asc_effect_raw", "Summer holidays", "Greens") +
#   ggtitle(NULL, "Summer holidays") + scale_asc +
#   theme(plot.margin = margin(t = 0.2, unit = "cm"))
# 
# p_Rt_2 <- plot_effect(plot_data, "school_holiday_Rt_effect_raw", "Summer holidays") +
#   ggtitle(NULL, "Summer holidays") + scale_Rt +
#   theme(plot.margin = margin(t = 0.2, b = 0.2, unit = "cm"))
# 
# 
# 
# 
# 
# ggplot(plot_data |> filter(ix_s == 1)) +
#   geom_point(aes(x = school_holiday_Rt_effect_raw, y = school_holiday_asc_effect_raw)) +
#   
#   facet_wrap(~location)

