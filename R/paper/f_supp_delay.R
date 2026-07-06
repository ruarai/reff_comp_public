source("R/functions.R")
source("R/paper/common.R")


delay_draws <- read_parquet("data/paper/delay_draws.parquet")


plot_data <- delay_draws |> 
  filter(location %in% locations_paper) |> 
  label_pathogen()

plot_data_summ <- plot_data |> 
  filter(days == 0)

p_mean <- ggplot(plot_data_summ) +
  
  ggdist::stat_interval(aes(x = delay_mean, y = location),
                        linewidth = 4.0) +
  
  scale_colour_manual(values = RColorBrewer::brewer.pal(4, "Blues")[2:4], name = "CrI Width") +
  
  facet_grid(cols = vars(pathogen_label)) +
  
  xlab("Delay mean") + ylab(NULL) +
  
  ggtitle("<b>A</b>. Infection-to-reporting delay mean") +
  
  plot_theme_paper +
  
  theme(panel.grid.major = element_gridline,
        legend.position = "none")


p_sd <- ggplot(plot_data_summ) +
  
  ggdist::stat_interval(aes(x = delay_sd, y = location),
                        linewidth = 4.0) +
  
  scale_colour_manual(values = RColorBrewer::brewer.pal(4, "Blues")[2:4], name = "CrI Width") +
  
  facet_grid(cols = vars(pathogen_label)) +
  
  xlab("Delay s.d.") + ylab(NULL) +
  
  ggtitle("<b>B</b>. Infection-to-reporting delay s.d.") +
  
  plot_theme_paper +
  
  theme(panel.grid.major = element_gridline,
        legend.position = "none")



plot_data_curve <- plot_data |> 
  group_by(location, pathogen, .draw) |> 
  mutate(delay_cdf = cumsum(delay_dist))

library(pammtools)

p_curve <- plot_data_curve |> 
  
  ggplot() +
  fastribbon::stat_fastribbon(aes(x = days, y = delay_cdf), geom = "stepribbon") +
  scale_fill_manual(values = RColorBrewer::brewer.pal(4, "Blues")[2:4], name = "CrI Width") +
  
  xlab("Delay (days)") + ylab(NULL) +
  
  facet_grid(cols = vars(pathogen_label), rows = vars(location)) +
  
  ggtitle("<b>C</b>. Infection-to-reporting delay distribution CMF") +
  
  plot_theme_paper +
  
  theme(panel.grid.major = element_gridline,
        legend.position = "bottom")

p_curve


(p_mean / p_sd / p_curve) +
  plot_layout(heights = c(1, 1, 3))

ggsave(
  "paper_figures/supp/f_delay.png",
  width = 8, height = 10,
  bg = "white"
)


