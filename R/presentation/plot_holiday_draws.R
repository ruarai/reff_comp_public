
public_holiday_draws <- read_parquet("data/paper/public_holiday_draws.parquet")

# Make a list of pre-2026 holidays (so we don't plot any that haven't been hit yet)
public_holidays <- get_public_holidays() |> 
  filter(location == "VIC") |> 
  filter(date < ymd("2025-12-01")) |> 
  distinct(location, public_holiday_name = holiday_name)

plot_data_public <- public_holiday_draws |> 
  filter(location == "VIC") |> 
  
  inner_join(public_holidays) |> 
  
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
  
  scale_x_log10(breaks = c(0.5, 1, 2), labels = \(x) scales::label_percent()(x - 1)) +
  coord_cartesian(xlim = exp(c(-1.1, 1.1))) +
  
  theme(panel.background = element_facet_background,
        panel.grid.major = element_gridline,
        legend.position = "none") +
  
  xlab("Effect") + ylab(NULL) +
  
  ggtitle("<b>A</b>. Effect of public holidays on <i>R</i>(<i>t</i>)")


p_asc_public <- plot_data_public |>
  
  ggplot() +
  
  ggdist::stat_interval(aes(x = exp(public_holiday_asc_effect_raw), y = public_holiday_name),
                        linewidth = 2.0) +
  
  geom_vline(xintercept = 1, linetype = "44") +
  
  facet_grid(rows = vars(location), cols = vars(pathogen_label),
             scales = "free_y") +
  
  scale_colour_manual(values = RColorBrewer::brewer.pal(4, "Greens")[2:4]) +
  
  plot_theme_paper +
  
  scale_x_log10(breaks = c(0.5, 1, 2), labels = \(x) scales::label_percent()(x - 1)) +
  coord_cartesian(xlim = exp(c(-1.1, 1.1))) +
  
  theme(panel.background = element_facet_background,
        panel.grid.major = element_gridline,
        legend.position = "none") +
  
  xlab("Effect") + ylab(NULL) +
  
  ggtitle("<b>B</b>. Effect of public holidays on case reporting")


p_pub <- p_Rt_public / p_asc_public
