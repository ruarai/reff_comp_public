source("R/functions.R")
source("R/paper/common.R")


delay_draws <- read_parquet("data/paper/delay_draws.parquet")


plot_data <- delay_draws |> 
  filter(location %in% locations_paper) |> 
  label_pathogen()

plot_data_summ <- plot_data |> 
  filter(days == 0)


n_draws_prior <- 2000
prior_dist <- tibble(
  delay_mean = rnorm(n_draws_prior, 9, 2.0),
  delay_sd = rlnorm(n_draws_prior, 0, 0.5)
) |> 
  filter(delay_mean > 3,
         delay_sd > 1)


prior_summ <- prior_dist |> 
  summarise(delay_mean_lower = quantile(delay_mean, 0.025),
            delay_mean_upper = quantile(delay_mean, 1 - 0.025),
            delay_sd_lower = quantile(delay_sd, 0.025),
            delay_sd_upper = quantile(delay_sd, 1 - 0.025))

p_mean <- ggplot(plot_data_summ) +
  
  ggdist::stat_interval(aes(x = delay_mean, y = location),
                        linewidth = 4.0) +
  
  geom_vline(xintercept = prior_summ$delay_mean_lower,
             linetype = "31", alpha = 0.5) +
  geom_vline(xintercept = prior_summ$delay_mean_upper,
             linetype = "31", alpha = 0.5) +
  
  scale_colour_manual(values = RColorBrewer::brewer.pal(4, "Blues")[2:4], name = "CrI Width") +
  
  facet_grid(cols = vars(pathogen_label)) +
  
  xlab("Delay mean") + ylab(NULL) +
  
  ggtitle("<b>A</b>. Infection-to-reporting delay mean") +
  
  plot_theme_paper +
  
  coord_cartesian(xlim = c(0, NA)) +
  
  theme(panel.grid.major = element_gridline,
        legend.position = "none")

p_mean

p_sd <- ggplot(plot_data_summ) +
  
  ggdist::stat_interval(aes(x = delay_sd, y = location),
                        linewidth = 4.0) +
  
  geom_vline(xintercept = prior_summ$delay_sd_lower,
             linetype = "31", alpha = 0.5) +
  geom_vline(xintercept = prior_summ$delay_sd_upper,
             linetype = "31", alpha = 0.5) +
  
  scale_colour_manual(values = RColorBrewer::brewer.pal(4, "Blues")[2:4], name = "CrI Width") +
  
  facet_grid(cols = vars(pathogen_label)) +
  
  xlab("Delay s.d.") + ylab(NULL) +
  
  ggtitle("<b>B</b>. Infection-to-reporting delay s.d.") +
  
  plot_theme_paper +
  
  coord_cartesian(xlim = c(0, 10)) +
  
  theme(panel.grid.major = element_gridline,
        legend.position = "none")


plot_data_curve <- plot_data |> 
  mutate(location = fct_rev(factor(location))) |> 
  group_by(location, pathogen, .draw) |> 
  mutate(delay_cdf = cumsum(delay_dist))


prior_curve <- prior_dist |> 
  mutate(
    sigma = sqrt(log(1 + (delay_sd ^ 2) / (delay_mean ^ 2))),
    mu = log(delay_mean) - 0.5 * sigma ^ 2
  ) |> 
  
  expand_grid(
    days = 0:21
  ) |> 
  
  mutate(delay_cdf = plnorm(days, mu, sigma))


library(pammtools)

p_curve <- ggplot() +
  fastribbon::stat_fastribbon(aes(x = days, y = delay_cdf),
                              plot_data_curve,
                              geom = "stepribbon") +
  
  
  fastribbon::stat_fastribbon(aes(x = days, y = delay_cdf),
                              prior_curve, .width = 0.95,
                              fill = "transparent", colour = "#000000aa",
                              linetype = "31",
                              geom = "stepribbon") +
  
  
  scale_fill_manual(values = RColorBrewer::brewer.pal(4, "Blues")[2:4], name = "CrI Width") +
  
  xlab("Delay (days)") + ylab(NULL) +
  
  facet_grid(cols = vars(pathogen_label), rows = vars(location)) +
  
  ggtitle("<b>C</b>. Infection-to-reporting delay distribution CMF") +
  
  plot_theme_paper +
  
  theme(panel.grid.major.x = element_gridline,
        legend.position = "bottom")

p_curve


(p_mean / p_sd / p_curve) +
  plot_layout(heights = c(1, 1, 3))

ggsave(
  "paper_figures/supp/f_delay.png",
  width = 8, height = 10,
  bg = "white"
)


