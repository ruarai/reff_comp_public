
source("R/functions.R")
source("R/presentation/plot_theme.R")

fit_files <- fs::dir_ls("saved_results/results_baseline/2025-10-09/", recurse = TRUE, glob = "*.rds")

holiday_draws <- map(
  fit_files,
  function(fit_file) {
    
    fit <- read_rds(fit_file)
    
    round_id <- ymd(str_extract(fit_file, "\\d{4}-\\d{2}-\\d{2}"))
    pathogen_location <- fs::path_file(fit_file) |> fs::path_ext_remove()
    pathogen <- str_split_fixed(pathogen_location, "_", 2)[1]
    location <- str_split_fixed(pathogen_location, "_", 2)[2]
    
    fit |> 
      
      spread_draws(school_holiday_Rt_effect_raw[ix_s],
                   school_holiday_asc_effect_raw[ix_s]) |> 
      mutate(location = !!location, pathogen = !!pathogen)
  },
  .progress = TRUE
) |> 
  bind_rows()

plot_data <- holiday_draws |> 
  mutate(school_holiday_name = c("Summer holidays", "Between-term holidays")[ix_s]) |> 
  filter(!(.chain == 3 & location == "NSW" & pathogen == "SARSCOV2"),
         location %in% c("NSW", "VIC", "WA", "ACT"))

p_Rt <- ggplot(plot_data) +
  
  ggdist::stat_slab(aes(x = exp(school_holiday_Rt_effect_raw), y = location),
                    colour = "black", trim = FALSE, expand = c(0.5, 1.5),
                    fill = ggokabeito::palette_okabe_ito(5),
                    density = "histogram",
                    scale = 1.5) +

  
  geom_vline(xintercept = 1, linetype = "44") +
  
  scale_x_log10(breaks = c(0.8, 1, 1.2)) +
  coord_cartesian(xlim = exp(c(-0.3, 0.3))) +
  
  plot_theme_pres_2 +
  
  ggokabeito::scale_fill_okabe_ito() +
  
  facet_grid(rows = vars(school_holiday_name), cols = vars(pathogen)) +
  
  ggtitle("Multiplicative effect of school holidays on R(t)") + 
  
  theme(legend.position = "none") +
  
  xlab(NULL) + ylab(NULL)


p_asc <- ggplot(plot_data) +
  
  ggdist::stat_slab(aes(x = exp(school_holiday_asc_effect_raw), y = location),
                    colour = "black", trim = FALSE, expand = c(0.5, 1.5),
                    fill = ggokabeito::palette_okabe_ito(3),
                    density = "histogram",
                    scale = 1.5) +
  
  
  geom_vline(xintercept = 1, linetype = "44") +
  
  scale_x_log10(breaks = c(0.75, 1, 1.3)) +
  coord_cartesian(xlim = exp(c(-0.5, 0.5))) +
  
  plot_theme_pres_2 +
  
  ggokabeito::scale_fill_okabe_ito() +
  
  facet_grid(rows = vars(school_holiday_name), cols = vars(pathogen)) +
  
  ggtitle("Multiplicative effect of school holidays on case asc.") + 
  
  theme(legend.position = "none") +
  
  xlab(NULL) + ylab(NULL)


p_Rt | p_asc


public_holiday_draws <- map(
  fit_files,
  function(fit_file) {
    
    fit <- read_rds(fit_file)
    
    round_id <- ymd(str_extract(fit_file, "\\d{4}-\\d{2}-\\d{2}"))
    pathogen_location <- fs::path_file(fit_file) |> fs::path_ext_remove()
    pathogen <- str_split_fixed(pathogen_location, "_", 2)[1]
    location <- str_split_fixed(pathogen_location, "_", 2)[2]
    
    public_holidays <- get_public_holidays() |> 
      filter(location == !!location) |> 
      pull(holiday_name) |> 
      unique()
    
    fit |> 
      
      spread_draws(public_holiday_Rt_effect_raw[ix_p],
                   public_holiday_asc_effect_raw[ix_p]) |> 
      mutate(location = !!location, pathogen = !!pathogen,
             holiday_name = public_holidays[ix_p])
  },
  .progress = TRUE
) |> 
  bind_rows()


plot_data <- public_holiday_draws |> 
  filter(!(.chain == 3 & location == "NSW" & pathogen == "SARSCOV2"),
         location == "NSW") |> 
  
  mutate(holiday_name = str_wrap(holiday_name, width = 20))


p_Rt <- ggplot(plot_data) +
  
  ggdist::stat_slab(aes(x = exp(public_holiday_Rt_effect_raw), y = holiday_name),
                    colour = "black", trim = FALSE, expand = c(0.5, 1.5),
                    fill = ggokabeito::palette_okabe_ito(5),
                    density = "histogram",
                    scale = 1.5) +
  
  
  geom_vline(xintercept = 1, linetype = "44") +
  
  scale_x_log10(breaks = c(0.5, 1, 2, 3)) +
  # coord_cartesian(xlim = exp(c(-0.3, 0.3))) +
  
  plot_theme_pres_2 +
  
  ggokabeito::scale_fill_okabe_ito() +
  
  facet_grid(cols = vars(pathogen)) +
  
  ggtitle("Multiplicative effect of public holidays on R(t)") + 
  
  theme(legend.position = "none") +
  
  xlab(NULL) + ylab(NULL)

p_asc <- ggplot(plot_data) +
  
  ggdist::stat_slab(aes(x = exp(public_holiday_asc_effect_raw), y = holiday_name),
                    colour = "black", trim = FALSE, expand = c(0.5, 1.5),
                    fill = ggokabeito::palette_okabe_ito(3),
                    density = "histogram",
                    scale = 1.5) +
  
  
  geom_vline(xintercept = 1, linetype = "44") +
  
  scale_x_log10(breaks = c(0.5, 1, 2, 3)) +
  # coord_cartesian(xlim = exp(c(-0.3, 0.3))) +
  
  plot_theme_pres_2 +
  
  ggokabeito::scale_fill_okabe_ito() +
  
  facet_grid(cols = vars(pathogen)) +
  
  ggtitle("Multiplicative effect of public holidays on case asc.") + 
  
  theme(legend.position = "none") +
  
  xlab(NULL) + ylab(NULL)

p_Rt | p_asc




plot_data <- public_holiday_draws |> 
  filter(!(.chain == 3 & location == "NSW" & pathogen == "SARSCOV2")) |> 
  
  mutate(holiday_name_short = str_wrap(holiday_name, width = 20),
         location = factor(location) |> fct_rev())

plot_data |> 
  ggplot() +
  ggdist::stat_interval(aes(x = exp(public_holiday_Rt_effect_raw), y = holiday_name)) +
  
  geom_vline(xintercept = 1, linetype = "44") +
  
  facet_grid(rows = vars(location), cols = vars(pathogen),
             scales = "free_y") +
  
  scale_colour_brewer() +
  
  xlab(NULL) + ylab(NULL) +
  
  plot_theme_pres_2 +
  
  ggtitle("Multiplicative effect of public holidays on R(t)") + 
  
  theme(legend.position = "none")


plot_data |> 
  ggplot() +
  ggdist::stat_interval(aes(x = exp(public_holiday_asc_effect_raw), y = location)) +
  
  geom_vline(xintercept = 1, linetype = "44") +
  
  facet_grid(rows = vars(holiday_name), cols = vars(pathogen),
             scales = "free_y") +
  
  scale_colour_brewer(palette = "Greens") +
  
  xlab(NULL) + ylab(NULL) +
  
  plot_theme_pres_2 +
  
  ggtitle("Multiplicative effect of public holidays on case asc.") + 
  
  theme(legend.position = "none",
        strip.text.y.right = element_text(angle = 0))


plot_data |> 
  filter(location != "NZ") |> 
  group_by(holiday_name) |> 
  filter(length(unique(location)) > 4) |> 
  ggplot() +
  ggdist::stat_interval(aes(x = exp(public_holiday_Rt_effect_raw), y = location),
                        linewidth = 2.0) +
  
  geom_vline(xintercept = 1, linetype = "44") +
  
  facet_grid(rows = vars(holiday_name_short), cols = vars(pathogen),
             scales = "free_y") +
  
  scale_colour_brewer(palette = "Blues") +
  
  xlab(NULL) + ylab(NULL) +
  
  plot_theme_pres_2 +
  
  ggtitle("Multiplicative effect of public holidays on R(t)") + 
  
  theme(legend.position = "none",
        strip.text.y.right = element_text(angle = 0),
        axis.text.y = element_text(size = 7))



plot_data |> 
  filter(location != "NZ") |> 
  group_by(holiday_name) |> 
  filter(length(unique(location)) > 4) |> 
  ggplot() +
  ggdist::stat_interval(aes(x = exp(public_holiday_asc_effect_raw), y = location),
                        linewidth = 2.0) +
  
  geom_vline(xintercept = 1, linetype = "44") +
  
  facet_grid(rows = vars(holiday_name_short), cols = vars(pathogen),
             scales = "free_y") +
  
  scale_colour_brewer(palette = "Greens") +
  
  xlab(NULL) + ylab(NULL) +
  
  plot_theme_pres_2 +
  ggtitle("Multiplicative effect of public holidays on case asc.") + 
  
  theme(legend.position = "none",
        strip.text.y.right = element_text(angle = 0),
        axis.text.y = element_text(size = 7))


