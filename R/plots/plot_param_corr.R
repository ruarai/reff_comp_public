
library(tidyverse)
library(cmdstanr)

fit <- read_rds("results/2026-04-29/RSV_NSW_case_incidence.rds")





parameter_names <- c(
  "scaling_annual", "scaling_noise", "scaling_weekly_asc", "ls_annual", "ls_noise",
  "log_R_0", "phi", "beta_annual", "beta_noise", "beta_offset_weekly", "r_0", "log_I_0",
  "school_holiday_Rt_effect_raw", "school_holiday_asc_effect_raw",
  "public_holiday_Rt_effect_raw", "public_holiday_asc_effect_raw",
  "delay_mean", "delay_sd"
)

draws <- fit$draws()

draws_sub <- posterior::subset_draws(draws, variable = str_c("^", parameter_names), regex = TRUE)

draws_chain <- posterior::as_draws_matrix(draws_sub)
draws_cor <- cor(draws_chain)

image(draws_cor)

high_cor_var <- draws_cor |> 
  as_tibble(rownames = "var_x") |> 
  pivot_longer(-var_x, names_to = "var_y") |> 
  filter(!(var_x %in% c("ls_annual_rad"))) |> 
  filter(var_x > var_y) |> 
  
  arrange(desc(abs(value))) |> 
  slice(1:40)

plots <- high_cor_var |> 
  pmap(function(var_x, var_y, value) {
    ggplot(tibble(x = as.numeric(draws_chain[, var_x]), y = as.numeric(draws_chain[, var_y]) )) +
      geom_point(aes(x, y),
                 size = 0.2, stroke = 0.2, alpha = 0.3) +
      ggtitle(str_c(var_x, "\n", var_y)) +
      
      xlab(NULL) + ylab(NULL) +
      
      theme_bw() +
      theme(plot.title = element_text(size = 8), plot.subtitle = element_text(size = 7),
            axis.text = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.margin = margin(0.1,0.1,0.1,0.1, unit = "cm"))
  })

plots[[1]]

wrap_plots(plots)




ggplot(tibble(
  x = as.numeric(draws_chain[, "scaling_weekly_asc"]),
  y = as.numeric(draws_chain[, "obs_offset_weekly[1]"]) )) +
  geom_point(aes(x, y),
             size = 0.2, stroke = 0.2, alpha = 0.3) +
  
  xlab(NULL) + ylab(NULL) +
  
  theme_bw() +
  theme(plot.title = element_text(size = 8), plot.subtitle = element_text(size = 7),
        axis.text = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(0.1,0.1,0.1,0.1, unit = "cm"))


