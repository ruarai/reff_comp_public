
source("R/functions.R")
source("R/paper/common.R")

require(cmdstanr)
file.remove("src/gp_renewal")
model <- cmdstan_model("src/gp_renewal.stan")

model$expose_functions(global = TRUE)

simulate_matern <- function(x, alpha_m, rho_m, seed = 0) {
  x <- x - mean(x)
  N <- length(x)
  
  M <- calculate_M_matern(c = 1.25, length_scale = 7, x = x)
  L <- max(abs(x)) * 1.25
  
  spd_matern <- diagSPD_Matern52(alpha = alpha_m, rho = rho_m, L = L, M = M)
  
  phi_matern <- PHI(N = N, M = M, L = L, x = x)
  
  set.seed(seed)
  beta_matern <- rnorm(n = M, mean = 0, sd = sqrt(spd_matern))
  
  c(phi_matern %*% beta_matern)
}

x <- seq(0, 365 * 3, by = 1)
alpha_m <- 1.0
rho_m <- 14.0


sim_ls <- map(
  1:10000,
  function(i) {
    simulate_matern(x, alpha_m, rho_m, seed = i)
  }
)

sim_mat <- do.call(rbind, sim_ls)
image(sim_mat)
cor_mat <- cor(sim_mat)

reshape2::melt(cor_mat, c("i", "j"), value.name = "cor") |> 
  as_tibble() |> 
  mutate(i = x[i], j = x[j]) |> 
  
  ggplot() +
  
  geom_tile(aes(i, j, fill = cor)) +
  
  scale_fill_viridis_c() +
  
  plot_theme_paper


fit <- read_rds("results/2026-04-29/flu_VIC_case_incidence.rds")

training_counts <- read_parquet("results/2026-04-29/flu_VIC_case_incidence_training_counts.parquet")

origin_date <- training_counts$origin_date[[1]]

plot_draws <- fit |> 
  spread_draws(f_noise[i]) |> 
  mutate(date = training_counts$date[i]) |> 
  filter(.draw %% 500 == 0)

plot_draws |> 
  ggplot() +
  geom_line(aes(x = date, y = f_noise, group = .draw, colour = factor(.draw)),
            alpha = 0.4) +
  
  geom_hline(yintercept = 0, linetype = "22") +
  geom_vline(xintercept = origin_date, linetype = "44") +
  
  plot_theme_paper +
  theme(legend.position = "none") +
  
  ggtitle(NULL, "f_noise trajectories (8 trajectories)")

draws_noise <- fit |> 
  spread_draws(f_noise[i]) |> 
  select(i, .draw, f_noise) |> 
  pivot_wider(names_from = .draw, values_from = f_noise,
              names_prefix = "draw_") |> 
  ungroup() |> 
  select(-i) |> 
  as.matrix() |> 
  t()

cor_noise <- cor(draws_noise)

image(cor_noise)

reshape2::melt(cor_noise, c("i", "j"), value.name = "corr") |> 
  as_tibble() |> 
  mutate(i = training_counts$date[i], j = training_counts$date[j]) |>
  
  ggplot() +
  
  geom_tile(aes(i, j, fill = corr)) +
  
  geom_hline(yintercept = origin_date, linetype = "44", alpha = 0.5,
             colour = "white") +
  geom_vline(xintercept = origin_date, linetype = "44", alpha = 0.5,
             colour = "white") +
  
  scale_fill_viridis_c() +
  
  coord_fixed() +
  xlab(NULL) + ylab(NULL) +
  
  ggtitle(NULL, "Correlation in f_noise (all trajectories)") +
  
  plot_theme_paper








