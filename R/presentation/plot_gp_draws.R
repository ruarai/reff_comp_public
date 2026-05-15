
source("R/functions.R")
source("R/presentation/plot_theme.R")

require(cmdstanr)
file.remove("src/gp_renewal")
model <- cmdstan_model("src/gp_renewal.stan")


model$expose_functions(global = TRUE)



simulate_matern <- function(x, alpha_m, rho_m, seed = 0, ls_M = rho_m) {
  x <- x - mean(x)
  N <- length(x)
  
  M <- calculate_M_matern(c = 1.25, length_scale = ls_M, x = x)
  L <- max(abs(x)) * 1.25
  
  spd_matern <- diagSPD_Matern52(alpha = alpha_m, rho = rho_m, L = L, M = M)
  
  phi_matern <- PHI(N = N, M = M, L = L, x = x)
  
  set.seed(seed)
  beta_matern <- rnorm(n = M, mean = 0, sd = sqrt(spd_matern))
  
  phi_matern %*% beta_matern
}

simulate_periodic <- function(x, alpha_p, rho_p, seed = 0) {
  x <- x - mean(x)
  N <- length(x)
  M <- 40
  period  <- 365.25
  w0      <- 2 * pi / period
  
  rho_p_rad <- rho_p * (2 * pi / 365.25)
  
  spd_periodic <- diagSPD_Periodic(alpha = alpha_p, rho = rho_p_rad, M = M)
  
  phi_periodic <- PHI_periodic(N = N, M = M, w0 = w0, x = x)
  
  set.seed(seed)
  beta_periodic <- rnorm(n = 2 * M, mean = 0, sd = sqrt(spd_periodic))
  
  phi_periodic %*% beta_periodic
}

tibble(x = 1:128) |> mutate(y = simulate_matern(x, 1.0, 7, ls_M = 7)) |> 
  ggplot() +
  geom_line(aes(x,y))

n_years <- 3



p_ex_matern <- tibble(t = 0:150) |>
  mutate(y_noise = simulate_matern(t, 0.2, 14, seed = 3)) |> 
  
  ggplot() +
  geom_line(aes(x = t, y = y_noise),
            linewidth = 0.75) +
  
  # scale_x_continuous(breaks = seq(0, 365 * n_years, by = 365),
  #                    labels = ~ . / 365) +
  
  scale_x_continuous(breaks = scales::breaks_extended())  +
  
  scale_y_continuous(breaks = c(-0.5, 0, 0.5)) +
  
  ylab("f") +
  
  xlab("Time (days)") +
  
  coord_cartesian(ylim = c(-0.5, 0.5)) +
  
  ggtitle("Example trajectory") +
  
  plot_theme_pres

p_ex_matern



p_ex_periodic <- tibble(t = 0:(365 * n_years)) |>
  mutate(y_periodic = simulate_periodic(t, 0.2, 28, seed = 0)) |> 
  
  ggplot() +
  geom_line(aes(x = t, y = y_periodic),
            linewidth = 0.75) +
  
  scale_x_continuous(breaks = seq(0, 365 * n_years, by = 365),
                     labels = ~ . / 365) +
  
  scale_y_continuous(breaks = c(-1, 0, 1)) +
  
  ylab("f") +
  
  xlab("Time (years)") +
  
  coord_cartesian(ylim = c(-1.75, 1.75)) +
  
  ggtitle("Example trajectory") +
  
  plot_theme_pres


p_covar_matern <- tibble(x = seq(0, 150, by = 0.1), l = 14, alpha = 1.0) |> 
  mutate(covar = alpha * (1 + (sqrt(5) * x) / l + (5 * x^2) / (3 * l^2)) * exp(-sqrt(5) * x / l)) |> 
  
  ggplot() +
  geom_line(aes(x = x, y = covar),
            linewidth = 0.75) +
  
  ylab("Covariance") + xlab("Time difference (days)") +
  
  scale_x_continuous(breaks = scales::breaks_extended())  +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  
  ggtitle("Covariance") +
  
  plot_theme_pres



p_covar_periodic <- tibble(x = seq(0, 365 * n_years, by = 0.1), omega = 2 * pi / 365.25, alpha = 1.0, l = 20 * (2 * pi / 365.25)) |> 
  mutate(covar = alpha * exp(- 2 * (sin(omega * x / 2)^2) / l)) |> 
  
  ggplot() +
  geom_line(aes(x = x, y = covar),
            linewidth = 0.75) +
  
  ylab("Covariance") + xlab("Time difference (years)") +
  
  scale_x_continuous(breaks = seq(0, 365 * n_years, by = 365),
                     labels = ~ . / 365)  +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  
  ggtitle("Covariance") +
  
  plot_theme_pres


p_covar_matern / p_ex_matern

p_covar_periodic / p_ex_periodic


plot_lengthscale_periodic <- function(l) {
  p <- tibble(t = 0:(365 * n_years)) |>
    mutate(y_periodic = simulate_periodic(t, 0.2, l, seed = 0)) |> 
    
    ggplot() +
    geom_line(aes(x = t, y = y_periodic),
              linewidth = 0.75) +
    
    scale_x_continuous(breaks = seq(0, 365 * n_years, by = 365),
                       labels = ~ . / 365) +
    
    scale_y_continuous(breaks = c(-1, 0, 1)) +
    
    ylab("f") +
    
    xlab("Time (years)") +
    
    coord_cartesian(ylim = c(-1.75, 1.75)) +
    
    plot_theme_pres +
    
    ggtitle(str_c("Length-scale = ", floor(l)))
  
  return(p)
}



plots <- map(
  seq(5, 60, by = 0.5),
  plot_lengthscale_periodic
)


plots_to_gif(plots, "test.gif", 10, 5, delay = 1 / 12)

plot_lengthscale_periodic(60)

ggsave("test.png", width = 10, height = 5)



x <- seq(0, 50, by = 1)

rho_m <- 10
alpha_m <- 0.2

x <- x - mean(x)
N <- length(x)

M <- calculate_M_matern(c = 1.25, length_scale = rho_m, x = x)
L <- max(abs(x)) * 1.25

spd_matern <- diagSPD_Matern52(alpha = alpha_m, rho = rho_m, L = L, M = M)

phi_matern <- PHI(N = N, M = M, L = L, x = x)

phi_matern |> 
  reshape2::melt(c("i", "y")) |> 
  mutate(t = x[i]) |> 
  as_tibble() |> 
  
  ggplot() +
  geom_line(aes(x = t, y = value, group = y)) +
  
  ylab(NULL) + xlab("x") +
  
  facet_wrap(~y, ncol = 4) +
  
  plot_theme_pres



tibble(
  rho_m = c(5, 10, 15, 25)
) |> 
  rowwise() |> 
  mutate(spd = list(diagSPD_Matern52(alpha = alpha_m, rho = rho_m, L = L, M = M))) |> 
  unnest(spd) |> 
  group_by(rho_m) |> 
  mutate(x = 1:n()) |> 
  ggplot() +
  geom_line(aes(x = x, y = spd, group = rho_m, colour = factor(rho_m))) +
  
  xlab("Basis function index") +
  
  ylab("Variance") +
  
  scale_colour_discrete(name = "Length-scale") +
  
  plot_theme_pres




