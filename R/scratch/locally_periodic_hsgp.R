
source("R/functions.R")

require(cmdstanr)
file.remove("src/gp_renewal")
model <- cmdstan_model("src/gp_renewal.stan")

model$expose_functions(global = TRUE)

simulate_lp <- function(x, alpha, rho_se, rho_per, period = 365.25, seed = 0, M = 100, c = 1.25) {
  x <- x - mean(x)
  N <- length(x)
  
  L <- max(abs(x)) * c
  
  phi <- PHI(N = N, M = M, L = L, x = x)
  
  spd_lp <- diagSPD_LP(
    alpha = 1.5, rho_se = rho_se, rho_per = rho_per * (2 * pi / period),
    period = period, L = L, M = M, log_bessel_cutoff = -10
  )
  
  set.seed(seed)
  beta_lp <- rnorm(n = M, mean = 0, sd = spd_lp)
  
  return(c(phi %*% beta_lp))
}


true_kernel <- function(tau, alpha, rho_se, rho_per, period = 365.25) {
  
  alpha * exp(
    -tau^2 / (2 * rho_se ^ 2) -
      
      2 * (sin(pi * tau / period)^2) / (rho_per^2)
  )
}

x <- seq(0, 365 * 3, by = 1)
x <- x - mean(x)
N <- length(x)

M <- 800
c <- 3
L <- max(abs(x)) * c

period <- 7

phi <- PHI(N = N, M = M, L = L, x = x)

alpha <- 1.0
rho_se <- period * 2
rho_per <- 2


spd_lp <- diagSPD_LP(
  alpha = alpha, rho_se = rho_se, rho_per = rho_per * (2 * pi / period),
  period = period, L = L, M = M, log_bessel_cutoff = -10
) |> c()

phi_scaled <- sweep(phi, 2, spd_lp^2, "*")
k_matrix <- tcrossprod(phi_scaled, phi)
plot_data <- expand_grid(
  x_1 = x,
  x_2 = x
) |> 
  mutate(k = as.vector(t(k_matrix)))

plot_data_true <- expand_grid(
  x_1 = x,
  x_2 = x
) |> 
  mutate(tau = abs(x_1 - x_2),
         k = true_kernel(tau, alpha ^ 2, rho_se, rho_per * (2 * pi / period), period)) |> 
  
  left_join(plot_data |> rename(k_est = k))

max(abs(plot_data_true$k - plot_data_true$k_est))


ggplot(plot_data_true |> filter(x_1 %in% x_1[seq(1, length(x_1), length.out = 6)])) +
  geom_line(aes(x = x_2, y = k, colour = "True")) +
  geom_line(aes(x = x_2, y = k_est, colour = "Approximate")) +
  
  facet_wrap(~x_1)

ggplot(plot_data_true) +
  geom_raster(aes(x_1, x_2, fill = k_est))


ggplot(plot_data_true) +
  geom_raster(aes(x_1, x_2, fill = k - k_est))



# x <- seq(0, 365, by = 1)
sim_data <- tibble(x = x + 365 / 2, y = simulate_lp(x, alpha, rho_se, rho_per, period, M = M, c = c, seed = 1))

ggplot(sim_data) +
  geom_line(aes(x = x, y = y))


ggplot(sim_data) +
  geom_line(aes(x = x %% 365, y = y, colour = factor(floor(x / 365))))


get_error <- function(x, rho_se, rho_per, M, c) {
  x <- x - mean(x)
  N <- length(x)
  
  L <- max(abs(x)) * c
  
  period <- 365
  
  phi <- PHI(N = N, M = M, L = L, x = x)
  
  alpha <- 1.0
  
  spd_lp <- diagSPD_LP(
    alpha = alpha, rho_se = rho_se, rho_per = rho_per * (2 * pi / period),
    period = period, L = L, M = M, log_bessel_cutoff = -10
  ) |> c()
  
  phi_scaled <- sweep(phi, 2, spd_lp^2, "*")
  k_matrix <- tcrossprod(phi_scaled, phi)
  
  plot_data_true <- expand_grid(
    x_1 = x,
    x_2 = x
  ) |> 
    mutate(tau = abs(x_1 - x_2),
           k_est = as.vector(t(k_matrix)),
           k = true_kernel(tau, alpha ^ 2, rho_se, rho_per * (2 * pi / period), period))
  
  max(abs(plot_data_true$k - plot_data_true$k_est))
}


# Determine what range of rho_per, rho_se is safe for our
# range of x
x <- seq(0, 365 * 5, by = 1)

plot_data <- expand_grid(
  rho_per = seq(2, 10, by = 1) * 7,
  rho_se = seq(1, 5, by = 1) * 365
) |> 
  rowwise() |> 
  mutate(e = get_error(x, rho_se, rho_per, 600, 5))


ggplot(plot_data) +
  geom_line(aes(x = rho_per, y = e, colour = factor(rho_se / 365))) +
  
  scale_y_log10() +
  
  theme_bw()


