


source("R/functions.R")
source("R/paper/common.R")

require(cmdstanr)
file.remove("src/gp_renewal")
model <- cmdstan_model("src/gp_renewal.stan")

model$expose_functions(global = TRUE)

true_kernel_matern52 <- function(x, ls_noise) {
  kernel <- function(tau, ls_noise) {
    (1 + sqrt(5) * tau / ls_noise + (5 * tau ^ 2) / (3 * ls_noise ^ 2)) * exp(-sqrt(5) * tau / ls_noise)
  }
  
  tbl <- expand_grid(
    x_1 = x,
    x_2 = x
  ) |> 
    mutate(k_true = kernel(abs(x_1 - x_2), ls_noise))
  
  return(tbl)
}

approx_kernel_matern52 <- function(x, ls_noise, c, M) {
  
  x_scaled <- x - mean(x)
  
  N <- length(x_scaled)
  L <- max(abs(x_scaled)) * c
  
  phi <- PHI(N = N, M = M, L = L, x = x_scaled)
  
  spd_matern52 <- c(diagSPD_Matern52(
    alpha = 1, ls_noise, L, M
  ))
  
  
  phi_scaled <- sweep(phi, 2, spd_matern52^2, "*")
  k_matrix <- tcrossprod(phi_scaled, phi)
  tbl <- expand_grid(
    x_1 = x,
    x_2 = x
  ) |> 
    mutate(k_approx = as.vector(t(k_matrix)))
  
  return(tbl)
}

x <- seq(0, 365 * 3, by = 1)

c <- 1.25
M <- calculate_M_matern(c = c, length_scale = 7.0, x = x)

plot_data_matern <- tibble(
  ls = 3:120
) |> 
  pmap(function(ls) {
    
    approximated_kernel <- approx_kernel_matern52(x, ls, c, M)
    true_kernel <- true_kernel_matern52(x, ls)
    
    true_kernel |> 
      left_join(approximated_kernel, by = join_by(x_1, x_2)) |> 
      summarise(epsilon = sum(abs(k_approx - k_true)),
                auc = sum(k_true),
                error = epsilon / auc,
                ls = ls)
    
  }, .progress = TRUE) |> 
  bind_rows()

p_err <- ggplot(plot_data_matern) +
  
  geom_vline(xintercept = 7, linetype = "44") +
  geom_vline(xintercept = 90, linetype = "44") +
  
  geom_line(aes(x = ls, y = error)) +
  
  xlab("Length-scale") + ylab("Relative error") +
  
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100)) +
  coord_cartesian(xlim = c(0, NA)) +
  
  scale_y_log10() +
  
  plot_theme_paper +
  
  theme(panel.grid.major = element_gridline) +
  
  ggtitle("<b>B</b>. Relative error in the HSGP approximation of the Matern 5/2 kernel")


plot_data_ex <- tibble(k = c(4, 7, 14)) |> 
  pmap(function(k) {
    approximated_kernel <- approx_kernel_matern52(x, k, c, M)
    true_kernel <- true_kernel_matern52(x, k)
    
    true_kernel |> 
      left_join(approximated_kernel) |> 
      mutate(k = k)
  }) |> 
  
  bind_rows() |> 
  filter(x_1 == 100) |> 
  
  mutate(k_label = str_c("Length-scale = ", k),
         k_label = factor(k_label, str_c("Length-scale = ", unique(k))))


p_ex <- ggplot(plot_data_ex) +
  
  geom_vline(xintercept = 100, linetype = "44") +
  
  geom_line(aes(x = x_2, y = k_approx, colour = "Approximated kernel"),
            linewidth = 1.0, alpha = 0.8) +
  geom_line(aes(x = x_2, y = k_true, colour = "True kernel")) +
  
  coord_cartesian(xlim = c(75, 125)) +
  
  facet_wrap(~k_label, ncol = 3) +
  
  ggokabeito::scale_colour_okabe_ito(name = "Kernel", order = c(5, 3)) +
  
  plot_theme_paper +
  
  xlab("<i>x</i><sub>2</sub>") + ylab("Covariance") +
  
  theme(panel.background = element_facet_background) +
  
  ggtitle("<b>A</b>. Example Matern 5/2 kernel HSGP approximations (<i>x</i><sub>1</sub> = 100)")







true_kernel_periodic <- function(x, ls_periodic) {
  
  kernel <- function(tau, ls_periodic) {
    period  <- 365.25
    w0      <- 2 * pi / period
    
    rho_p_rad <- ls_periodic * (2 * pi / 365.25)
    # \sigma^2_\text{ann} \exp \left(-\frac{2\sin ^ 2(\omega_0 \tau / 2)}{\ell_\text{ann}^2}\right),
    return(exp(-(2 * sin(w0 * tau / 2)^2)/rho_p_rad^2))
  }
  
  tbl <- expand_grid(
    x_1 = x,
    x_2 = x
  ) |> 
    mutate(k_true = kernel(abs(x_1 - x_2), ls_periodic))
  
  return(tbl)
}

approx_kernel_periodic <- function(x, ls_periodic, M) {
  
  x_centered <- x - mean(x)
  N <- length(x_centered)
  
  period  <- 365.25
  w0      <- 2 * pi / period
  
  
  rho_p_rad <- ls_periodic * (2 * pi / 365.25)
  spd_periodic <- diagSPD_Periodic(alpha = 1.0, rho = rho_p_rad, M = M)
  
  phi_periodic <- PHI_periodic(N = N, M = M, w0 = w0, x = x_centered)
  
  a <- 1 / rho_p_rad^2
  alpha <- 1.0
  
  # Compute the flat m=0 offset using the exponentially scaled Bessel function
  m0_offset <- alpha^2 * besselI(x = a, nu = 0, expon.scaled = TRUE)
  
  phi_scaled <- sweep(phi_periodic, 2, spd_periodic^2, "*")
  k_matrix <- tcrossprod(phi_scaled, phi_periodic)
  tbl <- expand_grid(
    x_1 = x,
    x_2 = x
  ) |> 
    mutate(k_approx = as.vector(t(k_matrix)) + m0_offset)
  
  return(tbl)
}


x <- seq(0, 365 * 3, by = 1)

M_per <- 40


plot_data_periodic <- tibble(
  ls = 3:120
) |> 
  pmap(function(ls) {
    
    approximated_kernel <- approx_kernel_periodic(x, ls, M_per)
    true_kernel <- true_kernel_periodic(x, ls)
    
    true_kernel |> 
      left_join(approximated_kernel, by = join_by(x_1, x_2)) |> 
      summarise(epsilon = sum(abs(k_approx - k_true)),
                auc = sum(k_true),
                error = epsilon / auc,
                ls = ls)
    
  }, .progress = TRUE) |> 
  bind_rows()


p_err_per <- ggplot(plot_data_periodic) +
  
  geom_vline(xintercept = 14, linetype = "44") +
  geom_vline(xintercept = 90, linetype = "44") +
  
  geom_line(aes(x = ls, y = error)) +
  
  xlab("Length-scale") + ylab("Relative error") +
  
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100)) +
  coord_cartesian(xlim = c(0, NA)) +
  
  scale_y_log10() +
  
  plot_theme_paper +
  
  theme(panel.grid.major = element_gridline) +
  
  ggtitle("<b>C</b>. Relative error in the HSGP approximation of the periodic kernel")


p_ex / p_err / p_err_per


ggsave(
  "paper_figures/supp/f_hsgp.png",
  width = 8, height = 8,
  bg = "white"
)




