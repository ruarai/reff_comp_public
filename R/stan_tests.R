

library(cmdstanr)
file.remove("src/gp_renewal")
model <- cmdstan_model("src/gp_renewal.stan")

x <- model$expose_functions(global = TRUE, verbose = FALSE)


plot(construct_delay_dist(1.5, 1.0, 21))


sum(construct_delay_dist(0.1, 1.5))


fit <- model$pathfinder(
  seed = 4, data = data_list,
  max_lbfgs_iters = 5000,
  num_elbo_draws = 30,
  tol_rel_grad = 1e-10, tol_rel_obj = 1e-10,
  history_size = 40
)

