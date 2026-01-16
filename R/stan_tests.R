

file.remove("src/gp_renewal")
model <- cmdstan_model("src/gp_renewal.stan")

x <- model$expose_functions(global = TRUE, verbose = FALSE)


plot(construct_delay_dist(15, 1.5))


sum(construct_delay_dist(0.1, 1.5))


fit <- model$pathfinder(
  seed = 3, data = data_list,
  max_lbfgs_iters = 2000
)
