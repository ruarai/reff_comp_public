// Using GP functions from EpiNow2
functions {
  #include gp_functions.stan
}


data {
  int<lower=1> N;
  int<lower=1> N_data;
  
  vector[N] x;
  
  array[N_data] int y;

  int<lower=1> M_annual;
  int<lower=1> M_noise; 
  
  real<lower=1.0> c_factor;
}

transformed data {
  matrix[N, 2 * M_annual] phi_annual = PHI_periodic(N, M_annual, 2 * pi() / 365.0, x);
  
  vector[N] x_centered = x - mean(x);
  real L = c_factor * max(abs(x_centered));
  
  matrix[N, M_noise] phi_noise = PHI(N, M_noise, L, x_centered);
}

parameters {
  real<lower=0> scaling_annual;
  
  real<lower=0> scaling_noise;
  
  real<lower=0> ls_annual;
  real<lower=0> ls_noise;
  
  real mu;
  real<lower=0> phi;
  
  vector[7] obs_offset_weekly;

  vector[2 * M_annual] beta_annual;
  vector[M_noise] beta_noise;
}

transformed parameters {
  vector[2 * M_annual] diag_SPD_annual = diagSPD_Periodic(scaling_annual, ls_annual, M_annual);
  vector[M_noise] diag_SPD_noise = diagSPD_Matern52(scaling_noise, ls_noise, L, M_noise);
  
  
  vector[N] f_annual = phi_annual * (diag_SPD_annual .* beta_annual);
  vector[N] f_noise  = phi_noise * (diag_SPD_noise .* beta_noise);
  
  vector[N] f_total = mu + f_annual + f_noise;
  
  for(i in 1:N) {
    f_total[i] += obs_offset_weekly[i % 7 + 1];
  }
}

model {
  scaling_annual ~ exponential(1.0);
  scaling_noise  ~ exponential(1.0);
  
  mu ~ std_normal();
  
  // Inv-gamma priors derived from:
  // find_inv_gamma_params(7.0, 28.0)
  // find_inv_gamma_params(14.0, 90.0)
  // find_inv_gamma_params(1.0, 3.0)
  // as in Betancourt approach.
  
  ls_noise ~ inv_gamma(11.83, 148.84);
  ls_annual ~ inv_gamma(6.81, 200.12);
  phi ~ gamma(2, 0.1);

  beta_annual ~ std_normal();
  beta_noise  ~ std_normal();
  
  obs_offset_weekly ~ std_normal();
  
  y ~ neg_binomial_2_log(f_total[1:N_data], phi);
}
