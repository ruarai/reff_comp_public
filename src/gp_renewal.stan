
// Using GP functions from EpiNow2
functions {
  #include gp_functions.stan
  
  #include qr_weekly.stan
}

data {
  int<lower=1> N;
  int<lower=1> N_data;
  
  array[N] int day_of_week;
  vector[N] school_holidays;
  
  vector[N] x;
  array[N_data] int y;
  
  
  int<lower=1> N_gi;
  vector[N_gi] gi_reversed;

  int<lower=1> M_annual;
  int<lower=1> M_noise; 
  
  real<lower=1.0> c_factor;
}

transformed data {
  matrix[N, 2 * M_annual] phi_annual = PHI_periodic(N, M_annual, 2 * pi() / 365.0, x);
  
  vector[N] x_centered = x - mean(x);
  real L = c_factor * max(abs(x_centered));
  
  matrix[N, M_noise] phi_noise = PHI(N, M_noise, L, x_centered);
  
  matrix[7, 6] weekly_Q = qr_weekly();
  
  
  // Infection time-series with burn-in period
  int N_total = N + N_gi;
}

parameters {
  real<lower=0> scaling_annual;
  real<lower=0> scaling_weekly;
  real<lower=0> scaling_noise;
  
  real<lower=0> ls_annual;
  real<lower=0> ls_noise;
  
  real<lower=0> phi;

  vector[2 * M_annual] beta_annual;
  vector[M_noise] beta_noise;
  
  vector[6] beta_weekly;
  
  real r_0;
  real log_I_0;
  
  real school_holiday_effect;
}

transformed parameters {
  vector[2 * M_annual] diag_SPD_annual = diagSPD_Periodic(scaling_annual, ls_annual, M_annual);
  vector[M_noise] diag_SPD_noise = diagSPD_Matern52(scaling_noise, ls_noise, L, M_noise);
  
  
  vector[N] f_annual = phi_annual * (diag_SPD_annual .* beta_annual);
  vector[N] f_noise  = phi_noise * (diag_SPD_noise .* beta_noise);
  
  vector[N] log_R = f_annual + f_noise + school_holiday_effect * school_holidays;
  
  // Implies each obs_offset_weekly_i ~ N(0, scaling_weekly)
  // with constraint of sum to zero
  vector[7] obs_offset_weekly = weekly_Q * (beta_weekly * scaling_weekly) * sqrt(7.0 / 6.0);
  
  vector[N] inf;
  
  {
    vector[N_total] inf_full = rep_vector(0.0, N_total);

    // Calculate initial infs as exp growth: I(t) = I_0 * exp(r * t)
    inf_full[1:N_gi] = exp(log_I_0 + r_0 * linspaced_vector(N_gi, 0, N_gi - 1));
    
    for(t in 1:N) {
      real R_t = exp(log_R[t]);
      
      real new_infections = R_t * dot_product(gi_reversed, segment(inf_full, t, N_gi));
      
      inf_full[t + N_gi] = new_infections;
    }
    
    inf = segment(inf_full, N_gi + 1, N);
  }
  
  vector[N] inf_weekly = inf .* exp(obs_offset_weekly[day_of_week]);
}

model {
  scaling_annual ~ exponential(1.0);
  scaling_weekly ~ exponential(1.0);
  scaling_noise  ~ exponential(1.0);
  
  school_holiday_effect ~ std_normal();
  
  // Inv-gamma priors derived from:
  // find_inv_gamma_params(7.0, 28.0)
  // find_inv_gamma_params(7.0, 90.0)
  // as in Betancourt approach.
  
  ls_noise ~ inv_gamma(11.83, 148.84);
  ls_annual ~ inv_gamma(3.86, 68.73);
  phi ~ gamma(2, 0.1);

  beta_annual ~ std_normal();
  beta_weekly ~ std_normal();
  beta_noise  ~ std_normal();
  
  log_I_0 ~ lognormal(3.2, 1);
  r_0 ~ normal(0, 0.2);
  
  y ~ neg_binomial_2(inf_weekly[1:N_data], phi);
}

generated quantities {
  array[N] int cases = neg_binomial_2_rng(inf_weekly, phi);
}
