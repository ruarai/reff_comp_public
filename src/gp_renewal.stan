
functions {
  // Using GP functions from EpiNow2
  #include gp_functions.stan
  
  #include qr_weekly.stan
  #include delay_dist.stan
  
}

data {
  int<lower=1> N;
  int<lower=1> N_data;
  
  int<lower=1> N_public_holiday_types;
  
  array[N] int day_of_week;
  array[N] int school_holiday_type;
  array[N] int public_holiday_type;
  
  vector[N] x;
  array[N_data] int y;
  
  vector[21] gi_reversed;

  int<lower=1> M_annual;
  int<lower=1> M_noise; 
  
  real<lower=1.0> c_factor;
}

transformed data {
  matrix[N, 2 * M_annual] phi_annual = PHI_periodic(N, M_annual, 2 * pi() / 365.25, x);
  
  vector[N] x_centered = x - mean(x);
  real L = c_factor * max(abs(x_centered));
  
  matrix[N, M_noise] phi_noise = PHI(N, M_noise, L, x_centered);
  
  matrix[7, 6] weekly_Q = qr_weekly();
  
  int N_delay = 21;
  
  int N_total = N + N_delay;
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
  
  // For school holidays and public holidays, type = 1 is no holiday
  // Only fit effects to type > 1
  vector[4] school_holiday_Rt_effect_raw;
  vector[4] school_holiday_asc_effect_raw;
  
  real public_holiday_Rt_mean;
  real public_holiday_asc_mean;
  real<lower=0> public_holiday_Rt_sd;
  real<lower=0> public_holiday_asc_sd;
  
  vector[N_public_holiday_types - 1] public_holiday_Rt_effect_raw;
  vector[N_public_holiday_types - 1] public_holiday_asc_effect_raw;
  
  real<lower=0, upper=21> delay_mean;
  real<lower=0> delay_sd;
}

transformed parameters {
  real ls_annual_rad = ls_annual * (2 * pi() / 365.25);

  // Use the radian lengthscale for the SPD
  vector[2 * M_annual] diag_SPD_annual = diagSPD_Periodic(scaling_annual, ls_annual_rad, M_annual);
  
  vector[M_noise] diag_SPD_noise = diagSPD_Matern52(scaling_noise, ls_noise, L, M_noise);
  
  
  vector[N] f_annual = phi_annual * (diag_SPD_annual .* beta_annual);
  vector[N] f_noise  = phi_noise * (diag_SPD_noise .* beta_noise);
  
  vector[5] school_holiday_Rt_effect = rep_vector(0.0, 5);
  school_holiday_Rt_effect[2:5] = school_holiday_Rt_effect_raw;
  
  vector[N_public_holiday_types] public_holiday_Rt_effect = rep_vector(0.0, N_public_holiday_types);
  public_holiday_Rt_effect[2:N_public_holiday_types] = 
    public_holiday_Rt_mean + public_holiday_Rt_sd * public_holiday_Rt_effect_raw;
  
  vector[N] log_R = f_annual + 
    f_noise + 
    school_holiday_Rt_effect[school_holiday_type] +
    public_holiday_Rt_effect[public_holiday_type];
  
  // Implies each obs_offset_weekly_i ~ N(0, scaling_weekly)
  // with constraint of sum to zero
  vector[7] obs_offset_weekly = weekly_Q * (beta_weekly * scaling_weekly) * sqrt(7.0 / 6.0);
  
  vector[N_delay] delay_dist = construct_delay_dist(delay_mean, delay_sd, N_delay);
  
  vector[N] inf;
  vector[N] inf_delayed = rep_vector(0.0, N);
  
  {
    vector[N_total] inf_full = rep_vector(0.0, N_total);

    // Calculate initial infs as exp growth: I(t) = I_0 * exp(r * t)
    inf_full[1:N_delay] = exp(log_I_0 + r_0 * linspaced_vector(N_delay, 0, N_delay - 1));
    
    for(t in 1:N) {
      real R_t = exp(log_R[t]);
      
      real new_infections = R_t * dot_product(gi_reversed, segment(inf_full, t, N_delay));
      
      inf_full[t + N_delay] = new_infections;
    }
    
    inf = segment(inf_full, N_delay + 1, N);
    
    for(k in 1:N_delay) {
      int lag = k - 1;
      inf_delayed += segment(inf_full, N_delay + 1 - lag, N) * delay_dist[k];
    }
  }
  
  vector[5] school_holiday_asc_effect = rep_vector(0.0, 5);
  school_holiday_asc_effect[2:5] = school_holiday_asc_effect_raw;
  
  vector[N_public_holiday_types] public_holiday_asc_effect = rep_vector(0.0, N_public_holiday_types);
  public_holiday_asc_effect[2:N_public_holiday_types] = 
    public_holiday_asc_mean + public_holiday_asc_sd * public_holiday_asc_effect_raw;
  
  vector[N] latent_cases = exp(
    log(inf_delayed) + 
    obs_offset_weekly[day_of_week] +
    school_holiday_asc_effect[school_holiday_type] +
    public_holiday_asc_effect[public_holiday_type]
  );
}

model {
  scaling_annual ~ exponential(1.0);
  scaling_weekly ~ exponential(1.0);
  scaling_noise  ~ exponential(1.0);
  
  school_holiday_Rt_effect_raw ~ std_normal();
  school_holiday_asc_effect_raw ~ std_normal();
  
  public_holiday_Rt_mean ~ std_normal();
  public_holiday_asc_mean ~ std_normal();
  public_holiday_Rt_sd ~ normal(0, 0.5);
  public_holiday_asc_sd ~ normal(0, 0.5);
  
  public_holiday_Rt_effect_raw ~ std_normal();
  public_holiday_asc_effect_raw ~ std_normal();
  
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
  
  delay_mean ~ lognormal(1, 1);
  delay_sd ~ lognormal(0, 1);
  
  
  
  y ~ neg_binomial_2(latent_cases[1:N_data], phi);
}

generated quantities {
  array[N] int cases = neg_binomial_2_rng(latent_cases, phi);
}
