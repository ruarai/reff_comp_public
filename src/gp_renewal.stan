
functions {
  // Using GP functions from EpiNow2
  #include gp_functions.stan
  
  #include helper_functions.stan
  
}

data {
  int<lower=1> N;
  int<lower=1> N_loc;
  
  int<lower=1> N_data;
  
  int<lower=1> N_public_holiday_types;
  
  array[N] int day_of_week;
  array[N_loc, N] int school_holiday_type;
  array[N_loc, N] int public_holiday_type;
  
  vector[N] x;
  array[N_loc, N_data] int y;
  
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
  int N_school_holiday_types = 4;
  
  int N_total = N + N_delay;
}

parameters {
  vector<lower=0>[N_loc] scaling_annual;
  vector<lower=0>[N_loc] scaling_noise;
  
  vector<lower=0>[N_loc] scaling_weekly_asc;
  
  vector<lower=0>[N_loc] ls_annual;
  vector<lower=0>[N_loc] ls_noise;
  
  vector<lower=0>[N_loc] phi;

  matrix[2 * M_annual, N_loc] beta_annual;
  matrix[M_noise, N_loc] beta_noise;
  
  matrix[6, N_loc] beta_weekly_asc;
  
  vector[N_loc] r_0;
  vector[N_loc] log_I_0;
  
  // For school holidays and public holidays, type = 1 is no holiday
  // Only fit effects to type > 1
  vector[N_school_holiday_types] mu_school_Rt; 
  vector<lower=0>[N_school_holiday_types] sigma_school_Rt;
  
  vector[N_school_holiday_types] mu_school_asc;
  vector<lower=0>[N_school_holiday_types] sigma_school_asc;

  matrix[N_school_holiday_types, N_loc] school_holiday_Rt_z;
  matrix[N_school_holiday_types, N_loc] school_holiday_asc_z;
  
  matrix[N_public_holiday_types, N_loc] public_holiday_Rt_effect_raw;
  matrix[N_public_holiday_types, N_loc] public_holiday_asc_effect_raw;
  
  vector<lower=0, upper=21>[N_loc] delay_mean;
  vector<lower=0>[N_loc] delay_sd;
}

transformed parameters {
  matrix[N, N_loc] log_R;
  matrix[N, N_loc] latent_cases_weekly;
  
  matrix[N, N_loc] f_annual;
  matrix[N, N_loc] f_noise;
  
  matrix[N_school_holiday_types, N_loc] school_holiday_Rt_effect_raw;
  matrix[N_school_holiday_types, N_loc] school_holiday_asc_effect_raw;

  for (i in 1:N_school_holiday_types) {
    school_holiday_Rt_effect_raw[i, ] = mu_school_Rt[i] + sigma_school_Rt[i] * school_holiday_Rt_z[i, ];
    school_holiday_asc_effect_raw[i, ] = mu_school_asc[i] + sigma_school_asc[i] * school_holiday_asc_z[i, ];
  }
  
  for(ix_loc in 1:N_loc) {
    real ls_annual_radians = ls_annual[ix_loc] * (2 * pi() / 365.25);
    vector[2 * M_annual] diag_SPD_annual = diagSPD_Periodic(scaling_annual[ix_loc], ls_annual_radians, M_annual);
    vector[M_noise] diag_SPD_noise = diagSPD_Matern52(scaling_noise[ix_loc], ls_noise[ix_loc], L, M_noise);
      
    f_annual[:, ix_loc] = phi_annual * (diag_SPD_annual .* beta_annual[:, ix_loc]);
    f_noise[:, ix_loc]  = phi_noise * (diag_SPD_noise .* beta_noise[:, ix_loc]);
    
    vector[N_school_holiday_types + 1] school_Rt_effect = get_effect(school_holiday_Rt_effect_raw, N_school_holiday_types, ix_loc);
    vector[N_school_holiday_types + 1] school_asc_effect = get_effect(school_holiday_asc_effect_raw, N_school_holiday_types, ix_loc);
    vector[N_public_holiday_types + 1] public_Rt_effect = get_effect(public_holiday_Rt_effect_raw, N_public_holiday_types, ix_loc);
    vector[N_public_holiday_types + 1] public_asc_effect = get_effect(public_holiday_asc_effect_raw, N_public_holiday_types, ix_loc);
    
    vector[7] obs_offset_weekly = weekly_Q * (col(beta_weekly_asc, ix_loc) * scaling_weekly_asc[ix_loc]) * sqrt(7.0 / 6.0);
    
    log_R[:, ix_loc] = f_annual[:, ix_loc] + 
                       f_noise[:, ix_loc] + 
                       school_Rt_effect[school_holiday_type[ix_loc]] +
                       public_Rt_effect[public_holiday_type[ix_loc]];
    
    vector[N_delay] delay_dist = construct_delay_dist(delay_mean[ix_loc], delay_sd[ix_loc], N_delay);
    vector[N] inf_delayed = rep_vector(0.0, N);
    vector[N_total] inf_full = rep_vector(0.0, N_total);

    // Initial infections for this location
    inf_full[1:N_delay] = exp(log_I_0[ix_loc] + r_0[ix_loc] * linspaced_vector(N_delay, -N_delay + 1, 0));
    
    for(t in 1:N) {
      real R_t = exp(log_R[t, ix_loc]);
      real new_infections = R_t * dot_product(gi_reversed, segment(inf_full, t, N_delay));
      inf_full[t + N_delay] = new_infections;
    }
    
    for(k in 1:N_delay) {
      int lag = k - 1;
      inf_delayed += segment(inf_full, N_delay + 1 - lag, N) * delay_dist[k];
    }
    
    vector[N] latent_cases = exp(
      log(inf_delayed) + 
      school_asc_effect[school_holiday_type[ix_loc]] +
      public_asc_effect[public_holiday_type[ix_loc]]
    );
    
    latent_cases_weekly[:, ix_loc] = latent_cases .* exp(obs_offset_weekly[day_of_week]);
  }
}
model {
  scaling_annual ~ exponential(1.0);
  scaling_noise  ~ exponential(1.0);
  
  scaling_weekly_asc ~ exponential(1.0);
  
  mu_school_Rt ~ normal(0, 0.5);
  sigma_school_Rt ~ exponential(2);
  
  mu_school_asc ~ normal(0, 0.5);
  sigma_school_asc ~ exponential(2);

  to_vector(school_holiday_Rt_z) ~ std_normal();
  to_vector(school_holiday_asc_z) ~ std_normal();
  
  to_vector(public_holiday_Rt_effect_raw) ~ normal(0, 0.25);
  to_vector(public_holiday_asc_effect_raw) ~ normal(0, 0.25);
  
  ls_noise ~ inv_gamma(3.858, 68.729);
  ls_annual ~ inv_gamma(6.809, 200.117);
  phi ~ gamma(2, 0.1);

  to_vector(beta_annual) ~ std_normal();
  to_vector(beta_noise)  ~ std_normal();
  
  to_vector(beta_weekly_asc) ~ std_normal();
  
  log_I_0 ~ normal(3.2, 1);
  r_0 ~ normal(0, 0.2);
  
  delay_mean ~ lognormal(1, 1);
  delay_sd ~ lognormal(0, 1);
  
  for (ix_loc in 1:N_loc) {
    y[ix_loc, :] ~ neg_binomial_2(latent_cases_weekly[1:N_data, ix_loc], phi[ix_loc]);
  }
}
generated quantities {
  array[N_loc, N] int cases; 
  
  for (ix_loc in 1:N_loc) {
    cases[ix_loc] = neg_binomial_2_rng(latent_cases_weekly[:, ix_loc], phi[ix_loc]);
  }
}
