

functions {
  // Using GP functions from EpiNow2
  #include gp_functions.stan
  
  #include helper_functions.stan
}

data {
  int<lower=1> N;
  int<lower=1> N_data;
  int<lower=1> N_loc;
  
  int<lower=1> N_public_holiday_types;
  int<lower=1> N_school_holiday_types;
  
  array[N] int day_of_week;
  array[N_loc, N] int school_holiday_type;
  array[N_loc, N] int public_holiday_type;
  vector[N] x;
  array[N_loc, N_data] int y;
  
  int<lower=1> N_years;
  array[N_years] int year_start;
  array[N_years] int year_end;
  
  vector[21] gi_reversed;

  int<lower=1> M_annual;
  int<lower=1> M_noise; 
  
  real<lower=1.0> c_factor;
  
  real<lower=0.0> sigma_log10_total_inf;
}

transformed data {
  // Construct the basis functions for HSGP
  vector[N] x_centered = x - mean(x);
  real L = c_factor * max(abs(x_centered));
  
  matrix[N, M_noise] phi_noise = PHI(N, M_noise, L, x_centered);
  matrix[N, 2 * M_annual] phi_annual = PHI_periodic(N, M_annual, 2 * pi() / 365.25, x);
  
  int K = 21;
  int N_total = N + K;
}

parameters {
  real<lower=0> scaling_annual;
  real<lower=0> scaling_noise;
  
  real<lower=0, upper=365> ls_annual;
  real<lower=0, upper=120> ls_noise;
  
  real<lower=0> scaling_weekly_asc;
  
  vector[N_loc] log_R_0;
  vector<lower=0, upper=200>[N_loc] phi;

  array[N_loc] vector[2 * M_annual] beta_annual;
  array[N_loc] vector[M_noise] beta_noise;
  
  array[N_loc] sum_to_zero_vector[7] beta_offset_weekly;
  
  vector[N_loc] r_0;
  vector[N_loc] log_I_0;
  
  // For school holidays and public holidays, type = 1 is no holiday
  // Only fit effects to type > 1
  array[N_loc] vector[N_school_holiday_types] school_holiday_Rt_effect_raw;
  
  array[N_loc] vector[N_public_holiday_types] public_holiday_Rt_effect_raw;
  
  vector<lower=0, upper=21>[N_loc] delay_mean;
  vector<lower=1.0>[N_loc] delay_sd;
  
  vector[N_years] log10_total_inf_year_mean;
  vector[N_loc - 1] log10_total_inf_offset_raw;
}

transformed parameters {
  // Calculate length-scale for annual GP in radians
  real ls_annual_rad = ls_annual * (2 * pi() / 365.25);
    
  // Calculate the weighting (SPD) vector for the periodic annual GP 
  // from ls_annual_rad and scaling_annual
  vector[2 * M_annual] diag_SPD_annual = diagSPD_Periodic(scaling_annual, ls_annual_rad, M_annual);
  
    
  // Calculate the weighting (SPD) vector for the Matern 5/2 noise GP
  // from ls_noise and scaling_noise
  vector[M_noise] diag_SPD_noise = diagSPD_Matern52(scaling_noise, ls_noise, L, M_noise);
  
  // Number of new infections by time of infection onset
  array[N_loc] vector[N] inf; 
  
  // Number of new infections convolved by delay_dist 
  // but not affected by any reporting/ascertainment effects
  array[N_loc] vector[N] inf_delayed;
  
  array[N_loc] vector[N] latent_cases_weekly_log;
  
  array[N_loc] vector[N_years] log10_total_inf_year; 
  
  vector[N_loc] log10_total_inf_offset = rep_vector(0.0, N_loc);
  log10_total_inf_offset[2:N_loc] = log10_total_inf_offset_raw;
  
  for(i in 1:N_loc) {
    // Calculate the GPs as the weighting vectors from the basis functions
    // and the weighting vectors
    vector[N] f_annual = phi_annual * (diag_SPD_annual .* beta_annual[i]);
    vector[N] f_noise  = phi_noise * (diag_SPD_noise .* beta_noise[i]);
    
    // Add a zero to the start of the school and public holiday effects
    // so that we can index into them using the holiday_type vectors
    // with type == 1 having an effect of zero
    vector[N_school_holiday_types + 1] school_holiday_Rt_effect = rep_vector(0.0, N_school_holiday_types + 1);
    school_holiday_Rt_effect[2:(N_school_holiday_types + 1)] = school_holiday_Rt_effect_raw[i];
    
    vector[N_public_holiday_types + 1] public_holiday_Rt_effect = rep_vector(0.0, N_public_holiday_types + 1);
    public_holiday_Rt_effect[2:(N_public_holiday_types + 1)] = public_holiday_Rt_effect_raw[i];
    
    // Calculate the value of log R(t) across t = 1 .. N
    vector[N] log_R = log_R_0[i] +
      f_annual +
      f_noise + 
      school_holiday_Rt_effect[school_holiday_type[i]] +
      public_holiday_Rt_effect[public_holiday_type[i]];
    
    // Calculate the delay distribution between infection onset and reporting
    vector[K] delay_dist = construct_delay_dist(delay_mean[i], delay_sd[i], K);
    
    {
      // The full time-series of infections including M days prior to t = 1
      vector[N_total] inf_full = rep_vector(0.0, N_total);
  
      // Calculate initial infections as exp growth: I(t) = I_0 * exp(r * t)
      inf_full[1:K] = exp(log_I_0[i] + r_0[i] * linspaced_vector(K, -K + 1, 0));
      
      for(t in 1:N) {
        // Perform the renewal equation step
        real R_t = exp(log_R[t]);
        
        // Use a dot-product here to make the sum over the generation interval
        // a single auto-diff operation
        real new_infections = R_t * dot_product(gi_reversed, segment(inf_full, t, K));
        
        inf_full[t + K] = new_infections;
      }
      
      // We only care about t = 1..N, so drop the first K days
      inf[i] = segment(inf_full, K + 1, N);
      
      // Calculate the delayed infections by creating an N x K matrix
      // where each column k is infections delayed by k - 1 days
      matrix[N, K] inf_lagged;
      for(k in 1:K) {
         inf_lagged[:, k] = segment(inf_full, K + 1 - (k - 1), N);
      }
      
      // The above construction allows us to calculate the full convolution
      // in a single matrix multiplication:
      inf_delayed[i] = inf_lagged * delay_dist;
    }
    
    for (j in 1:N_years) { 
    	log10_total_inf_year[i][j] = log10(sum(inf[i][year_start[j]:year_end[j]])); 
    }
    
    vector[N] latent_cases_log = log(inf_delayed[i] + 1e-4);
    
    // Calculate the day-of-week effect (non-centred parametrisation)
    vector[7] obs_offset_weekly = beta_offset_weekly[i] * scaling_weekly_asc;
    
    
    
    // Add in the day-of-week effect
    latent_cases_weekly_log[i] = obs_offset_weekly[day_of_week] + latent_cases_log;
  }
}

model {
  scaling_annual ~ exponential(1.0);
  scaling_noise  ~ exponential(1.0);
  
  scaling_weekly_asc ~ exponential(1.0);
  
  // Inv-gamma priors derived from:
  // find_inv_gamma_params(7.0, 90.0)
  // find_inv_gamma_params(14.0, 90.0)
  // as in Betancourt approach (Section 3.2.3)
  // https://betanalpha.github.io/assets/case_studies/gaussian_processes.html
  
  ls_noise ~ inv_gamma(3.858, 68.729);
  ls_annual ~ inv_gamma(6.809, 200.117);
  phi ~ gamma(2, 0.1);
  
  log_R_0 ~ normal(0, 0.5);
  
  
  log_I_0 ~ normal(3.2, 1);
  r_0 ~ normal(0, 0.2);
  
  delay_mean ~ lognormal(2.0, 0.5);
  delay_sd ~ lognormal(1, 0.5);
  
  log10_total_inf_year_mean ~ normal(4, 3);
  log10_total_inf_offset_raw ~ normal(0, 2);
  
  for(i in 1:N_loc) {
      
    beta_offset_weekly[i] ~ normal(0, 1);
    
    school_holiday_Rt_effect_raw[i] ~ normal(0, 1);
    public_holiday_Rt_effect_raw[i] ~ normal(0, 0.25);
  
    beta_annual[i] ~ normal(0, 1);
    beta_noise[i]  ~ normal(0, 1);
  
    log10_total_inf_year[i] ~ normal(log10_total_inf_year_mean + log10_total_inf_offset[i], sigma_log10_total_inf);
    
    
    // Data likelihood
    y[i] ~ neg_binomial_2_log(latent_cases_weekly_log[i][1:N_data], phi[i]);
  }
}

generated quantities {
  // Calculate out case incidence for backcasting and forecasting
  array[N_loc, N] int cases;
  
  for(i in 1:N_loc) {
    cases[i] = neg_binomial_2_log_rng(latent_cases_weekly_log[i], phi[i]);
  }
}
