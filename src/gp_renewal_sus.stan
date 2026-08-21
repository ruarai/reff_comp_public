functions {
  // Using GP functions from EpiNow2
  #include gp_functions.stan
  
  #include helper_functions.stan
}

data {
  int<lower=1> N;
  int<lower=1> N_data;
  
  int<lower=1> N_public_holiday_types;
  int<lower=1> N_school_holiday_types;
  
  array[N] int day_of_week;
  array[N] int school_holiday_type;
  array[N] int public_holiday_type;
  vector[N] x;
  array[N_data] int y;
  
  int<lower=1> max_window_size;
  vector[max_window_size] past_infections;
  
  vector[21] gi_reversed;

  int<lower=1> M_annual;
  int<lower=1> M_noise; 
  
  real<lower=1.0> c_factor;
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
  // These are standard deviations of the GPs
  real<lower=0> scaling_annual; // \sigma_\text{ann}
  real<lower=0> scaling_noise; // \sigma_\text{noise}
  
  // Standard deviation of day of week term
  real<lower=0> scaling_weekly_asc;
  
  // Length-scales
  real<lower=0, upper=365> ls_annual;
  real<lower=0, upper=120> ls_noise;
  
  // logR(t) intercept term
  real log_R_0;
  
  // NB dispersion term
  real<lower=0, upper=200> phi;

  vector[2 * M_annual] beta_annual;
  vector[M_noise] beta_noise;
  
  sum_to_zero_vector[7] beta_offset_weekly;
  
  real r_0;
  real log_I_0;
  
  real<lower=0.0> cumulative_infections_effect;
  
  // For school holidays and public holidays, type = 1 is no holiday
  // Only fit effects to type > 1
  vector[N_school_holiday_types] school_holiday_Rt_effect_raw;
  vector[N_school_holiday_types] school_holiday_asc_effect_raw;
  
  vector[N_public_holiday_types] public_holiday_Rt_effect_raw;
  vector[N_public_holiday_types] public_holiday_asc_effect_raw;
  
  real<lower=3, upper=21> delay_mean;
  real<lower=1.0> delay_sd;
  
  real<lower=0, upper=max_window_size> sus_window_size;
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
  
  // Calculate the GPs as the weighting vectors from the basis functions
  // and the weighting vectors
  vector[N] f_annual = phi_annual * (diag_SPD_annual .* beta_annual);
  vector[N] f_noise  = phi_noise  * (diag_SPD_noise .* beta_noise);
  
  // Add a zero to the start of the school and public holiday effects
  // so that we can index into them using the holiday_type vectors
  // with type == 1 having an effect of zero
  vector[N_school_holiday_types + 1] school_holiday_Rt_effect = rep_vector(0.0, N_school_holiday_types + 1);
  school_holiday_Rt_effect[2:(N_school_holiday_types + 1)] = school_holiday_Rt_effect_raw;
  
  vector[N_public_holiday_types + 1] public_holiday_Rt_effect = rep_vector(0.0, N_public_holiday_types + 1);
  public_holiday_Rt_effect[2:(N_public_holiday_types + 1)] = public_holiday_Rt_effect_raw;
  
  // Calculate the value of log R(t) across t = 1 .. N
  vector[N] log_R = log_R_0 +
    f_annual +
    f_noise + 
    school_holiday_Rt_effect[school_holiday_type] +
    public_holiday_Rt_effect[public_holiday_type];
  
  // Calculate the delay distribution between infection onset and reporting
  vector[K] delay_dist = construct_delay_dist(delay_mean, delay_sd, K);
  
  // Number of new infections by time of infection onset
  vector[N] inf; 
  
  // Number of new infections convolved by delay_dist 
  // but not affected by any reporting/ascertainment effects
  vector[N] inf_delayed;
  
  vector[N] cumulative_effect;
  
  {
    // The full time-series of infections including M days prior to t = 1
    vector[N_total] inf_full = rep_vector(0.0, N_total);

    // Calculate initial infections as exp growth: I(t) = I_0 * exp(r * t)
    inf_full[1:K] = exp(log_I_0 + r_0 * linspaced_vector(K, -K + 1, 0));
    
    real running_inf = sum(past_infections);
    
    for(t in 1:N) {
      cumulative_effect[t] = cumulative_infections_effect * running_inf;
      
      real R_t = exp(log_R[t] - cumulative_effect[t]);
      
      real new_infections = R_t * dot_product(gi_reversed, segment(inf_full, t, K));
      inf_full[t + K] = new_infections;
      
      // Update the rolling sum for the next iteration
      running_inf += new_infections;
      
      if (t <= N_sus_window) {
        running_inf -= past_infections[t];
      } else {
        running_inf -= inf_full[t - N_sus_window + K];
      }
    }
    
    // We only care about t = 1..N, so drop the first K days
    inf = segment(inf_full, K + 1, N);
    
    // Calculate the delayed infections by creating an N x K matrix
    // where each column k is infections delayed by k - 1 days
    matrix[N, K] inf_lagged;
    for(k in 1:K) {
       inf_lagged[:, k] = segment(inf_full, K + 1 - (k - 1), N);
    }
    
    // The above construction allows us to calculate the full convolution
    // in a single matrix multiplication:
    inf_delayed = inf_lagged * delay_dist;
    
    inf_delayed = fmax(0.0, inf_delayed + 1e-5);
  }
  
  // As above, create vectors allowing for indexing of ascertainment effect by holiday_type
  vector[N_school_holiday_types + 1] school_holiday_asc_effect = rep_vector(0.0, N_school_holiday_types + 1);
  school_holiday_asc_effect[2:(N_school_holiday_types + 1)] = school_holiday_asc_effect_raw;
  
  vector[N_public_holiday_types + 1] public_holiday_asc_effect = rep_vector(0.0, N_public_holiday_types + 1);
  public_holiday_asc_effect[2:(N_public_holiday_types + 1)] = public_holiday_asc_effect_raw;
  
  // Calculate the latent cases time-series on the log scale with asc. effects
  // (but not day of week effects, as it is nice to look at them separately)
  vector[N] latent_cases_log = log(inf_delayed) + 
    school_holiday_asc_effect[school_holiday_type] +
    public_holiday_asc_effect[public_holiday_type];

  // Calculate the day-of-week effect (non-centred parametrisation)
  vector[7] obs_offset_weekly = beta_offset_weekly * scaling_weekly_asc;
  
  // Add in the day-of-week effect
  vector[N] latent_cases_weekly_log = obs_offset_weekly[day_of_week] + latent_cases_log;
}

model {
  scaling_annual ~ exponential(1.0);
  scaling_noise   ~ exponential(1.0);
  
  scaling_weekly_asc ~ exponential(1.0);
  
  school_holiday_Rt_effect_raw ~ std_normal();
  school_holiday_asc_effect_raw ~ std_normal();
  
  public_holiday_Rt_effect_raw ~ normal(0, 0.25);
  public_holiday_asc_effect_raw ~ normal(0, 0.25);
  
  cumulative_infections_effect ~ normal(0, 0.001);
  
  // Inv-gamma priors derived from:
  // find_inv_gamma_params(7.0, 90.0)
  // find_inv_gamma_params(14.0, 90.0)
  // as in Betancourt approach (Section 3.2.3)
  // https://betanalpha.github.io/assets/case_studies/gaussian_processes.html
  
  ls_noise ~ inv_gamma(3.858, 68.729);
  ls_annual ~ inv_gamma(6.809, 200.117);
  phi ~ gamma(2, 0.1);

  beta_annual ~ std_normal();
  beta_noise   ~ std_normal();
  
  log_R_0 ~ normal(0, 0.5);
  
  beta_offset_weekly ~ std_normal();
  
  log_I_0 ~ normal(3.2, 1);
  r_0 ~ normal(0, 0.2);
  
  delay_mean ~ normal(9, 2.5);
  delay_sd ~ lognormal(0, 0.5);
  
  // Data likelihood
  y ~ neg_binomial_2_log(latent_cases_weekly_log[1:N_data], phi);
}

generated quantities {
  // Calculate out case incidence for backcasting and forecasting
  array[N] int cases = neg_binomial_2_log_rng(latent_cases_weekly_log, phi);
}