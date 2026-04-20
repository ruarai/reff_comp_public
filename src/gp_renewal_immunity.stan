

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
  array[N] int day_of_year;
  array[N] int school_holiday_type;
  array[N] int public_holiday_type;
  vector[N] x;
  array[N_data] int y;
  
  vector[366] humidity;
  
  
  vector[21] gi_reversed;

  // int<lower=1> M_annual;
  int<lower=1> M_noise; 
  
  real<lower=1.0> c_factor;
  
  real<lower=0.0> population_size;
  
  real init_factor;
  real log_R_0;
  
  int k;
}

transformed data {
  // matrix[N, 2 * M_annual] phi_annual = PHI_periodic(N, M_annual, 2 * pi() / 365.25, x);
  
  vector[N] x_centered = x - mean(x);
  real L = c_factor * max(abs(x_centered));
  
  matrix[N, M_noise] phi_noise = PHI(N, M_noise, L, x_centered);
  
  int N_delay = 21;
  
  int N_total = N + N_delay;
}

parameters {
  real<lower=0.0> scaling_noise;
  
  real<lower=0.0> scaling_weekly_asc;
  
  real<lower=0.0> ls_noise;
  
  real<lower=0.0> phi;

  vector[M_noise] beta_noise;
  
  sum_to_zero_vector[7] beta_offset_weekly;
  
  real r_0;
  real log_I_0;
  
  real scaling_humidity;
  
  // For school holidays and public holidays, type = 1 is no holiday
  // Only fit effects to type > 1
  vector[N_school_holiday_types] school_holiday_Rt_effect_raw;
  vector[N_school_holiday_types] school_holiday_asc_effect_raw;
  
  vector[N_public_holiday_types] public_holiday_Rt_effect_raw;
  vector[N_public_holiday_types] public_holiday_asc_effect_raw;
  
  real<lower=0, upper=21> delay_mean;
  real<lower=1.0> delay_sd;
  
  vector<lower=0, upper=1>[k] initial_immunity_prop;
  real<lower=0,upper=0.01> alpha;
  
  real<lower=0.225, upper=2> immunity_half_life_200;
}

transformed parameters {
  vector[M_noise] diag_SPD_noise = diagSPD_Matern52(scaling_noise, ls_noise, L, M_noise);
  
  
  vector[N] f_noise  = phi_noise * (diag_SPD_noise .* beta_noise);
  
  vector[N_school_holiday_types + 1] school_holiday_Rt_effect = rep_vector(0.0, N_school_holiday_types + 1);
  school_holiday_Rt_effect[2:(N_school_holiday_types + 1)] = school_holiday_Rt_effect_raw;
  
  vector[N_public_holiday_types + 1] public_holiday_Rt_effect = rep_vector(0.0, N_public_holiday_types + 1);
  public_holiday_Rt_effect[2:(N_public_holiday_types + 1)] = public_holiday_Rt_effect_raw;
  
  vector[N] humidity_effect = humidity[day_of_year] * scaling_humidity;
  
  vector[N] log_R = log_R_0 +
    humidity_effect +
    f_noise +
    school_holiday_Rt_effect[school_holiday_type] +
    public_holiday_Rt_effect[public_holiday_type];
  
  vector[7] obs_offset_weekly = beta_offset_weekly * scaling_weekly_asc;
  
  vector[N_delay] delay_dist = construct_delay_dist(delay_mean, delay_sd, N_delay);
  
  vector[N] inf;
  vector[N] inf_delayed;
  
  vector[N] R_t;
  
  real immunity_half_life = immunity_half_life_200 * 200;
  real lambda = exp(k * log(0.5) / immunity_half_life);
  
  vector[k] accumulated_immunity = initial_immunity_prop * population_size * init_factor;
  vector[k] immunity_next;
  
  vector[N] immune_effect;
  
  {
    vector[N_total] inf_full = rep_vector(0.0, N_total);
    
    // TODO skip this 
    vector[N_delay] rev_delay_dist = reverse(delay_dist);

    // Calculate initial infs as exp growth: I(t) = I_0 * exp(r * t)
    inf_full[1:N_delay] = exp(log_I_0 + r_0 * linspaced_vector(N_delay, -N_delay + 1, 0));
    
    for(t in 1:N) {
      immune_effect[t] = alpha * sum(accumulated_immunity);
      R_t[t] = exp(log_R[t] - immune_effect[t]);
      
      real new_infections = R_t[t] * dot_product(gi_reversed, segment(inf_full, t, N_delay));
      
      inf_full[t + N_delay] = new_infections;
    
      immunity_next[1:(k - 1)] = accumulated_immunity[1:(k - 1)] * lambda 
                                 + accumulated_immunity[2:k] * (1 - lambda);
                              
      immunity_next[k] = accumulated_immunity[k] * lambda + new_infections;
      accumulated_immunity = immunity_next;
    }
    
    inf = segment(inf_full, N_delay + 1, N);
    
    for(t in 1:N) {
      inf_delayed[t] = dot_product(inf_full[t + 1 : t + N_delay], rev_delay_dist);
    }
  }
  
  vector[N_school_holiday_types + 1] school_holiday_asc_effect = rep_vector(0.0, N_school_holiday_types + 1);
  school_holiday_asc_effect[2:(N_school_holiday_types + 1)] = school_holiday_asc_effect_raw;
  
  vector[N_public_holiday_types + 1] public_holiday_asc_effect = rep_vector(0.0, N_public_holiday_types + 1);
  public_holiday_asc_effect[2:(N_public_holiday_types + 1)] = public_holiday_asc_effect_raw;
  
  vector[N] latent_cases_log = log(inf_delayed + 1e-3) + 
    school_holiday_asc_effect[school_holiday_type] +
    public_holiday_asc_effect[public_holiday_type];
  
  vector[N] latent_cases_weekly_log = obs_offset_weekly[day_of_week] + latent_cases_log;
}

model {
  scaling_noise  ~ normal(0, 0.1);
  
  scaling_weekly_asc ~ exponential(1.0);
  
  school_holiday_Rt_effect_raw ~ normal(0, 0.5);
  school_holiday_asc_effect_raw ~ normal(0, 0.5);
  
  public_holiday_Rt_effect_raw ~ normal(0, 0.25);
  public_holiday_asc_effect_raw ~ normal(0, 0.25);
  
  scaling_humidity ~ normal(0, 0.5);
  
  // Inv-gamma priors derived from:
  // find_inv_gamma_params(7.0, 90.0)
  // as in Betancourt approach.
  
  ls_noise ~ inv_gamma(3.858, 68.729);
  phi ~ gamma(2, 0.1);

  beta_noise  ~ std_normal();
  
  initial_immunity_prop ~ beta(1, 10);
  alpha ~ lognormal(-13, 1.0);
  
  
  // find_inv_gamma_params(60 / 200, 365 * 2 / 200)
  immunity_half_life_200 ~ inv_gamma(4.008, 3.017);
  
  beta_offset_weekly ~ std_normal();
  
  log_I_0 ~ normal(3.2, 1);
  r_0 ~ normal(0, 0.2);
  
  delay_mean ~ lognormal(2.0, 0.5);
  delay_sd ~ lognormal(1, 0.5);
  
  y ~ neg_binomial_2_log(latent_cases_weekly_log[1:N_data], phi);
}

generated quantities {
  array[N] int cases = neg_binomial_2_log_rng(latent_cases_weekly_log, phi);
}
