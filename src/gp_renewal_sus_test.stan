functions {
  
/**
  Stan code taken from EpiNow2 (Sam Abbott, Sebastian Funk)
  
  https://github.com/epiforecasts/EpiNow2/blob/main/inst/stan/functions/gaussian_process.stan

*/

/**
 * These functions implement approximate Gaussian processes for Stan using
 * Hilbert space methods. The functions are based on the following:
 * - https://avehtari.github.io/casestudies/Motorcycle/motorcycle_gpcourse.html (Section 4)
 * - https://doi.org/10.1007/s11222-022-10167-2
 */

/**
  * Spectral density for periodic kernel
  *
  * @param alpha Scaling parameter
  * @param rho Length scale parameter
  * @param M Number of basis functions
  * @param indices Pre-computed 1:M vector sequence
  * @return A vector of spectral densities
  *
  * @ingroup estimates_smoothing
  */
vector diagSPD_Periodic(real alpha, real rho, int M, vector indices) {
  real a = inv_square(rho);
  vector[M] q = exp(
    log(alpha) + 0.5 *
      (log(2) - a + to_vector(log_modified_bessel_first_kind(indices, a)))
  );
  return append_row(q, q);
}

/**
  * Spectral density for 5/2 Matern kernel
  *
  * @param alpha Scaling parameter
  * @param rho Length scale parameter
  * @param L Length of the interval
  * @param M Number of basis functions
  * @param indices Pre-computed 1:M vector sequence
  * @return A vector of spectral densities
  *
  * @ingroup estimates_smoothing
  */
vector diagSPD_Matern52(real alpha, real rho, real L, int M, vector indices) {
  real factor = 16 * pow(sqrt(5) / rho, 5);
  vector[M] denom =
    3 * pow((sqrt(5) / rho)^2 + pow((pi() / 2 / L) * indices, 2), 3);
  return alpha * sqrt(factor * inv(denom));
}

matrix PHI(int N, int M, real L, vector x, vector indices) {
  matrix[N, M] phi = sin(
    diag_post_multiply(
      rep_matrix(pi() / (2 * L) * (x + L), M), indices
    )
  ) / sqrt(L);
  return phi;
}

/**
  * Basis functions for periodic Gaussian Process
  *
  * @param N Number of data points
  * @param M Number of basis functions
  * @param w0 Fundamental frequency
  * @param x Vector of input data
  * @param indices Pre-computed 1:M vector sequence
  * @return A matrix of basis functions
  *
  * @ingroup estimates_smoothing
  */
matrix PHI_periodic(int N, int M, real w0, vector x, vector indices) {
  matrix[N, M] mw0x = diag_post_multiply(
    rep_matrix(w0 * x, M), indices
  );
  return append_col(cos(mw0x), sin(mw0x));
}

vector diagSPD_LP(
  real alpha, real rho_se, real rho_per, real period,
  real L, int M, real log_bessel_cutoff, vector indices
) {
  vector[M] omega = (pi() / 2 / L) * indices;

  // b = 1 / rho_per^2
  real b = inv(square(rho_per));
  
  real scaling_factor = sqrt(2 * pi()) * rho_se;

  vector[M] spd = rep_vector(0.0, M);
    
  real d_coeff = -0.5 * square(rho_se);
  
  // Avoid mean offset by starting at 1
  for(n in 1:64) {
    real log_bessel_scaled = log_modified_bessel_first_kind(n, b) - b;
    
    // If the log-bessel value gets too small, we can safely stop early
    if(log_bessel_scaled < log_bessel_cutoff) break;
    
    vector[M] d_pos = square(omega - 2 * pi() * n / period);
    
    spd += exp(log_bessel_scaled + d_coeff * d_pos);
    
    // For n > 0, we count both sides of the infinite sum
    if(n > 0) {
      vector[M] d_neg = square(omega + 2 * pi() * n / period);
      spd += exp(log_bessel_scaled + d_coeff * d_neg);
    }
  }

  return alpha * sqrt(scaling_factor * spd);
}

vector construct_delay_dist(real delay_mean, real delay_sd, int N) {
  vector[N] delay_pdf;
  
  real sigma = sqrt(log(1 + (delay_sd ^ 2) / (delay_mean ^ 2)));
  real mu = log(delay_mean) - 0.5 * sigma ^ 2;
  
  for(i in 1:N) {
    delay_pdf[i] = lognormal_cdf(i | mu, sigma) - lognormal_cdf(i - 1 | mu, sigma);
  }
  
  delay_pdf = delay_pdf / sum(delay_pdf);
  
  return delay_pdf;
}

matrix qr_weekly() {
  matrix[7, 7] A = diag_matrix(rep_vector(1.0, 7));
  
  A[, 1] = rep_vector(1.0, 7);
  
  matrix[7, 7] Q_full = qr_Q(A);
  
  return Q_full[, 2:7];
}

vector get_effect(matrix effect_raw, int N, int ix_loc) {
  vector[N + 1] effect;
  effect[1] = 0.0;
  effect[2:(N + 1)] = col(effect_raw, ix_loc);
  
  return effect;
}

real sigmoid(real x, real k, real x0, real S) {
  return S * inv_logit(k * (x - x0));
}

int which_max(vector y) {
  return sort_indices_asc(y)[size(y)];
}

real soft_argmax(vector v, real tau, vector indices) {
  int N = num_elements(v);
  vector[N] w = softmax(v / tau);
  return dot_product(w, indices);
}
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
  
  
  vector[21] gi_reversed;

  int<lower=1> M_annual;
  int<lower=1> M_noise; 
  
  real<lower=1.0> c_factor;
}

transformed data {
  // Pre-computed index vectors replacing linspaced_vector
  vector[M_annual] idx_M_annual;
  for (i in 1:M_annual) {
    idx_M_annual[i] = i;
  }

  vector[M_noise] idx_M_noise;
  for (i in 1:M_noise) {
    idx_M_noise[i] = i;
  }

  // Construct the basis functions for HSGP
  vector[N] x_centered = x - mean(x);
  real L = c_factor * max(abs(x_centered));
  
  matrix[N, M_noise] phi_noise = PHI(N, M_noise, L, x_centered, idx_M_noise);
  matrix[N, 2 * M_annual] phi_annual = PHI_periodic(N, M_annual, 2 * pi() / 365.25, x, idx_M_annual);
  
  int K = 21;
  int N_total = N + K;

  // Pre-computed vector sequences
  vector[K] k_seq;
  for (i in 1:K) {
    k_seq[i] = -K + i;
  }
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
  
  // For school holidays and public holidays, type = 1 is no holiday
  // Only fit effects to type > 1
  vector[N_school_holiday_types] school_holiday_Rt_effect_raw;
  vector[N_school_holiday_types] school_holiday_asc_effect_raw;
  
  vector[N_public_holiday_types] public_holiday_Rt_effect_raw;
  vector[N_public_holiday_types] public_holiday_asc_effect_raw;
  
  real<lower=3, upper=21> delay_mean;
  real<lower=1.0> delay_sd;
}

transformed parameters {
  // Calculate length-scale for annual GP in radians
  real ls_annual_rad = ls_annual * (2 * pi() / 365.25);
  
  // Calculate the weighting (SPD) vector for the periodic annual GP 
  // from ls_annual_rad and scaling_annual
  vector[2 * M_annual] diag_SPD_annual = diagSPD_Periodic(scaling_annual, ls_annual_rad, M_annual, idx_M_annual);
  
  // Calculate the weighting (SPD) vector for the Matern 5/2 noise GP
  // from ls_noise and scaling_noise
  vector[M_noise] diag_SPD_noise = diagSPD_Matern52(scaling_noise, ls_noise, L, M_noise, idx_M_noise);
  
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
  
  {
    // The full time-series of infections including M days prior to t = 1
    vector[N_total] inf_full = rep_vector(0.0, N_total);

    // Calculate initial infections as exp growth: I(t) = I_0 * exp(r * t)
    inf_full[1:K] = exp(log_I_0 + r_0 * k_seq);
    
    for(t in 1:N) {
      // Perform the renewal equation step
      real R_t = exp(log_R[t]);
      
      // Use a dot-product here to make the sum over the generation interval
      // a single auto-diff operation
      real new_infections = R_t * dot_product(gi_reversed, segment(inf_full, t, K));
      
      inf_full[t + K] = new_infections;
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
  scaling_noise  ~ exponential(1.0);
  
  scaling_weekly_asc ~ exponential(1.0);
  
  school_holiday_Rt_effect_raw ~ std_normal();
  school_holiday_asc_effect_raw ~ std_normal();
  
  public_holiday_Rt_effect_raw ~ normal(0, 0.25);
  public_holiday_asc_effect_raw ~ normal(0, 0.25);
  
  // Inv-gamma priors derived from:
  // find_inv_gamma_params(7.0, 90.0)
  // find_inv_gamma_params(14.0, 90.0)
  // as in Betancourt approach (Section 3.2.3)
  // https://betanalpha.github.io/assets/case_studies/gaussian_processes.html
  
  ls_noise ~ inv_gamma(3.858, 68.729);
  ls_annual ~ inv_gamma(6.809, 200.117);
  phi ~ gamma(2, 0.1);

  beta_annual ~ std_normal();
  beta_noise  ~ std_normal();
  
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
