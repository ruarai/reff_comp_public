

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

real soft_argmax(vector v, real tau) {
  int N = num_elements(v);
  vector[N] indices = linspaced_vector(N, 1, N);
  vector[N] w = softmax(v / tau);
  return dot_product(w, indices);
}