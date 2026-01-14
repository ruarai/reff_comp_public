


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