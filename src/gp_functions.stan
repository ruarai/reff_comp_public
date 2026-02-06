
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
    * @return A vector of spectral densities
    *
    * @ingroup estimates_smoothing
    */
  vector diagSPD_Periodic(real alpha, real rho, int M) {
    real a = inv_square(rho);
    vector[M] indices = linspaced_vector(M, 1, M);
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
    * @return A vector of spectral densities
    *
    * @ingroup estimates_smoothing
    */
  vector diagSPD_Matern52(real alpha, real rho, real L, int M) {
    vector[M] indices = linspaced_vector(M, 1, M);
    real factor = 16 * pow(sqrt(5) / rho, 5);
    vector[M] denom =
      3 * pow((sqrt(5) / rho)^2 + pow((pi() / 2 / L) * indices, 2), 3);
    return alpha * sqrt(factor * inv(denom));
  }
  
  
  matrix PHI(int N, int M, real L, vector x) {
    matrix[N, M] phi = sin(
      diag_post_multiply(
        rep_matrix(pi() / (2 * L) * (x + L), M), linspaced_vector(M, 1, M)
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
    * @return A matrix of basis functions
    *
    * @ingroup estimates_smoothing
    */
  matrix PHI_periodic(int N, int M, real w0, vector x) {
    matrix[N, M] mw0x = diag_post_multiply(
      rep_matrix(w0 * x, M), linspaced_vector(M, 1, M)
    );
    return append_col(cos(mw0x), sin(mw0x));
  }

}
