
matrix qr_weekly() {
  matrix[7, 7] A = diag_matrix(rep_vector(1.0, 7));
  
  A[, 1] = rep_vector(1.0, 7);
  
  matrix[7, 7] Q_full = qr_Q(A);
  
  return Q_full[, 2:7];
}