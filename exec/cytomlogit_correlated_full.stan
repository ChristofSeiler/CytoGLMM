/*
 * Logistic normal multinomial regression with correlated errors
 * Author: Christof Seiler
 */
data {
  int<lower=1> n; // num of cells
  int<lower=1> d; // num of markers
  int<lower=1> p; // num of explanatory variables (including intercept)
  int<lower=0> Y[n,d]; // observed cell counts
  vector[p] X[n]; // design matrix
}
parameters {
  matrix[d,p] A;
  vector<lower=0>[d] sigma;
  vector[d] theta[n];
  cholesky_factor_corr[d] L_Omega; // cholesky factor of donor random effects corr matrix
}
model {
  matrix[d,d] L_Sigma;
  to_vector(A) ~ normal(0,5);
  sigma ~ cauchy(0,5);
  L_Omega ~ lkj_corr_cholesky(4);
  L_Sigma = diag_pre_multiply(sigma, L_Omega);
  for (i in 1:n) {
    theta[i] ~ multi_normal_cholesky(A * X[i], L_Sigma);
    Y[i] ~ multinomial(softmax(theta[i]));
  }
}
generated quantities {
  matrix[d,d] Cor;
  Cor = multiply_lower_tri_self_transpose(L_Omega);
}
