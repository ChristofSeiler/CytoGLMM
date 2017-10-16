/*
 * Low-dimensional covariance multinomial regression
 * Author: Christof Seiler
 */
data {
  int<lower=1> n; // num of cells
  int<lower=1> d; // num of markers
  int<lower=1> p; // num of explanatory variables (including intercept)
  int<lower=1> k; // num donors
  int<lower=1,upper=k> donor[n]; // donor indicator
  int<lower=0> Y[n,d]; // observed cell counts
  vector[p] X[n]; // design matrix
}
parameters {
  //real gamma[n];
  matrix[d,p] A;
  vector<lower=0>[d] L_sigma;
  cholesky_factor_corr[d] L_Omega;
  vector<lower=0>[d] z[k];
  vector[d] theta[n];
}
model {
  matrix[d, d] L_Sigma;
  //gamma ~ normal(0,1);
  to_vector(A) ~ normal(0,1);
  L_Omega ~ lkj_corr_cholesky(4);
  L_sigma ~ cauchy(0, 2.5);
  L_Sigma = diag_pre_multiply(L_sigma, L_Omega);
  for (i in 1:n) {
    theta[i] ~ multi_normal_cholesky(A * X[i] + z[donor[i]], L_Sigma);
    Y[i] ~ multinomial(softmax(theta[i]));
  }
}
generated quantities {
  corr_matrix[d] Omega;
  Omega = multiply_lower_tri_self_transpose(L_Omega);
}
