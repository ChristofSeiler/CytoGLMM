/*
 * Normal-multinomial regression with correlated errors
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
  matrix[d,p] A;
  vector<lower=0>[d] sigma; // donor random effects std
  cholesky_factor_corr[d] L; // cholesky factor of donor random effects corr matrix
  vector[d] z[k]; // donor random effects
}
transformed parameters {
  vector[d] u[k]; // donor random effects
  {
    matrix[d,d] Sigma; // donor random effects cov matrix
    Sigma = diag_pre_multiply(sigma, L);
    for(j in 1:k)
      u[j] = Sigma * z[j];
  }
}
model {
  to_vector(A) ~ normal(0,10);
  L ~ lkj_corr_cholesky(1.0);
  sigma ~ cauchy(0,5);
  for (j in 1:k)
    z[j] ~ normal(0,1);
  for (i in 1:n)
    Y[i] ~ multinomial(softmax(A * X[i] + u[donor[i]]));
}
generated quantities {
  matrix[d,d] Cor;
  Cor = multiply_lower_tri_self_transpose(L);
}
