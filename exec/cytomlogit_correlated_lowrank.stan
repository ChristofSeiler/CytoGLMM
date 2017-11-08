/*
 * Logistic normal multinomial regression with low-rank correlated errors
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
  real gamma[n];
  vector[d] b;
}
model {
  to_vector(A) ~ normal(0,5);
  sigma ~ cauchy(0,5);
  gamma ~ normal(0,1);
  b ~ normal(0,5);
  for (i in 1:n) {
    theta[i] ~ normal(A * X[i] + gamma[i] * b, sigma);
    Y[i] ~ multinomial(softmax(theta[i]));
  }
}
