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
  vector[p] X1[n]; // design matrix for mean
  vector[k] X2[n]; // design matrix for covariance
}
parameters {
  real gamma[n];
  matrix[d,p] A;
  matrix[d,k] B;
  vector<lower=0>[d] sigma;
  vector<lower=0>[d] z[k];
  vector[d] theta[n];
}
model {
  gamma ~ normal(0,1);
  to_vector(A) ~ normal(0,1);
  to_vector(B) ~ normal(0,1);
  for (i in 1:n) {
    theta[i] ~ normal(A * X1[i] + gamma[i] * B * X2[i] + z[donor[i]], sigma);
    Y[i] ~ multinomial(softmax(theta[i]));
  }
}
