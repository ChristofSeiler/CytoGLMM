/*
 * Normal-multinomial regression
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
  vector<lower=0>[d] sigma;
  vector<lower=0>[d] z[k];
  vector[d] theta[n];
}
model {
  to_vector(A) ~ normal(0,1);
  sigma ~ cauchy(0,5);
  for (i in 1:n) {
    theta[i] ~ normal(A * X[i] + z[donor[i]], sigma);
    Y[i] ~ multinomial(softmax(theta[i]));
  }
}
