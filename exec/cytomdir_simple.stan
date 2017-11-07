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
  vector[d] z[k];
  simplex[d] theta[n];
}
model {
  vector[d] alpha;
  to_vector(A) ~ normal(0,5);
  for (j in 1:k)
    z[j] ~ normal(0,5);
  for (i in 1:n) {
    alpha = exp(A * X[i] + z[donor[i]]);
    theta[i] ~ dirichlet(alpha);
    Y[i] ~ multinomial(theta[i]);
  }
}
