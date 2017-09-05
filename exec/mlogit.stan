/*
 * Multi-Logit Regression
 * Author: Christof Seiler
 */
data { 
  int<lower=2> K; // num of bins
  int<lower=0> N; // num of cells
  int<lower=1> D; // num of explanatory variables (including intercept)
  int<lower=1,upper=K> y[N]; // cell expression
  vector[D] x[N]; // explanatory variables
}
transformed data {
  row_vector[D] zeros;
  zeros = rep_row_vector(0, D);
}
parameters {
  matrix[K - 1, D] beta_raw;
}
transformed parameters {
  matrix[K, D] beta;
  beta = append_row(beta_raw, zeros);
}
model {
  for (n in 1:N)
    y[n] ~ categorical(softmax(beta * x[n]));
}
