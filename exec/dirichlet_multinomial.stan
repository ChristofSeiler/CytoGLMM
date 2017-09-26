/*
 * Multi-Logit Regression
 * Author: Christof Seiler
 */
data { 
  //int<lower=2> K; // num of bins
  int<lower=1> n; // num of cells
  int<lower=1> d; // num of markers
  int<lower=1> p; // num of explanatory variables (including intercept)
  int<lower=0> Y[n,d]; // observed cell counts
  matrix[n,p] X; // design matrix
  //vector[d] alpha;
}
transformed data {
  vector[d] alpha;
  //row_vector[d] alpha;
  alpha = rep_vector(1, d);
}
parameters {
  //matrix[K - 1, D] beta_raw;
  simplex[d] theta;
}
transformed parameters {
  //matrix[K, D] beta;
  //beta = append_row(beta_raw, zeros);
}
model {
  theta ~ dirichlet(alpha);
  for (i in 1:n)
    Y[i] ~ multinomial(theta);
  //Y[n] ~ multinomial(softmax(beta * x[n]));
}
