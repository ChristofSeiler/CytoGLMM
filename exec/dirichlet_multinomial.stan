/*
 * Dirichlet multinomial regression
 * Author: Christof Seiler
 */
data { 
  //int<lower=2> K; // num of bins
  int<lower=1> n; // num of cells
  int<lower=1> d; // num of markers
  int<lower=1> p; // num of explanatory variables (including intercept)
  //int<lower=0>[d] Y[n]; // observed cell counts
  int<lower=0> Y[n,d]; // observed cell counts
  //matrix[n,p] X; // design matrix
  vector[p] X[n]; // design matrix
  //vector[d] alpha;
}
//transformed data {
  //vector[d] alpha;
  //alpha = rep_vector(1, d);
  //row_vector[d] alpha;
//}
parameters {
  //simplex[d] theta;
  //matrix[K - 1, D] beta_raw;
  matrix[d,p] beta;
}
//transformed parameters {
  //matrix[K, D] beta;
  //beta = append_row(beta_raw, zeros);
//}
model {
  vector[d] alpha;
  for (j in 1:d)
    beta[j, ] ~ normal(0, 5);
  //theta ~ dirichlet(alpha);
  //for (i in 1:n)
  //  Y[i] ~ multinomial(theta);
  for (i in 1:n) {
    alpha = softmax(beta * X[i]);
    Y[i] ~ multinomial(alpha);    
  }
}
