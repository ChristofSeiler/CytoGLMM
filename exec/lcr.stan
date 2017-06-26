data {
  int<lower=1> D;                 // num donors
  int<lower=1> R;                 // num latent class
  int<lower=1> J;                 // num markers (polytomous variables)
  int<lower=1> P;                 // num predictors
  int<lower=1> K;                 // num of bins
  vector<lower=0>[K] alpha;       // class-conditional prior
  int<lower=1> N;                 // num observations
  int<lower=1,upper=D> donor[N];  // donor indicator
  vector[P] x[N];                 // predictor
  int<lower=1,upper=K> y[N,J];    // response
}
parameters {
  simplex[K] pi[R,J];         // class-conditional probability
  vector[R] eta[N];     // class-membership probability
  vector<lower=0>[R] sigma_e; // class-membership variance
  matrix[R,P] beta;           // regression coefficients
  //matrix[R,2] z[D];           // donor random effects
  //vector<lower=0>[2] sigma_z;
  //vector<lower=0>[R] sigma_b;
}
transformed parameters {
  simplex[R] theta[N];
  for (n in 1:N)
    theta[n] = softmax(eta[n]);
}
model {
  // prior
  for (r in 1:R) {
    for (j in 1:J)
      pi[r,j] ~ dirichlet(alpha);
  }
  //for (d in 1:D) {
  //  for (r in 1:R)
  //    z[d,r] ~ normal(0,sigma_z);
  //}
  //sigma_b ~ cauchy(0, 0.5);
  //for (r in 1:R)
  //  beta[r,2] ~ double_exponential(0, sigma_b[r]);
  for (n in 1:N)
    eta[n] ~ normal(beta * x[n], sigma_e);
    //eta[n] ~ normal(beta * x[n] + z[donor[n]] * x[n], sigma_e);
  // likelihood
  for (n in 1:N) {
    real target_class[R];
    for (r in 1:R) {
      target_class[r] = log(theta[n,r]);
      for (j in 1:J)
        target_class[r] = target_class[r] + log(pi[r,j][y[n,j]]);
    }
    target += log_sum_exp(target_class);
  }
}
generated quantities {
  real log_lik[N];
  for (n in 1:N) {
    real target_class[R];
    for (r in 1:R) {
      target_class[r] = log(theta[n,r]);
      for (j in 1:J)
        target_class[r] = target_class[r] + log(pi[r,j][y[n,j]]);
    }
    log_lik[n] = log_sum_exp(target_class);
  }
}
