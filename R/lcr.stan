data {
  int<lower=1> D;                             // num donors
  int<lower=1> R;                             // num latent class
  int<lower=1> J;                             // num markers (polytomous variables)
  int<lower=1> K;                             // num of bins
  vector<lower=0>[K] alpha;                   // class-conditional prior
  int<lower=1> N_train;                       // (training) num observations
  int<lower=1,upper=D> donor_train[N_train];  // (training) donor indicator
  vector[2] x_train[N_train];                 // (training) predictor
  int<lower=1,upper=K> y_train[N_train,J];    // (training) response
  int<lower=1> N_test;                        // (test) num observations
  vector[2] x_test[N_test];                   // (test) predictor
  int<lower=1,upper=K> y_test[N_test,J];      // (test) response
}
parameters {
  simplex[K] pi[R,J];       // class-conditional probability
  vector[R] eta[N_train];         // class-membership probability
  vector<lower=0>[R] sigma_e; // class-membership variance
  matrix[R,2] beta;         // regression coefficients
  matrix[R,2] z[D];         // donor random effects
  vector<lower=0>[2] sigma_z;
  //vector<lower=0>[R] sigma_b;
}
transformed parameters {
  simplex[R] theta[N_train];
  for (n in 1:N_train)
    theta[n] = softmax(eta[n]);
}
model {
  // prior
  for (r in 1:R) {
    for (j in 1:J)
      pi[r,j] ~ dirichlet(alpha);
  }
  for (d in 1:D) {
    for (r in 1:R)
      z[d,r] ~ normal(0,sigma_z);
  }
  //sigma_b ~ cauchy(0, 0.5);
  //for (r in 1:R)
  //  beta[r,2] ~ double_exponential(0, sigma_b[r]);
  for (n in 1:N_train)
    eta[n] ~ normal(beta * x_train[n] + z[donor_train[n]] * x_train[n], sigma_e);
  // likelihood
  for (n in 1:N_train) {
    real target_class[R];
    for (r in 1:R) {
      target_class[r] = log(theta[n,r]);
      for (j in 1:J)
        target_class[r] = target_class[r] + log(pi[r,j][y_train[n,j]]);
    }
    target += log_sum_exp(target_class);
  }
}
generated quantities {
  real log_lik[N_test];
  for (n in 1:N_test) {
    matrix[R,2] z_test;
    vector[R] mu;
    vector[R] eta_test;
    vector[R] theta_test;
    real target_class[R];
    for (r in 1:R) {
      z_test[r,1] = normal_rng(0,sigma_z[1]);
      z_test[r,2] = normal_rng(0,sigma_z[2]);
    }
    mu = beta * x_test[n] + z_test * x_test[n];
    for (r in 1:R) {
      eta_test[r] = normal_rng(mu[r], sigma_e[r]);
    }
    theta_test = softmax(eta_test);
    for (r in 1:R) {
      target_class[r] = log(theta_test[r]);
      for (j in 1:J)
        target_class[r] = target_class[r] + log(pi[r,j][y_test[n,j]]);
    }
    log_lik[n] = log_sum_exp(target_class);
  }
}
