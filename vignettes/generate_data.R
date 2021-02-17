library("MASS")

generate_data = function() {

  # simulation parameters
  n_samples = 16
  n_donors = n_samples/2
  n_cells = 1000
  n_markers = 10
  rho_b = rho_u = 0.1
  sigma_b = sigma_u = 1
  beta_treatment = log(24.7)
  beta_control = log(22.9)
  n_true = 3

  # patient information
  donor = rep(1:n_donors, each = n_cells)
  donor %<>% rep(2)
  condition = c(
    rep("treatment", length(donor)/2),
    rep("control", length(donor)/2)
  )
  df = tibble(donor, condition)

  # generate protein counts
  protein_names = paste0("m", str_pad(1:n_markers, width = 2, pad = "0"))
  rcov = function(rho, sigma) {
    corr = rho^toeplitz(0:(n_markers-1))
    sigma_vec = rep(sigma, n_markers)
    diag(sigma_vec) %*% corr %*% diag(sigma_vec)
  }
  rcov_block = function(rho, sigma) {
    corr = diag(1, nrow = n_markers)
    corr_act = rho^toeplitz(0:(n_true-1))
    corr_notact = rho^toeplitz(0:(n_markers-n_true-1))
    corr[1:n_true,1:n_true] = corr_act
    corr[(n_true+1):n_markers,(n_true+1):n_markers] = corr_notact
    sigma_vec = rep(sigma, n_markers)
    diag(sigma_vec) %*% corr %*% diag(sigma_vec)
  }
  Sigma_b = rcov(rho_b, sigma_b) # cell level variability
  Sigma_u = rcov(rho_u, sigma_u) # donor level variability
  b = mvrnorm(n = nrow(df), mu = rep(0, n_markers), Sigma_b)
  u = mvrnorm(n = n_donors, mu = rep(0, n_markers), Sigma_u)
  u = u[donor, ]
  beta = matrix(beta_control, nrow = nrow(b), ncol = n_markers)
  beta[,1:n_true] = ifelse(condition == "treatment", beta_treatment, beta_control)
  log_lambda = beta + b + u
  lambda = exp(log_lambda)
  y = rpois(length(lambda), lambda)
  dim(y) = dim(lambda)
  colnames(y) = protein_names
  df %<>% bind_cols(as_tibble(y))
  df$condition %<>% factor(levels = c("control", "treatment"))
  df

}
