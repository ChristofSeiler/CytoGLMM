#' Latent class regression with mixed effects using Stan.
#'
#' @import rstan
#' @import magrittr
#' @export
#'
dirichlet_multinomial = function(df_samples_subset) {

  stop("under developement")

  # load stan model from file
  #file = system.file("exec", "dirichlet_multinomial.stan", package = "CytoGLMM")
  file = "/Users/cseiler/github_projects/CytoGLMM/exec/dirichlet_multinomial.stan"
  model = stan_model(file = file,
                     model_name = "dm_model")

  # subsample during debugging
  #df_samples_subset %<>% group_by(file_name) %>% sample_n(100)
  df_samples_subset %<>% dplyr::filter(donor == "110007" | donor == "326004")
  df_samples_subset$treatment %<>% relevel(ref = "placebo")
  with(df_samples_subset,table(file_name,treatment))
  y = df_samples_subset %>% pull("CD16")
  x = model.matrix(~ treatment, data = df_samples_subset)
  stan_data = list(y = as.numeric(y),
                   x = x,
                   K = nlevels(y),
                   N = length(y),
                   D = ncol(x))
  fit = rstan::sampling(model,
                        data = stan_data,
                        iter = 1000,
                        chains = 2,
                        cores = 2,
                        seed = 0xdada)
  #beta = rstan::extract(fit)[["beta"]]
  #beta_mean = apply(X = beta, c(2,3), mean)
  #beta_sd = apply(X = beta, c(2,3), sd)
  fit
}
