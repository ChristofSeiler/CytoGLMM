#' Latent class regression with mixed effects using Stan
#'
#' @import rstan
#' @import rprojroot
#' @export
#'
lcr <- function(df_samples_binned,
                protein_names,
                condition,
                num_latent_classes,
                num_bins = 8,
                subsample = 10000,
                seed = 1) {
  # load stan model from file
  file = find_package_root_file("R", "lcr.stan")
  model = stan_model(file = file,
                     model_name = "lca_model")
  # bining
  count_range = range(df_samples[,protein_names])
  bin_breaks = seq(count_range[1],
                   count_range[2],
                   diff(count_range)/num_bins)
  bin_breaks[1] = -Inf
  bin_breaks[length(bin_breaks)] = Inf
  df_samples_binned = df_samples
  for(name in protein_names)
    df_samples_binned[,name] = cut(df_samples[,name],
                                   breaks = bin_breaks,
                                   labels = 1:num_bins)
  # subsample to speed up computations
  set.seed(seed)
  subids = sample(nrow(df_samples_binned),subsample)
  df_samples_binned = df_samples_binned[subids,]
  # prepare data for stan
  x = model.matrix(as.formula(paste("~",condition)), data = df_samples_binned)
  donor = as.numeric(df_samples_binned$donor)
  as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
  y = df_samples_binned[,protein_names]
  for(i in 1:ncol(y))
    y[,i] = as.numeric.factor(y[,i])
  y = as.matrix(y)
  N = nrow(y)
  R = as.integer(num_latent_classes) # num of latent classes
  J = ncol(y)
  K = length(levels(df_samples_binned[,protein_names[1]]))
  D = length(levels(df_samples_binned$donor)) # num of donors
  stan_data = list(N = N,
                   R = R,
                   J = J,
                   K = K,
                   D = D,
                   alpha = rep(1/K,K),
                   x = x,
                   y = y,
                   donor = donor)
  # variational inference
  fit_vb = vb(model,
              data = stan_data,
              seed = seed,
              pars = c("pi","beta"))
  fit_vb
}
