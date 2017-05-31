#' Latent class regression with mixed effects using Stan
#'
#' @import rstan
#' @import magrittr
#' @export
#'
lcr <- function(df_samples,
                protein_names,
                condition,
                num_latent_classes,
                num_bins = 8,
                subsample_size = 10000,
                seed = 1,
                hmc = FALSE,
                cores = 2) {
  # load stan model from file
  file = system.file("exec", "lcr.stan", package = "CytoGLMM")
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
  #subsample_ids = sample(x = nrow(df_samples),
  #                       size = subsample_size,
  #                       replace = FALSE)
  cell_n = round(subsample_size/length(levels(df_samples_binned$donor)))
  cell_n_min = min(table(df_samples_binned$donor))
  if(cell_n_min < cell_n) cell_n = cell_n_min
  subsample_ids = lapply(levels(df_samples_binned$donor),function(donor_id) {
    all_ids = which(df_samples_binned$donor == donor_id)
    sample(x = all_ids,
           size = cell_n,
           replace = FALSE)
  }) %>% unlist
  df_samples_binned = df_samples_binned[subsample_ids,]
  # prepare data for stan
  x = model.matrix(as.formula(paste("~",condition)), data = df_samples_binned)
  donor = as.numeric(df_samples_binned$donor)
  as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
  y = df_samples_binned[,protein_names]
  for(i in 1:ncol(y))
    y[,i] = as.numeric.factor(y[,i])
  y = as.matrix(y)
  # other dimensions
  R = as.integer(num_latent_classes) # num of latent classes
  J = ncol(y)
  K = length(levels(df_samples_binned[,protein_names[1]]))
  D = length(levels(df_samples_binned$donor)) # num of donors
  stan_data = list(N = nrow(y),
                   R = R,
                   J = J,
                   K = K,
                   D = D,
                   alpha = rep(1/K,K),
                   x = x,
                   y = y,
                   donor = donor)
  fit = NULL
  if(!hmc) {
    fit = vb(model,
             data = stan_data,
             seed = seed,
             pars = "log_lik")
  } else {
    fit = sampling(model,
                   data = stan_data,
                   iter = 1000,
                   chains = cores,
                   cores = cores,
                   seed = seed,
                   pars = c("pi","beta"))
  }
  fit
}
