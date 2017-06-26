#' Latent class regression with mixed effects using Stan
#'
#' @import rstan
#' @import magrittr
#' @export
#'
lcr_cv <- function(df_samples,
                   protein_names,
                   formula,
                   num_latent_classes,
                   num_bins = 8,
                   subsample_size = 10000,
                   seed = 1) {
  # load stan model from file
  file = system.file("exec", "lcr_cv.stan", package = "CytoGLMM")
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
  x = model.matrix(formula, data = df_samples_binned)
  #x = model.matrix(as.formula(paste("~",condition)), data = df_samples_binned)
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
  n_donor = table(donor) %>% length

  compute_elpd = function(donor_test_ind) {
    cat(paste0("num_latent_classes: ",num_latent_classes,", test donor: ",donor_test_ind))
    # split data into training and test part
    donor_train_ind = (1:n_donor)[-donor_test_ind]
    donor_train = donor[donor %in% donor_train_ind] %>% factor(.,labels = 1:8) %>% as.numeric
    x_train = x[donor %in% donor_train_ind,]
    y_train = y[donor %in% donor_train_ind,]
    x_test =  x[donor %in% donor_test_ind,]
    y_test =  y[donor %in% donor_test_ind,]
    stan_data = list(N_train = nrow(y_train),
                     N_test = nrow(y_test),
                     R = R,
                     J = J,
                     K = K,
                     D = length(donor_train_ind),
                     alpha = rep(1/K,K),
                     x_train = x_train,
                     x_test = x_test,
                     y_train = y_train,
                     y_test = y_test,
                     donor_train = donor_train)
    fit = vb(model,
             data = stan_data,
             seed = seed,
             pars = c("log_lik"))
    log_lik = rstan::extract(fit,pars = "log_lik")[[1]]
    log(colMeans(exp(log_lik)))
  }
  elpd = lapply(1:n_donor,function(donor_test_ind)
    compute_elpd(donor_test_ind)) %>% do.call(c,.) %>% sum
  elpd
  # fit = vb(model,
  #          data = stan_data,
  #          seed = seed,
  #          pars = c("pi","beta","log_lik"))
}
