#' Covariance Dirichlet multinomial model using Stan
#'
#' @import rstan
#' @import magrittr
#' @import batchtools
#' @export
#'
covdm = function(df_samples_subset,
                 donors,
                 protein_names,
                 condition,
                 num_boot = 100,
                 cell_n_max = 1000,
                 seed = 0xdada) {

  # subsample cells
  # (to speed up computations we subsample at cell level,
  # the results won't change much because the major
  # variability happens at donor level)
  set.seed(seed)
  df_samples_subset = apply(donors,1,function(donor_cells) {
    df_donor = df_samples_subset %>%
      dplyr::filter(donor == donor_cells["donor"])
    if(nrow(df_donor) > cell_n_max)
      df_donor %<>% sample_n(size = cell_n_max)
    df_donor
  }) %>% bind_rows

  # prepare cluster script
  cells_total = nrow(donors)*cell_n_max
  # data from previous runs
  time = c(149,152,276,287) # min
  cells = c(84000,86000,168000,172000) # cells
  memory = c(9100,9500,17600,18000) # MB
  expected_walltime = predict(lm(time ~ cells),
                              data.frame(cells = cells_total)) %>% ceiling
  expected_walltime = expected_walltime + 120 # add two hours
  expected_mem = predict(lm(memory ~ cells),
                         data.frame(cells = cells_total)) %>% ceiling
  expected_mem = expected_mem + 4000 # add 4GBs
  cat("requested walltime:",expected_walltime,"min\n")
  cat("requested mem:",expected_mem,"MB")

  # # OLD: use BiocParallel (which uses BatchJobs internally)
  # slurm_settings = system.file("exec", "slurm.tmpl", package = "CytoGLMM")
  # param = BatchJobsParam(workers = num_boot,
  #                        resources = list(ntasks = 1,
  #                                         ncpus = 1,
  #                                         mem = expected_mem,
  #                                         walltime = expected_walltime),
  #                        cluster.functions = makeClusterFunctionsSLURM(slurm_settings),
  #                        log = TRUE,
  #                        logdir = ".",
  #                        progressbar = TRUE,
  #                        cleanup = FALSE,
  #                        stop.on.error = FALSE,
  #                        seed = 0xdada)
  #
  # # run in parallel on cluster
  # dm_model_list = bptry({
  #     bplapply(seq(num_boot),
  #              run_vb,
  #              BPPARAM = param,
  #              df_samples_subset = df_samples_subset,
  #              donors = donors,
  #              protein_names = protein_names,
  #              condition = condition)
  #   })

  # replace BiocParallel,BatchJobs with batchtools
  # NEW: use batchtools
  slurm_settings = system.file("exec", "slurm_batchtools.tmpl", package = "CytoGLMM")
  current_time = Sys.time() %>%
    str_replace_all(":","") %>%
    str_replace_all("-| ","_")
  reg = makeRegistry(file.dir = paste0("registry_",current_time),
                     packages = c("dplyr","magrittr","rstan"))
  reg$cluster.functions = makeClusterFunctionsSlurm(slurm_settings)
  saveRegistry(reg)
  batchMap(fun = run_vb,
           seed = seq(num_boot),
           more.args = list(df_samples_subset = df_samples_subset,
                            donors = donors,
                            protein_names = protein_names,
                            condition = condition),
           reg = reg)
  submitJobs(resources = list(ncpus = 1,
                              memory = expected_mem,
                              walltime = expected_walltime,
                              partition = "hns,normal",
                              measure.memory = TRUE))
  waitForJobs(reg = reg,sleep = 300)
  getStatus()
  findErrors()
  getErrorMessages()
  dm_model_list = reduceResultsList(missing.val = NULL, reg = reg)

  dm_model_list
}

# cluster function
run_vb = function(seed,
                  df_samples_subset,
                  donors,
                  protein_names,
                  condition) {

  # need to load it here because this will be run on the cluster
  #library("dplyr")
  #library("magrittr")
  #library("rstan")
  print(seed)
  set.seed(seed)

  # load stan model from file
  stan_file = system.file("exec", "covdm5.stan", package = "CytoGLMM")
  model = rstan::stan_model(file = stan_file, model_name = "covdm_model")

  # cases bootstrap
  # (sample with replacement at donor level)
  # cases bootstrap
  # (sample with replacement at donor level)
  boot_donors = NULL
  if(seed == 1) { # first seed is reserved for original bootstrap
    boot_donors =  donors$donor
  } else {
    boot_donors = sample(donors$donor,replace = TRUE)
  }
  df_boot = lapply(boot_donors,function(boot_donor)
    df_samples_subset %>%
      dplyr::filter(donor == boot_donor)
  ) %>% bind_rows %>% droplevels

  # prepare data for rstan
  Y = df_boot %>%
    select(protein_names) %>%
    as.matrix
  X1 = model.matrix(as.formula(paste("~",condition)), data = df_boot)
  X2 = model.matrix(~ donor, data = df_boot)
  donor = df_boot$donor %>% as.factor %>% as.numeric
  n = nrow(Y)
  d = ncol(Y)
  p = ncol(X1)
  k = length(table(donor))
  stan_data = list(n = n,
                   d = d,
                   p = p,
                   donor = donor,
                   Y = Y,
                   X1 = X1,
                   X2 = X2,
                   k = k)

  # # maximum likelihood estimate
  # init = list(
  #   gamma = rnorm(n),
  #   A = matrix(rnorm(d*p),nrow = d,ncol = p),
  #   B = matrix(rnorm(d*p),nrow = d,ncol = p)
  #   )
  # fit = rstan::optimizing(model,
  #                         data = stan_data,
  #                         as_vector = FALSE,
  #                         # init = init,
  #                         verbose = TRUE)
  # fit$par$theta = NULL
  # fit

  # # sample from model using variatonal inference
  fit = rstan::vb(model,
                  iter = 2000,
                  output_samples = 100,
                  #pars = c("A","sigma","z"),
                  pars = c("A","B","sigma","z"),
                  #pars = c("A","z","L_sigma","Omega"),
                  data = stan_data,
                  seed = 0xdada)
  A = rstan::extract(fit)[["A"]] %>% apply(c(2,3),median)
  B = rstan::extract(fit)[["B"]] %>% apply(c(2,3),median)
  z = rstan::extract(fit)[["z"]] %>% apply(c(2,3),median)
  sigma = rstan::extract(fit)[["sigma"]] %>% apply(2,median)
  #L_sigma = rstan::extract(fit)[["L_sigma"]] %>% apply(2,median)
  #Omega = rstan::extract(fit)[["Omega"]] %>% apply(c(2,3),median)
  par = NULL
  par$A = A
  par$B = B
  par$z = z
  par$sigma = sigma
  #par$L_sigma = L_sigma
  #par$Omega = Omega
  res = NULL
  res$par = par
  res

  # # sample using HMC
  # fit = rstan::sampling(model,
  #                       data = stan_data,
  #                       iter = 1000,
  #                       chains = 1,
  #                       cores = 1,
  #                       seed = 0xdada)
  # A = rstan::extract(fit)[["A"]] %>% apply(c(2,3),median)
  # B = rstan::extract(fit)[["B"]] %>% apply(c(2,3),median)
  # sigma = rstan::extract(fit)[["sigma"]] %>% apply(c(2,3),median)
  # par = NULL
  # par$A = A
  # par$B = B
  # par$sigma = sigma
  # res = NULL
  # res$par = par
  # res

  # Gibbs sampler
  # library("BayesLogit")
  # y_all = Y / rowSums(Y)
  # y = y_all[,-d]
  # fit = mlogit(y = y,
  #              X = X,
  #              n = rowSums(Y),
  #              samp = 1,
  #              burn = 1)

}
