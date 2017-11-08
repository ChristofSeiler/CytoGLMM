#' Normal-multinomial regression model using Stan
#'
#' @import rstan
#' @import magrittr
#' @import stringr
#' @import batchtools
#' @export
#'
cytomlogit = function(df_samples_subset,
                      protein_names,
                      condition,
                      unpaired = TRUE,
                      num_boot = 100,
                      cell_n_max = 1000,
                      seed = 0xdada,
                      partition = "normal") {

  donors = df_samples_subset %>%
    group_by_("donor",condition) %>%
    tally()

  # subsample cells
  # (to speed up computations we subsample at cell level,
  # the results won't change much because the major
  # variability happens at donor level)
  if(min(donors$n) < cell_n_max) stop("some donors have less than cell_n_max number of cells")
  set.seed(seed)
  df_samples_subset %<>%
    group_by(donor,treatment) %>%
    sample_n(size = cell_n_max) %>%
    ungroup

  # prepare cluster script
  cells_total = nrow(donors)*cell_n_max
  # data from previous runs
  time = c(149,152,276,287) # min
  cells = c(84000,86000,168000,172000) # cells
  memory = c(9100,9500,17600,18000) # MB
  expected_walltime = predict(lm(time ~ cells),
                              data.frame(cells = cells_total)) %>% ceiling
  expected_walltime = expected_walltime + 120 # add 2 hours
  expected_mem = predict(lm(memory ~ cells),
                         data.frame(cells = cells_total)) %>% ceiling
  expected_mem = expected_mem + 4000 # add 4GBs
  if(!unpaired) expected_mem = expected_mem + 8000 # paired analysis need more memory
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
  reg$cluster.functions = makeClusterFunctionsSlurm(slurm_settings,
                                                    scheduler.latency = 120,
                                                    fs.latency = 120)
  batchMap(fun = run_vb,
           seed = seq(num_boot),
           more.args = list(df_samples_subset = df_samples_subset,
                            protein_names = protein_names,
                            condition = condition,
                            unpaired = unpaired),
           reg = reg)
  submitJobs(resources = list(ncpus = 1,
                              memory = expected_mem,
                              walltime = expected_walltime,
                              partition = partition,
                              measure.memory = TRUE),
             sleep = 300,
             reg = reg)
  waitForJobs(sleep = 300, reg = reg)
  # if(!"1" %in% findDone()$job.id)
  #   stop("original bootstrap failed")
  model_fit_list = reduceResultsList(missing.val = NULL, reg = reg)

  # return cytomlogit object
  fit = NULL
  fit$model_fit_list = model_fit_list
  fit$df_samples_subset = df_samples_subset
  fit$protein_names = protein_names
  fit$condition = condition
  fit$num_boot = num_boot
  fit$cell_n_max = cell_n_max
  fit$seed = seed
  class(fit) = "cytomlogit"
  fit

}

# cluster function
run_vb = function(seed,
                  df_samples_subset,
                  protein_names,
                  condition,
                  unpaired = TRUE) {

  # need to load it here because this will be run on the cluster
  #library("dplyr")
  #library("magrittr")
  #library("rstan")
  print(seed)
  set.seed(seed)

  # load stan model from file
  stan_file = system.file("exec", "cytomlogit_simple.stan", package = "CytoGLMM")
  model = rstan::stan_model(file = stan_file, model_name = "cytomlogit")

  # cases bootstrap
  # (sample with replacement at donor level)
  # if(seed == 1) { # first seed is reserved for original bootstrap
  #   df_boot = df_samples_subset
  # } else {
  # }
  donor_boot = NULL
  if(unpaired) {
    donor_boot = df_samples_subset %>%
      group_by_("donor",condition) %>%
      tally() %>%
      group_by_(condition) %>%
      sample_frac(replace = TRUE) %>%
      ungroup
  } else {
    donor_boot = df_samples_subset %>%
      group_by(donor) %>%
      tally() %>%
      sample_frac(replace = TRUE)
  }
  df_boot = inner_join(donor_boot,
                       df_samples_subset,
                       by = "donor",
                       suffix = c("",".y")) %>% droplevels

  # prepare data for rstan
  # df_boot %<>% mutate(total = df_boot %>%
  #                       select_at(protein_names) %>%
  #                       rowSums)
  Y = df_boot %>%
    select(protein_names) %>%
    as.matrix
  X = model.matrix(as.formula(paste("~",condition)), data = df_boot)
  #X = model.matrix(as.formula(paste("~",condition,"* total")), data = df_boot)
  donor = df_boot$donor %>% as.factor %>% as.numeric
  n = nrow(Y)
  d = ncol(Y)
  p = ncol(X)
  k = length(table(donor))
  stan_data = list(n = n,
                   d = d,
                   p = p,
                   donor = donor,
                   k = k,
                   Y = Y,
                   X = X)

  # # maximum likelihood estimate
  # init = list(
  #   #gamma = rnorm(n),
  #   A = matrix(rnorm(d*p),nrow = d,ncol = p)
  #   #B = matrix(rnorm(d*p),nrow = d,ncol = p)
  #   )
  # fit = rstan::optimizing(model,
  #                         data = stan_data,
  #                         as_vector = FALSE,
  #                         # init = init,
  #                         verbose = TRUE)
  # # fit$par$theta = NULL
  # fit

  # sample from model using variatonal inference
  fit = rstan::vb(model,
                  iter = 2000,
                  output_samples = 100,
                  pars = c("A","z","sigma"),
                  #pars = c("A","u","sigma","Cor"),
                  #pars = c("A","z","sigma","b"),
                  data = stan_data,
                  seed = 0xdada)
                  #,adapt_engaged = FALSE,
                  #eta = 1)
  A = rstan::extract(fit)[["A"]] %>% apply(c(2,3),median)
  colnames(A) = colnames(X)
  z = rstan::extract(fit)[["z"]] %>% apply(c(2,3),median)
  sigma = rstan::extract(fit)[["sigma"]] %>% apply(2,median)
  #b = rstan::extract(fit)[["b"]] %>% apply(2,median)
  par = NULL
  par$A = A
  par$z = z
  par$sigma = sigma
  #par$b = b
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

  # # Gibbs sampler
  # library("BayesLogit")
  # y_all = Y / rowSums(Y)
  # y = y_all[,-d]
  # fit = mlogit(y = y,
  #              X = X,
  #              n = rowSums(Y),
  #              samp = 1,
  #              burn = 1)

}
