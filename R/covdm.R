#' Covariance Dirichlet multinomial model using Stan
#'
#' @import rstan
#' @import magrittr
#' @import BiocParallel
#' @import BatchJobs
#' @export
#'
covdm = function(df_samples_subset,
                 donors,
                 protein_names,
                 condition,
                 num_boot = 100) {

  # prepare cluster script
  slurm_settings = system.file("exec", "slurm.tmpl", package = "CytoGLMM")
  param = BatchJobsParam(workers = num_boot,
                         resources = list(ntasks=1,ncpus=1,mem=16000,walltime=360),
                         cluster.functions = makeClusterFunctionsSLURM(slurm_settings),
                         log = TRUE,
                         logdir = ".",
                         progressbar = TRUE,
                         cleanup = FALSE,
                         stop.on.error = FALSE,
                         seed = 0xdada)

  # cluster function
  run_vb = function(seed,
                    df_samples_subset,
                    donors,
                    protein_names,
                    condition) {

    # need to load it here because this will be run on the cluster
    library("dplyr")
    library("magrittr")
    library("rstan")
    print(seed)
    set.seed(seed)

    # load stan model from file
    stan_file = system.file("exec", "covdm.stan", package = "CytoGLMM")
    model = rstan::stan_model(file = stan_file, model_name = "covdm_model")

    # prepare bootstrap sample
    boot_donors = sample(donors$donor,replace = TRUE)
    df_boot = lapply(boot_donors,function(boot_donor)
      df_samples_subset %>%
        dplyr::filter(donor == boot_donor) %>%
        sample_n(min(donors$n))
    ) %>% bind_rows
    Y = df_boot %>%
      select(protein_names) %>%
      as.matrix
    X = model.matrix(as.formula(paste("~",condition)), data = df_boot)
    donor = as.numeric(df_samples_subset$donor)
    n = nrow(Y)
    d = ncol(Y)
    p = ncol(X)
    k = length(table(donor))
    stan_data = list(n = n,
                     d = d,
                     p = p,
                     donor = donor,
                     Y = Y,
                     X = X,
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
    #                         init = init,
    #                         verbose = TRUE)

    # sample from model using variatonal inference
    fit = rstan::vb(model,
                    output_samples = 100,
                    pars = c("A","B","sigma"),
                    data = stan_data,
                    seed = 0xdada)

    # # sample using HMC
    # fit = rstan::sampling(model,
    #                       data = stan_data,
    #                       iter = 1000,
    #                       chains = 2,
    #                       cores = 2,
    #                       seed = 0xdada)

    names(fit)
    object.size(fit) %>% print(units = "MB")
    fit
  }

  # run in parallel on cluster
  dm_model_list = bptry({
      bplapply(seq(num_boot),
               run_vb,
               BPPARAM = param,
               df_samples_subset = df_samples_subset,
               donors = donors,
               protein_names = protein_names,
               condition = condition)
    })
  dm_model_list
}
