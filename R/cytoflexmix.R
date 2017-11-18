#' Logistic mixture regression
#'
#' @import magrittr
#' @import stringr
#' @import flexmix
#' @import BiocParallel
#' @export
#'
cytoflexmix = function(df_samples_subset,
                       protein_names,
                       condition,
                       group = "donor",
                       cell_n_min = Inf,
                       cell_n_subsample = 0,
                       ks = 1:10,
                       seed = 0xdada,
                       num_cores = 4) {

  set.seed(seed)

  # some error checks
  cyto_check(cell_n_subsample = cell_n_subsample,
             cell_n_min = cell_n_min,
             protein_names = protein_names)

  # are the samples paired?
  unpaired = is_unpaired(df_samples_subset,
                         condition = condition,
                         group = group)

  # remove donors with low cell count
  df_samples_subset = remove_samples(df_samples_subset,
                                     condition = condition,
                                     group = group,
                                     cell_n_min = cell_n_min)

  # subsample cells
  if(cell_n_subsample > 0) {
    df_samples_subset %<>%
      group_by_(group,condition) %>%
      sample_n(size = cell_n_subsample) %>%
      ungroup
  }

  # create formulas
  df_samples_subset %<>% mutate_at(.vars = condition,.funs = as.factor)
  df_samples_subset %<>% mutate(
    xtreatment = ifelse(pull(df_samples_subset,condition) ==
                          levels(pull(df_samples_subset,condition))[1],yes = 0,no = 1)
  )
  #varying_formula = paste0("cbind(xtreatment,1-xtreatment) ~ 1 | donor")
  #fixed_formula = paste("~",paste(protein_names, collapse = " + "))
  varying_formula = paste0("cbind(xtreatment,1-xtreatment) ~ (",
                           paste(protein_names, collapse = " + "),
                           ") | ",group)

  # find best number of cluster
  param = MulticoreParam(workers = num_cores,
                         tasks = length(ks),
                         progressbar = TRUE,
                         RNGseed = seed)
  flexmixfits = bplapply(ks,
                         function(k) {
                           stepFlexmix(as.formula(varying_formula),
                                       data = df_samples_subset,
                                       model = FLXMRglm(family = "binomial"),
                                       #model = FLXMRglmfix(family = "binomial",
                                       #                    fixed = as.formula(fixed_formula)),
                                       k = k,
                                       nrep = 5)
                           },
                         BPPARAM = param)

  # return cytoflexmix object
  fit = NULL
  fit$flexmixfits = flexmixfits
  fit$df_samples_subset = df_samples_subset
  fit$protein_names = protein_names
  fit$condition = condition
  fit$group = group
  fit$cell_n_min = cell_n_min
  fit$cell_n_subsample = cell_n_subsample
  fit$seed = seed
  fit$ks = ks
  fit$num_cores = num_cores
  class(fit) = "cytoflexmix"
  fit

}
