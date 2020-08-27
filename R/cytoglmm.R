#' Estimate random effects model
#'
#' @import magrittr
#' @import stringr
#' @export
#'
cytoglmm = function(df_samples_subset,
                    protein_names,
                    condition,
                    group = "donor",
                    covariate_names = NULL,
                    cell_n_min = Inf,
                    cell_n_subsample = 0,
                    seed = 0xdada,
                    num_cores = 4) {

  set.seed(seed)

  # some error checks
  cyto_check(cell_n_subsample = cell_n_subsample,
             cell_n_min = cell_n_min,
             protein_names = protein_names)
  if(sum(make.names(covariate_names) != covariate_names) > 0)
    stop("cleanup your covariates names (don't use special characters)")

  # are the samples paired?
  unpaired = is_unpaired(df_samples_subset,
                         condition = condition,
                         group = group)

  # remove donors with low cell count
  df_samples_subset = remove_samples(df_samples_subset,
                                     condition = condition,
                                     group = group,
                                     unpaired = unpaired,
                                     cell_n_min = cell_n_min)

  # subsample cells
  if(cell_n_subsample > 0) {
    df_samples_subset %<>%
      group_by_(group,condition) %>%
      sample_n(size = cell_n_subsample) %>%
      ungroup
  }

  glmmfit = glmm_moment(df_samples = df_samples_subset,
                        protein_names = protein_names,
                        response = condition,
                        group = group,
                        covariate_names = covariate_names,
                        num_cores = num_cores)

  # return cytoglmm object
  fit = NULL
  fit$glmmfit = glmmfit
  fit$df_samples_subset = df_samples_subset
  fit$protein_names = protein_names
  fit$condition = condition
  fit$group = group
  fit$covariate_names = covariate_names
  fit$cell_n_min = cell_n_min
  fit$cell_n_subsample = cell_n_subsample
  fit$seed = seed
  fit$num_cores = num_cores
  class(fit) = "cytoglmm"
  fit

}
