#' Fit GLM with bootstrap resampling
#'
#' @import magrittr
#' @import stringr
#' @import BiocParallel
#' @export
#'
#' @param df_samples_subset Data frame or tibble with proteins counts,
#'   cell condition, and group information
#' @param protein_names A vector of column names of protein to use in the analysis
#' @param condition The column name of the condition variable
#' @param group The column name of the group variable
#' @param covariate_names The column names of covariates
#' @param cell_n_min Remove samples that are below this cell counts threshold
#' @param cell_n_subsample Subsample samples to have this maximum cell count
#' @param num_boot Number of bootstrap samples
#' @param num_cores Number of computing cores
#'
#' @return A list of class \code{cytoglm} containing
#'   \item{tb_coef}{coefficent table}
#'   \item{df_samples_subset}{possibly subsampled df_samples_subset table}
#'   \item{protein_names}{input protein names}
#'   \item{condition}{input condition variable}
#'   \item{group}{input group names}
#'   \item{covariate_names}{input covariates}
#'   \item{cell_n_min}{input cell_n_min}
#'   \item{cell_n_subsample}{input cell_n_subsample}
#'   \item{unpaired}{true if unpaired samples were provided as input}
#'   \item{num_boot}{input num_boot}
#'   \item{num_cores}{input num_cores}
#'   \item{formula_str}{formula use in the regression model}
#'
#' @examples
#' set.seed(23)
#' df = generate_data()
#' protein_names = names(df)[3:12]
#' df = dplyr::mutate_at(df, protein_names, function(x) asinh(x/5))
#' glm_fit = CytoGLMM::cytoglm(df,
#'                             protein_names = protein_names,
#'                             condition = "condition",
#'                             group = "donor",
#'                             num_boot = 10) # just for docs, in practice >=1000
#' glm_fit
cytoglm = function(df_samples_subset,
                   protein_names,
                   condition,
                   group = "donor",
                   covariate_names = NULL,
                   cell_n_min = Inf,
                   cell_n_subsample = 0,
                   num_boot = 100,
                   num_cores = 1) {

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

  # formula
  formula_str = paste(condition,"~",
                      paste(c(protein_names, covariate_names),
                            collapse = " + "))

  # bootstrap
  bs = function(i) {
    # bootstrap sample
    df_boot = df_samples_subset
    df_boot %<>% group_by(.data[[ group ]], .data[[ condition ]])
    df_boot %<>% slice_sample(prop = 1, replace = TRUE)
    if(!unpaired) {
      df_boot %<>% group_by(.data[[ group ]])
    }
    df_boot %<>%
      group_split() %>%
      sample(replace = TRUE) %>%
      bind_rows()

    # logistic regression
    fit_glm = glm(formula = formula_str,
                  family = binomial(),
                  data = df_boot)
    tibble(protein_name = protein_names,
           coeff = fit_glm$coefficients[protein_names],
           run = i)
  }
  bpparam = MulticoreParam(workers = num_cores)
  tb_coef = bplapply(seq_len(num_boot), bs, BPPARAM = bpparam) %>% bind_rows()

  # return cytoglm object
  fit = NULL
  fit$tb_coef = tb_coef
  fit$df_samples_subset = df_samples_subset
  fit$protein_names = protein_names
  fit$condition = condition
  fit$group = group
  fit$covariate_names = covariate_names
  fit$cell_n_min = cell_n_min
  fit$cell_n_subsample = cell_n_subsample
  fit$unpaired = unpaired
  fit$num_boot = num_boot
  fit$num_cores = num_cores
  fit$formula_str = formula_str
  class(fit) = "cytoglm"
  fit

}
