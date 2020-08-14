#' Logistic regression with cases bootstrap
#'
#' @import magrittr
#' @import stringr
#' @import parallel
#' @import BiocParallel
#' @export
#'
cytoglm = function(df_samples_subset,
                   protein_names,
                   condition,
                   group = "donor",
                   covariate_names = NULL,
                   cell_n_min = Inf,
                   cell_n_subsample = 0,
                   num_boot = 100,
                   seed = 0xdada,
                   num_cores = parallel::detectCores()) {

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

  # formula
  formula_str = paste(condition,"~",
                      paste(c(protein_names, covariate_names),
                            collapse = " + "))

  # bootstrap
  bs = function(seed) {
    set.seed(seed)
    # bootstrap sample
    df_boot = df_samples_subset
    df_boot %<>% group_by(.data[[ group ]], .data[[ condition ]])
    df_boot %<>% slice_sample(prop = 1, replace = TRUE)
    if(unpaired) {
      df_boot %<>% group_by(.data[[ group ]], .data[[ condition ]])
    } else {
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
           run = seed)
  }
  bpparam = MulticoreParam(workers = num_cores)
  tb_coef = bplapply(1:num_boot, bs, BPPARAM = bpparam) %>% bind_rows()

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
  fit$seed = seed
  fit$num_cores = num_cores
  fit$formula_str = formula_str
  class(fit) = "cytoglm"
  fit

}
