#' Group-specific fixed effects model
#'
#' @import magrittr
#' @import stringr
#' @import parallel
#' @import flexmix
#' @import cowplot
#' @import caret
#' @import speedglm
#' @export
#'
cytogroup = function(df_samples_subset,
                     protein_names,
                     condition,
                     group = "donor",
                     cell_n_min = Inf,
                     cell_n_subsample = 0,
                     seed = 0xdada) {

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

  # create formula
  df_samples_subset %<>% mutate_at(.vars = group,.funs = as.factor)
  donor_dummy = class2ind(pull(df_samples_subset,group))
  colnames(donor_dummy) = paste0("X",colnames(donor_dummy))
  df_samples_subset %<>% bind_cols(as.tibble(donor_dummy))
  pnames = paste(protein_names, collapse = " + ")
  dnames = paste(colnames(donor_dummy), collapse = " + ")
  formula_str = paste0(condition," ~ (",pnames,") * (",dnames,")")

  # logistic regression
  #test = model.matrix(as.formula(formula_str), df_samples_subset)
  # another option is gpuGlm from package gputools (issues with installation)
  groupfit = speedglm(formula = formula_str,
                      family = binomial(),
                      data = df_samples_subset)

  # return cytoglmm object
  fit = NULL
  fit$groupfit = groupfit
  fit$df_samples_subset = df_samples_subset
  fit$protein_names = protein_names
  fit$condition = condition
  fit$group = group
  fit$cell_n_min = cell_n_min
  fit$cell_n_subsample = cell_n_subsample
  fit$seed = seed
  class(fit) = "cytogroup"
  fit

}
