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
                     cell_n_min = Inf,
                     cell_n_subsample = 0,
                     seed = 0xdada) {

  # some error checks
  if(cell_n_subsample > cell_n_min) stop("cell_n_subsample is larger than cell_n_min")
  if(sum(str_detect(protein_names,"/")) > 0) stop("protein names cannot contain '/'")
  starts_with_number = sapply(protein_names,
         function(x) str_locate(x, "[0-9]")[1] == 1)
  if(sum(starts_with_number,na.rm = TRUE) > 0) stop("protein names cannot start with numbers")
  if(sum(make.names(protein_names) != protein_names) > 0)
    stop("cleanup your protein names (don't use special characters)")
  set.seed(seed)

  # are the samples paired?
  cell_count = table(df_samples_subset$donor,pull(df_samples_subset,condition))
  unpaired = TRUE
  if(sum(apply(cell_count,1,min) > 0) == nrow(cell_count))
    unpaired = FALSE

  # remove donors with low cell count
  include = NULL
  if(cell_n_min < Inf) {
    if(unpaired) {
      exclude = which(rowSums(cell_count) < cell_n_min)
    } else {
      exclude = which(apply(cell_count,1,min) < cell_n_min)
    }
    include = rownames(cell_count)[rownames(cell_count) %nin% names(exclude)]
  } else {
    include = rownames(cell_count)
  }
  df_samples_subset %<>%
    dplyr::filter(donor %in% include) %>%
    droplevels()

  # subsample cells
  if(cell_n_subsample > 0) {
    df_samples_subset %<>%
      group_by_("donor",condition) %>%
      sample_n(size = cell_n_subsample) %>%
      ungroup
  }

  # create formula
  df_samples_subset$donor %<>% as.factor
  donor_dummy = class2ind(df_samples_subset$donor)
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
  fit$cell_n_min = cell_n_min
  fit$cell_n_subsample = cell_n_subsample
  fit$seed = seed
  class(fit) = "cytogroup"
  fit

}
