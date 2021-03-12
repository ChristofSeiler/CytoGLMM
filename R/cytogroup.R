#' Group-specific fixed effects model
#'
#' @import magrittr
#' @importFrom stats binomial
#' @importFrom speedglm speedglm
#' @importFrom caret class2ind
#' @export
#'
#' @param df_samples_subset Data frame or tibble with proteins counts,
#'   cell condition, and group information
#' @param protein_names A vector of column names of protein to use in the analysis
#' @param condition The column name of the condition variable
#' @param group The column name of the group variable
#' @param cell_n_min Remove samples that are below this cell counts threshold
#' @param cell_n_subsample Subsample samples to have this maximum cell count
#'
#' @return A list of class \code{cytoglm} containing
#'   \item{groupfit}{\code{\link[speedglm]{speedglm}} object}
#'   \item{df_samples_subset}{possibly subsampled df_samples_subset table}
#'   \item{protein_names}{input protein names}
#'   \item{condition}{input condition variable}
#'   \item{group}{input group names}
#'   \item{cell_n_min}{input cell_n_min}
#'   \item{cell_n_subsample}{input cell_n_subsample}
#'
#' @examples
#' set.seed(23)
#' df <- generate_data()
#' protein_names <- names(df)[3:12]
#' df <- dplyr::mutate_at(df, protein_names, function(x) asinh(x/5))
#' group_fit <- CytoGLMM::cytogroup(df,
#'                                  protein_names = protein_names,
#'                                  condition = "condition",
#'                                  group = "donor")
#' group_fit
cytogroup <- function(df_samples_subset,
                      protein_names,
                      condition,
                      group = "donor",
                      cell_n_min = Inf,
                      cell_n_subsample = 0) {

  # some error checks
  cyto_check(cell_n_subsample = cell_n_subsample,
             cell_n_min = cell_n_min,
             protein_names = protein_names)

  # are the samples paired?
  unpaired <- is_unpaired(df_samples_subset,
                          condition = condition,
                          group = group)

  # remove donors with low cell count
  df_samples_subset <- remove_samples(df_samples_subset,
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

  # create formula
  df_samples_subset %<>% mutate_at(.vars = group,.funs = as.factor)
  donor_dummy <- class2ind(pull(df_samples_subset,group))
  colnames(donor_dummy) <- paste0("X",colnames(donor_dummy))
  df_samples_subset %<>% bind_cols(as.tibble(donor_dummy))
  pnames <- paste(protein_names, collapse = " + ")
  dnames <- paste(colnames(donor_dummy), collapse = " + ")
  formula_str <- paste0(condition," ~ (",pnames,") * (",dnames,")")

  # logistic regression
  #test <- model.matrix(as.formula(formula_str), df_samples_subset)
  # another option is gpuGlm from package gputools (issues with installation)
  groupfit <- speedglm(formula = formula_str,
                       family = binomial(),
                       data = df_samples_subset)

  # return cytoglmm object
  fit <- NULL
  fit$groupfit <- groupfit
  fit$df_samples_subset <- df_samples_subset
  fit$protein_names <- protein_names
  fit$condition <- condition
  fit$group <- group
  fit$cell_n_min <- cell_n_min
  fit$cell_n_subsample <- cell_n_subsample
  class(fit) <- "cytogroup"
  fit

}
