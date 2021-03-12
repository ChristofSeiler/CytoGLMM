#' Logistic mixture regression
#'
#' @import magrittr
#' @import stringr
#' @import flexmix
#' @import BiocParallel
#' @importFrom stats as.formula
#' @export
#'
#' @param df_samples_subset Data frame or tibble with proteins counts,
#'   cell condition, and group information
#' @param protein_names A vector of column names of protein to use in the
#'   analysis
#' @param condition The column name of the condition variable
#' @param group The column name of the group variable
#' @param cell_n_min Remove samples that are below this cell counts threshold
#' @param cell_n_subsample Subsample samples to have this maximum cell count
#' @param ks A vector of cluster sizes
#' @param num_cores Number of computing cores
#'
#' @return A list of class \code{cytoglm} containing
#'   \item{flexmixfits}{list of \code{\link[flexmix]{flexmix}} objects}
#'   \item{df_samples_subset}{possibly subsampled df_samples_subset table}
#'   \item{protein_names}{input protein names}
#'   \item{condition}{input condition variable}
#'   \item{group}{input group names}
#'   \item{cell_n_min}{input cell_n_min}
#'   \item{cell_n_subsample}{input cell_n_subsample}
#'   \item{ks}{input ks}
#'   \item{num_cores}{input num_cores}
#'
#' @examples
#' set.seed(23)
#' df <- generate_data()
#' protein_names <- names(df)[3:12]
#' df <- dplyr::mutate_at(df, protein_names, function(x) asinh(x/5))
#' mix_fit <- CytoGLMM::cytoflexmix(df,
#'                                  protein_names = protein_names,
#'                                  condition = "condition",
#'                                  group = "donor",
#'                                  ks = 2)
#' mix_fit
cytoflexmix <- function(df_samples_subset,
                        protein_names,
                        condition,
                        group = "donor",
                        cell_n_min = Inf,
                        cell_n_subsample = 0,
                        ks = seq_len(10),
                        num_cores = 1) {

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

  # create formulas
  df_samples_subset %<>% mutate_at(.vars = condition,.funs = as.factor)
  df_samples_subset %<>% mutate(
    xtreatment = ifelse(pull(df_samples_subset,condition) ==
                          levels(pull(df_samples_subset,condition))[1],
                        yes = 0,no = 1)
  )
  varying_formula <- paste0("cbind(xtreatment,1-xtreatment) ~ (",
                            paste(protein_names, collapse = " + "),
                            ") | ",group)

  # find best number of cluster
  param <- MulticoreParam(workers = num_cores,
                          tasks = length(ks),
                          progressbar = FALSE)
  flexmixfits <- bplapply(ks,
                          function(k) {
                            stepFlexmix(as.formula(varying_formula),
                                        data = df_samples_subset,
                                        model = FLXMRglm(family = "binomial"),
                                        k = k,
                                        nrep = 5)
                            },
                          BPPARAM = param)

  # return cytoflexmix object
  fit <- NULL
  fit$flexmixfits <- flexmixfits
  fit$df_samples_subset <- df_samples_subset
  fit$protein_names <- protein_names
  fit$condition <- condition
  fit$group <- group
  fit$cell_n_min <- cell_n_min
  fit$cell_n_subsample <- cell_n_subsample
  fit$ks <- ks
  fit$num_cores <- num_cores
  class(fit) <- "cytoflexmix"
  fit

}
