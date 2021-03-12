#' Fit GLMM with method of moments
#'
#' @import magrittr
#' @import stringr
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
#' @param num_cores Number of computing cores
#'
#' @return A list of class \code{cytoglm} containing
#'   \item{glmmfit}{\code{\link[mbest]{mbest}} object}
#'   \item{df_samples_subset}{possibly subsampled df_samples_subset table}
#'   \item{protein_names}{input protein names}
#'   \item{condition}{input condition variable}
#'   \item{group}{input group names}
#'   \item{covariate_names}{input covariates}
#'   \item{cell_n_min}{input cell_n_min}
#'   \item{cell_n_subsample}{input cell_n_subsample}
#'   \item{num_cores}{input num_cores}
#'
#' @examples
#' set.seed(23)
#' df <- generate_data()
#' protein_names <- names(df)[3:12]
#' df <- dplyr::mutate_at(df, protein_names, function(x) asinh(x/5))
#' glmm_fit <- CytoGLMM::cytoglmm(df,
#'                                protein_names = protein_names,
#'                                condition = "condition",
#'                                group = "donor")
#' glmm_fit
cytoglmm <- function(df_samples_subset,
                     protein_names,
                     condition,
                     group = "donor",
                     covariate_names = NULL,
                     cell_n_min = Inf,
                     cell_n_subsample = 0,
                     num_cores = 1) {

  # suppress printout of warning messages in R package mbest
  logging::setLevel("ERROR", container = "mbest.mhglm")
  logging::setLevel("ERROR", container = "mbest.mhglm.fit")

  # some error checks
  cyto_check(cell_n_subsample = cell_n_subsample,
             cell_n_min = cell_n_min,
             protein_names = protein_names)
  if(sum(make.names(covariate_names) != covariate_names) > 0)
    stop("cleanup your covariates names (don't use special characters)")

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

  glmmfit <- glmm_moment(df_samples = df_samples_subset,
                         protein_names = protein_names,
                         response = condition,
                         group = group,
                         covariate_names = covariate_names,
                         num_cores = num_cores)

  # return cytoglmm object
  fit <- NULL
  fit$glmmfit <- glmmfit
  fit$df_samples_subset <- df_samples_subset
  fit$protein_names <- protein_names
  fit$condition <- condition
  fit$group <- group
  fit$covariate_names <- covariate_names
  fit$cell_n_min <- cell_n_min
  fit$cell_n_subsample <- cell_n_subsample
  fit$num_cores <- num_cores
  class(fit) <- "cytoglmm"
  fit

}
