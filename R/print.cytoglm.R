#' Extact and print bootstrap GLM fit
#'
#' @importFrom methods is
#' @aliases print.cytoglm
#' @method print cytoglm
#' @export
#'
#' @param x A \code{cytoglm} class
#' @param ... Other parameters
#' @return NULL.
#'
#' @examples
#' set.seed(23)
#' df <- generate_data()
#' protein_names <- names(df)[3:12]
#' df <- dplyr::mutate_at(df, protein_names, function(x) asinh(x/5))
#' glm_fit <- CytoGLMM::cytoglm(df,
#'                              protein_names = protein_names,
#'                              condition = "condition",
#'                              group = "donor",
#'                              num_boot = 10) # in practice >=1000
#' print(glm_fit)
print.cytoglm <- function(x, ...) {

  if(!is(x, "cytoglm"))
    stop("Input needs to be a cytoglm object computed by cytoglm function.")

  cat("\n#######################\n")
  if(x$unpaired)  {
    cat("## unpaired anlaysis ##")
  } else {
    cat("## paired analysis ####")
  }
  cat("\n#######################\n\n")
  cat("number of bootstrap samples:",x$num_boot,"\n\n")
  cat("number of cells per group and condition:")
  cell_count = table(pull(x$df_samples_subset,x$group),
                     pull(x$df_samples_subset,x$condition))
  print(cell_count)

  cat("\nproteins included in the analysis:\n",x$protein_names,"\n\n")
  cat("condition compared:",x$condition,"\n")
  cat("grouping variable:",x$group,"\n")
  if(!is.null(x$covariate_names))
    cat("controlled covariates:",x$covariate_names,"\n")

}
