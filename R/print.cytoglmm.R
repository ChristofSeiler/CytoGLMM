#' Extact and print GLMM fit
#'
#' @aliases print.cytoglmm
#' @method print cytoglmm
#' @export
#'
#' @param x A \code{cytoglmm} class
#' @param ... Other parameters
#'
#' @examples
#' set.seed(23)
#' df = generate_data()
#' protein_names = names(df)[3:12]
#' df = dplyr::mutate_at(df, protein_names, function(x) asinh(x/5))
#' glmm_fit = CytoGLMM::cytoglmm(df,
#'                               protein_names = protein_names,
#'                               condition = "condition",
#'                               group = "donor")
#' print(glmm_fit)
print.cytoglmm = function(x, ...) {

  if(!is(x, "cytoglmm"))
    stop("Input needs to be a cytoglmm object computed by cytoglmm function.")

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
