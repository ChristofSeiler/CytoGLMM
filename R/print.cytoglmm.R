#' Extact and print GLMM fit
#'
#' @aliases print.cytoglmm
#' @method print cytoglmm
#'
#' @export
#'
print.cytoglmm = function(fit) {

  if(class(fit) != "cytoglmm")
    stop("Input needs to be a cytoglmm object computed by cytoglmm function.")

  cat("number of cells per group and condition:")
  cell_count = table(pull(fit$df_samples_subset,fit$group),
                     pull(fit$df_samples_subset,fit$condition))
  print(cell_count)

  cat("\nproteins included in the analysis:\n",fit$protein_names,"\n\n")
  cat("condition compared:",fit$condition,"\n")
  cat("grouping variable:",fit$group,"\n")
  if(!is.null(fit$covariate_names))
    cat("controlled covariates:",fit$covariate_names,"\n")

}
