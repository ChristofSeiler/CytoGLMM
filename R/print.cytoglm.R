#' Extact and print bootstrap GLM fit
#'
#' @aliases print.cytoglm
#' @method print cytoglm
#' @export
#'
#' @param fit A \code{cytoglm} class
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
#' print(glm_fit)
print.cytoglm = function(fit) {

  if(!is(fit, "cytoglm"))
    stop("Input needs to be a cytoglm object computed by cytoglm function.")

  cat("\n#######################\n")
  if(fit$unpaired)  {
    cat("## unpaired anlaysis ##")
  } else {
    cat("## paired analysis ####")
  }
  cat("\n#######################\n\n")
  cat("number of bootstrap samples:",fit$num_boot,"\n\n")
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
