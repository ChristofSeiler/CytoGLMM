#' Extact and calculate p-values of GLMM fit
#'
#' @aliases summary.cytoglmm
#' @method summary cytoglmm
#'
#' @import tibble
#' @import magrittr
#' @import dplyr
#' @export
#'
#' @param fit A \code{cytoglmm} class
#' @param method Multiple comparison adjustment method
#' @param ... Other parameters
#' @return \code{\link[tibble]{tibble}} data frame
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
#' summary(glmm_fit)
summary.cytoglmm = function(fit, method = "BH", ...) {

  if(!is(fit, "cytoglmm"))
    stop("Input needs to be a cytoglmm object computed by cytoglmm function.")

  pvalues_unadj = summary(fit$glmmfit)$coefficients[-1,4]
  pvalues_adj = p.adjust(pvalues_unadj,method = method)
  df_pvalues = tibble(protein_name = names(pvalues_unadj),
                      pvalues_unadj,
                      pvalues_adj)
  df_pvalues$protein_name = as.character(df_pvalues$protein_name)
  df_pvalues = df_pvalues[order(df_pvalues$pvalues_unadj),]
  rownames(df_pvalues) = NULL
  df_pvalues

}
