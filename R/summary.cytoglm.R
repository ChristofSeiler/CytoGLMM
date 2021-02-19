#' Extact and calculate p-values of bootstrap GLM fit
#'
#' @aliases summary.cytoglm
#' @method summary cytoglm
#'
#' @import tibble
#' @import magrittr
#' @import dplyr
#' @export
#'
#' @param fit A \code{cytoglm} class
#' @param method Multiple comparison adjustment method
#' @return \code{\link[tibble]{tibble}} data frame
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
#' summary(glm_fit)
summary.cytoglm = function(fit, method = "BH") {

  if(!is(fit, "cytoglm"))
    stop("Input needs to be a cytoglm object computed by cytoglm function.")

  # calculate p-values from bootstrap distribution
  df_pvalues = fit$tb_coef %>%
    group_by(protein_name) %>%
    summarize(pvalues_unadj = 2*min(mean(coeff < 0), mean(coeff > 0)))
  df_pvalues %<>% mutate(pvalues_unadj = if_else(
    condition = pvalues_unadj == 0,
    true = 2*1/fit$num_boot,
    false = pvalues_unadj))
  df_pvalues %<>% mutate(pvalues_adj = p.adjust(pvalues_unadj,
                                                method = method))
  df_pvalues = df_pvalues[order(df_pvalues$pvalues_unadj),]
  df_pvalues

}
