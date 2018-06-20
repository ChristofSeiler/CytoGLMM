#' Extact and calculate p-values of bootstrap GLM fit
#'
#' @import tibble
#' @import magrittr
#' @import dplyr
#' @export
#'
summary.cytoglm = function(fit) {

  if(class(fit) != "cytoglm")
    stop("Input needs to be a cytoglm object computed by cytoglm function.")

  # calculate p-values from bootstrap distribution
  df_pvalues = fit$tb_coef %>%
    group_by(protein_name) %>%
    summarize(pvalues_unadj = 2*min(mean(coeff < 0), mean(coeff > 0)))
  df_pvalues %<>% mutate(pvalues_unadj = if_else(
    condition = pvalues_unadj == 0,
    true = 1/fit$num_boot,
    false = pvalues_unadj))
  df_pvalues %<>% mutate(pvalues_adj = p.adjust(pvalues_unadj,
                                                method = "BH"))
  df_pvalues = df_pvalues[order(df_pvalues$pvalues_unadj),]
  df_pvalues

}
