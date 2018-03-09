#' Extact and plot noise term
#'
#' @import tibble
#' @import magrittr
#' @import dplyr
#' @export
#'
summary.cytoglmm = function(fit) {

  if(class(fit) != "cytoglmm")
    stop("Input needs to be a cytoglmm object computed by cytoglmm function.")

  pvalues_unadj = summary(fit$glmmfit)$coefficients[-1,4]
  pvalues_adj = p.adjust(pvalues_unadj,method = "BH")
  df_pvalues = tibble(protein_names = names(pvalues_unadj),
                      pvalues_unadj,
                      pvalues_adj)
  df_pvalues$protein_names = as.character(df_pvalues$protein_names)
  df_pvalues = df_pvalues[order(df_pvalues$pvalues_unadj),]
  rownames(df_pvalues) = NULL
  df_pvalues

}
