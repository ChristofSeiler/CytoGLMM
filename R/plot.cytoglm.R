#' Extact and plot for experimental condition
#'
#' @import ggplot2
#' @import tibble
#' @import magrittr
#' @import dplyr
#' @import cowplot
#' @export
#'
plot.cytoglm = function(fit) {

  if(class(fit) != "cytoglm")
    stop("Input needs to be a cytoglm object computed by cytoglm function.")

  # some jobs may fail (because of computing cluster instabilities)
  if(nrow(fit$tb_coef) == 0)
    stop("no results available")

  xlab_str = fit$df_samples_subset %>%
    pull(fit$condition) %>%
    levels %>%
    paste(collapse = " <-> ")
  tb_A = extract_A(fit$model_fit_list,protein_names = fit$protein_names,column = 2)
  plot_coeff(tb_A,"Differential Expression",xlab_str)

}
