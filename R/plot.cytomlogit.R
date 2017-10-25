#' Extact and plot noise term
#'
#' @import ggplot2
#' @import tibble
#' @import magrittr
#' @import dplyr
#' @import cowplot
#' @export
#'
plot.cytomlogit = function(fit) {

  if(class(fit) != "cytomlogit")
    stop("Input needs to be a cytomlogit object computed by cytomlogit function.")

  # some jobs may fail (because of computing cluster instabilities)
  if(length(fit$model_fit_list) == 0)
    stop("no jobs completed successfully")

  xlab_str = fit$df_samples_subset %>%
    pull(fit$condition) %>%
    levels %>%
    paste(collapse = " <-> ")
  tb_A = extract_A(fit$model_fit_list,protein_names = fit$protein_names)
  plot_coeff(tb_A,"Differential Expression",xlab_str)

}
