#' Extact and plot noise term
#'
#' @import ggplot2
#' @import tibble
#' @import magrittr
#' @import dplyr
#' @import cowplot
#' @export
#'
plot_b = function(fit) {

  if(class(fit) != "cytomlogit")
    stop("Input needs to be a cytomlogit object computed by cytomlogit function.")

  # some jobs may fail (because of computing cluster instabilities)
  if(length(fit$model_fit_list) == 0)
    stop("no jobs completed successfully")

  tb_B = extract_B(fit$model_fit_list,
                   ref_run = 1,
                   protein_names = fit$protein_names)
  plot_coeff(tb_B,"Residual Co-Expression","b")

}
