#' Plot bootstraped coefficients
#'
#' @aliases plot.cytoglm
#' @method plot cytoglm
#'
#' @import ggplot2
#' @import tibble
#' @import magrittr
#' @import dplyr
#' @import cowplot
#' @export
#'
plot.cytoglm = function(fit, order = FALSE, separate = FALSE) {

  if(!is(fit, "cytoglm"))
    stop("Input needs to be a cytoglm object computed by cytoglm function.")

  # some jobs may fail (because of computing cluster instabilities)
  if(nrow(fit$tb_coef) == 0)
    stop("no results available")

  xlab_str = fit$df_samples_subset %>%
    pull(fit$condition) %>%
    levels %>%
    paste(collapse = " <-> ")

  plot_coeff(tb = fit$tb_coef,
             title_str = "Bootstraps",
             title_str_right = "cytoglm",
             xlab_str = xlab_str,
             order = order,
             separate = separate)

}
