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
#' @param fit A \code{cytoglm} class
#' @param order Order the markers according to the mangintute of the coefficients
#' @param separate create two separate \code{\link[ggplot2]{ggplot2}} objects
#' @return \code{\link[ggplot2]{ggplot2}} object
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
#' plot(glm_fit)
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
