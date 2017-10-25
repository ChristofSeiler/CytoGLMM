#' Extact and plot noise term
#'
#' @import ggplot2
#' @import tibble
#' @import magrittr
#' @import dplyr
#' @import cowplot
#' @export
#'
plot_sigma = function(fit) {

  if(class(fit) != "cytomlogit")
    stop("Input needs to be a cytomlogit object computed by cytomlogit function.")

  # some jobs may fail (because of computing cluster instabilities)
  if(length(fit$model_fit_list) == 0)
    stop("no jobs completed successfully")

  # extract coefficients
  tb = lapply(seq_along(fit$model_fit_list),function(run) {
    model_fit = fit$model_fit_list[[run]]
    tibble(protein_name = fit$protein_names,
           coeff = model_fit$par$sigma,
           run = run)
  }) %>% bind_rows

  plot_coeff(tb,
             title_str = "Protein-Specific Dispersion",
             xlab_str = "sigma",hline = 1)

}
