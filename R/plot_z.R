#' Extact donor random effect z
#'
#' @import ggplot2
#' @import tibble
#' @import magrittr
#' @import dplyr
#' @import cowplot
#' @export
#'
plot_z = function(dm_model_list,
                  protein_names = protein_names) {

  # some jobs may fail (because of computing cluster instabilities)
  if(length(dm_model_list) == 0)
    stop("no jobs completed successfully")

  # extract coefficients
  tb = lapply(seq_along(dm_model_list),function(run) {
    fit = dm_model_list[[run]]
    z = fit$par$z
    tibble(protein_name = protein_names,
           coeff = apply(z,2,median),
           run = run)
  }) %>% bind_rows

  plot_coeff(tb,
             title_str = "Donor-Specific Expression",
             xlab_str = "z",hline = 0)

}
