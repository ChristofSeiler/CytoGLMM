#' Extact and plot noise term
#'
#' @import ggplot2
#' @import tibble
#' @import magrittr
#' @import dplyr
#' @import cowplot
#' @export
#'
plot_b = function(dm_model_list,
                  protein_names = protein_names) {

  # some jobs may fail (because of computing cluster instabilities)
  if(length(dm_model_list) == 0)
    stop("no jobs completed successfully")

  tb_B = extract_B(dm_model_list,
                   ref_run = 1,
                   protein_names = protein_names)
  plot_coeff(tb_B,"Residual Co-Expression","b")

}
