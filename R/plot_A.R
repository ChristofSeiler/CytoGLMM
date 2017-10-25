#' Extact and plot noise term
#'
#' @import ggplot2
#' @import tibble
#' @import magrittr
#' @import dplyr
#' @import cowplot
#' @export
#'
plot_A = function(df_samples_subset,
                  dm_model_list,
                  protein_names = protein_names) {

  # some jobs may fail (because of computing cluster instabilities)
  if(length(dm_model_list) == 0)
    stop("no jobs completed successfully")

  xlab_str = df_samples_subset %>%
    pull(condition) %>%
    levels %>%
    paste(collapse = " <-> ")
  tb_A = extract_A(dm_model_list,protein_names = protein_names)
  plot_coeff(tb_A,"Differential Expression",xlab_str)

}
