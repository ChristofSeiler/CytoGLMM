#' Extact coefficient matrix A
#'
#' @import tibble
#' @import magrittr
#' @import dplyr
#' @import rstan
#' @export
#'
extract_A = function(dm_model_list,
                     protein_names) {

  # some jobs may fail (because of computing cluster instabilities)
  if(length(dm_model_list) == 0)
    stop("no jobs completed successfully")

  # collect from result list
  tb = lapply(seq_along(dm_model_list),function(run) {
    fit = dm_model_list[[run]]
    tibble(protein_name = protein_names,
           coeff = fit$par$A[,2],
           run = run)
    }) %>% bind_rows

}
