#' Extact coefficient matrix B
#'
#' @import tibble
#' @import magrittr
#' @import dplyr
#' @import rstan
#' @export
#'
extract_B = function(dm_model_list,ref_run,protein_names) {

  # some jobs may fail (because of computing cluster instabilities)
  if(length(dm_model_list) == 0)
    stop("no jobs completed successfully")

  # measure distance between two B matrices
  distance = function(ref,target) sum((target-ref)^2)

  # extract from list
  fit = dm_model_list[[ref_run]]
  B_ref = fit$par$b

  # collect from result list
  tb = lapply(seq_along(dm_model_list),function(run) {
    fit = dm_model_list[[run]]
    B_target = fit$par$b
    sign_flip = 1
    if(distance(B_ref,B_target) > distance(B_ref,-B_target))
      sign_flip = -1
    B_target = sign_flip*B_target
    tibble(protein_name = protein_names,
           coeff = B_target,
           run = run)
    }) %>% bind_rows

}
