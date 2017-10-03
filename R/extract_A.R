#' Extact coefficient matrix A
#'
#' @import tibble
#' @import magrittr
#' @import dplyr
#' @import rstan
#' @import BiocParallel
#' @export
#'
extract_A = function(dm_model_list,
                     protein_names) {

  # some jobs may fail (not sure why yet)
  jobs_ok = which(bpok(dm_model_list))

  # collect from result list
  tb = lapply(jobs_ok,function(seed) {
    fit = dm_model_list[[seed]]
    post = rstan::extract(fit)[["A"]]
    column_condition = post[,,2]
    tibble(protein_name = protein_names,
           coeff = apply(column_condition,2,median),
           min = apply(column_condition,2,quantile,probs = 0.05),
           max = apply(column_condition,2,quantile,probs = 0.95),
           seed = seed)
  }) %>% bind_rows
}
