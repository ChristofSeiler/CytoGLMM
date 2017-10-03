#' Extact coefficient matrix B
#'
#' @import tibble
#' @import magrittr
#' @import dplyr
#' @import rstan
#' @import BiocParallel
#' @export
#'
extract_B = function(dm_model_list,job_id,protein_names) {

  # collect tables from result list
  jobs_ok = which(bpok(dm_model_list))

  # extract from stan object
  fit = dm_model_list[[job_id]]
  post = rstan::extract(fit)[["B"]]
  B_ref = apply(post,c(2,3),median)

  # align other estimates
  lapply(jobs_ok,function(seed) {
    # extract from stan object
    fit = dm_model_list[[seed]]
    post = rstan::extract(fit)[["B"]]
    B_target = apply(post,c(2,3),median)
    distance = function(ref,target) sum((target-ref)^2)
    sign_flip = 1
    if(distance(B_ref,B_target) > distance(B_ref,-B_target))
      sign_flip = -1
    column_condition = sign_flip*post[,,2]
    tibble(protein_name = protein_names,
           coeff = apply(column_condition,2,median),
           min = apply(column_condition,2,quantile,probs = 0.05),
           max = apply(column_condition,2,quantile,probs = 0.95),
           seed = seed)
  }) %>% bind_rows
}
