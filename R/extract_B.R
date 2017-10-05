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
  if(length(jobs_ok) == 0) stop("no jobs completed successfully")
  
  # measure distance between two B matrices
  distance = function(ref,target) sum(abs(target-ref))
  
  if(class(dm_model_list[[1]]) == "stanfit") {
    
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
    
  } else {
    
    B_ref = fit$par$B
    # collect from result list
    tb = lapply(jobs_ok,function(seed) {
      fit = dm_model_list[[seed]]
      B_target = fit$par$B
      sign_flip = 1
      if(distance(B_ref,B_target) > distance(B_ref,-B_target))
        sign_flip = -1
      B_target = sign_flip*B_target
      tibble(protein_name = protein_names,
             coeff = B_target[,2],
             seed = seed)
    }) %>% bind_rows
    
  }

}
