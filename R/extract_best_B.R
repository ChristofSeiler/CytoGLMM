#' Extact best coefficient matrix B
#'
#' @import tibble
#' @import BiocParallel
#' @export
#'
extract_best_B = function(dm_model_list,protein_names) {

  # collect tables from result list
  jobs_ok = which(bpok(dm_model_list))
  tb_list = lapply(jobs_ok,function(job_id) extract_B(dm_model_list,
                                                      job_id,
                                                      protein_names))

  # measure error for each reference coefficient
  ## measure of scale using sample variance
  align_error = sapply(tb_list,function(tb) var(tb$coeff))
  # robust measure of scale using interquartile range
  #align_error = sapply(tb_list,function(tb)
  #  diff(quantile(tb$coeff,probs = c(0.25,0.75))))

  # select minimum
  tb_list[[which.min(align_error)]]
}
