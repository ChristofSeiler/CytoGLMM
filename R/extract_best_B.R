#' Extact best coefficient matrix B
#'
#' @import tibble
#' @export
#'
extract_best_B = function(dm_model_list,protein_names) {

  # some jobs may fail (because of computing cluster instabilities)
  if(length(dm_model_list) == 0)
    stop("no jobs completed successfully")

  # collect tables from result list
  tb_list = lapply(seq_along(dm_model_list),function(run) extract_B(dm_model_list,
                                                                    run,
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
