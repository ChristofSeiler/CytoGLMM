#' Extact and plot noise term
#'
#' @import ggplot2
#' @import tibble
#' @import magrittr
#' @import dplyr
#' @import cowplot
#' @export
#'
plot_sigma = function(dm_model_list,
                  protein_names = protein_names) {

  # some jobs may fail (because of computing cluster instabilities)
  if(length(dm_model_list) == 0)
    stop("no jobs completed successfully")
  
  if(class(dm_model_list[[1]]) == "stanfit") {
    stop("plotting for HMC samples not implemented yet")
  }
  
  # collect from result list
  tb = lapply(seq_along(dm_model_list),function(run) {
    fit = dm_model_list[[run]]
    sigma = fit$par$sigma
    tibble(protein_name = protein_names,
           sigma = sigma,
           run = run)
  }) %>% bind_rows
  
  # combine two plots
  tb$run %<>% as.factor
  ggplot(tb,aes(x = protein_name,y = sigma,color = protein_name)) +
    geom_jitter(size = 1,height = 0.001,alpha = 0.5) +
    coord_flip() +
    ylab("posterior median") +
    theme(legend.position="none") +
    ggtitle("sigma")
  
}
