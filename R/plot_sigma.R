#' Extact coefficient matrix A
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

  tb = NULL
  if(class(dm_model_list[[1]]) == "stanfit") {

    # collect from result list
    tb = lapply(seq_along(dm_model_list),function(run) {
      fit = dm_model_list[[run]]
      post = rstan::extract(fit)[["sigma"]]
      sigma = apply(post,c(2,3),median)
      tibble(protein_name = protein_names,
             sigma_median = apply(sigma,2,median),
             sigma_min = apply(sigma,2,function(x) quantile(x,probs = 0.05)),
             sigma_max = apply(sigma,2,function(x) quantile(x,probs = 0.95)),
             run = run)
    }) %>% bind_rows

  } else {

    # collect from result list
    tb = lapply(seq_along(dm_model_list),function(run) {
      fit = dm_model_list[[run]]
      sigma = fit$par$sigma
      tibble(protein_name = protein_names,
             sigma_median = apply(sigma,2,median),
             sigma_min = apply(sigma,2,function(x) quantile(x,probs = 0.05)),
             sigma_max = apply(sigma,2,function(x) quantile(x,probs = 0.95)),
             run = run)
    }) %>% bind_rows

  }

  # combine two plots
  tb$run %<>% as.factor
  p1 = ggplot(tb,aes(x = protein_name,y = sigma_median,color = protein_name)) +
    geom_jitter(size = 1,height = 0.001,alpha = 0.5) +
    coord_flip() +
    ylab("median") +
    theme(legend.position="none") +
    ggtitle("sigma")
  p2 = ggplot(tb,aes(x = protein_name,y = sigma_max-sigma_min,color = protein_name)) +
    geom_jitter(size = 1,height = 0.001,alpha = 0.5) +
    coord_flip() +
    ylab("95% - 5% quantile") +
    theme(legend.position="none") +
    ggtitle("")
  plot_grid(p1,p2,align = "h")

}
