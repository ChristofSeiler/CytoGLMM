#' Plot model selection to choose number optimal number of clusters
#'
#' @import ggplot2
#' @import tibble
#' @import magrittr
#' @import dplyr
#' @import cowplot
#' @import flexmix
#' @export
#'
plot_model_selection = function(fit,k = NULL) {

  if(!is(fit, "cytoflexmix"))
    stop("Input needs to be a cytoflexmix object computed by cytoflexmix function.")

  # plot selection criteria
  tb_sel = tibble(
    id = seq(fit$flexmixfits),
    k = vapply(fit$flexmixfits, function(fit) fit@components %>% length, numeric(1)),
    BIC = vapply(fit$flexmixfits, BIC, numeric(1)),
    AIC = vapply(fit$flexmixfits, AIC, numeric(1))
  )
  # select best model
  best_id = tb_sel$id[which.min(tb_sel$BIC)]
  if(!is.null(k)) best_id = k
  pmodel = ggplot(tb_sel %>% gather(criterion,value,-c(k,id)),
                  aes(k,value, color = criterion)) +
    geom_vline(xintercept = tb_sel$k[best_id],color = "darkgray") +
    geom_hline(yintercept = tb_sel$BIC[best_id],color = "darkgray") +
    geom_point() +
    geom_line() +
    scale_x_continuous(breaks = fit$ks) +
    xlab("number of clusters") +
    ggtitle("Model Selection")

  # # plot cluster sizes
  # ct = table(pull(fit$df_samples_subset,fit$group),
  #            fit$flexmixfits[[best_id]]@cluster)
  # ct[ct > 0] = 1
  # tb_size = tibble(comp = as.factor(colnames(ct)),
  #                  size = colSums(ct))
  # psize = ggplot(tb_size,
  #                aes(size,comp,color = comp,label = size)) +
  #   geom_point(size = 2) +
  #   scale_x_continuous(breaks = tb_size$size) +
  #   xlab("cluster size") +
  #   ylab("cluster label") +
  #   theme(legend.position="none") +
  #   ggtitle("Cluster Assignment")

  plot_grid(pmodel)

}
