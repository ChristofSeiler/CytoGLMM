#' Plot all components of mixture regression
#'
#' @import ggplot2
#' @import tibble
#' @import magrittr
#' @import dplyr
#' @import cowplot
#' @import flexmix
#' @export
#'
plot.cytoflexmix = function(fit,k = NULL) {

  if(class(fit) != "cytoflexmix")
    stop("Input needs to be a cytoflexmix object computed by cytoflexmix function.")

  # plot selection criteria
  tb_sel = tibble(
    id = seq(fit$flexmixfits),
    k = sapply(fit$flexmixfits,function(fit) fit@components %>% length),
    BIC = sapply(fit$flexmixfits,BIC),
    #AIC = sapply(fit$flexmixfits,AIC)
    )
  # select best model
  best_id = tb_sel$id[which.min(tb_sel$BIC)]
  if(!is.null(k)) best_id = k
  pmodel = ggplot(tb_sel %>% gather(criterion,value,-c(k,id)),
                  aes(k,value)) +
    geom_vline(xintercept = tb_sel$k[best_id],color = "darkgray") +
    geom_hline(yintercept = tb_sel$BIC[best_id],color = "darkgray") +
    geom_point() +
    geom_line() +
    scale_x_continuous(breaks = fit$ks) +
    ggtitle("Model Selection (BIC)") #+
    #theme(legend.position = c(0.15, 0.15))
    #theme(legend.position = "top")

  # plot cluster sizes
  ct = table(pull(fit$df_samples_subset,fit$group),
             fit$flexmixfits[[best_id]]@cluster)
  ct[ct > 0] = 1
  tb_size = tibble(comp = as.factor(colnames(ct)),
                   size = colSums(ct))
  psize = ggplot(tb_size,
                 aes(size,comp,color = comp,label = size)) +
    geom_point(size = 2) +
    #geom_text() +
    scale_x_continuous(breaks = tb_size$size) +
    ggtitle("Cluster Size") +
    theme(legend.position="none")

  # plot component-wise coefficients
  xlab_str = fit$df_samples_subset %>%
    pull(fit$condition) %>%
    levels %>%
    paste(collapse = " <-> ")

  best_refit = refit(fit$flexmixfits[[best_id]],method = "mstep")
  alpha = 0.05
  ci = qnorm(1-alpha/2)
  tb_coeff_all = lapply(seq(tb_sel$k[best_id]),function(comp) {
    summ = summary(best_refit@components[[1]][[comp]])
    tb_coeff = summ$coefficient
    tb_coeff = tb_coeff[fit$protein_names,]
    tb_coeff %<>%
      as.data.frame %>%
      rownames_to_column(var = "protein_name") %>%
      as.tibble
    tb_coeff %<>%
      mutate(high = tb_coeff$Estimate+ci*tb_coeff$`Std. Error`,
             low = tb_coeff$Estimate-ci*tb_coeff$`Std. Error`)
    tb_coeff %<>% mutate(comp = comp)
    tb_coeff
  }) %>% bind_rows
  tb_coeff_all$comp %<>% as.factor

  peffects = ggplot(tb_coeff_all, aes(x = Estimate, y = protein_name, color = comp)) +
    geom_vline(xintercept = 0,color = "red") +
    geom_point(size = 2) +
    geom_errorbarh(aes(xmin = low, xmax = high)) +
    ggtitle("Fixed Mixture Effects") +
    xlab(xlab_str) +
    theme(legend.position="none",
          axis.title.y = element_blank())

  pleft = plot_grid(psize, pmodel,
                    nrow = 2, rel_heights = c(0.4,0.6), align = "v")

  plot_grid(pleft, NULL, peffects,
            ncol = 3,
            rel_widths = c(0.45,0.05,0.5))

}
