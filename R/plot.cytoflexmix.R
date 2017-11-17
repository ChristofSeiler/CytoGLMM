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
plot.cytoflexmix = function(fit, model_selection = FALSE) {

  if(class(fit) != "cytoflexmix")
    stop("Input needs to be a cytoflexmix object computed by cytoflexmix function.")

  # select best model
  BICs = sapply(fit$flexmixfits,BIC)
  AICs = sapply(fit$flexmixfits,AIC)
  best = which.min(BICs)

  # plot selection criteria
  pmodel = ggplot(tibble(k = fit$ks,
                         BIC = BICs,
                         AIC = AICs) %>% gather(criterion,value,-k),
                  aes(k,value,color = criterion)) +
    geom_vline(xintercept = best,color = "darkgray") +
    geom_hline(yintercept = BICs[best],color = "darkgray") +
    geom_point() +
    geom_line() +
    scale_x_continuous(breaks = fit$ks) +
    ggtitle("Model Selection") +
    theme(legend.position = c(0.8, 0.8))

  # plot cluster sizes
  ct = table(fit$df_samples_subset$donor,
             fit$flexmixfits[[best]]@cluster)
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

  best_refit = refit(fit$flexmixfits[[best]],method = "mstep")
  alpha = 0.05
  ci = qnorm(1-alpha/2)
  tb_coeff_all = lapply(seq(fit$ks[best]),function(comp) {
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
  peffects = ggplot(tb_coeff_all, aes(x = protein_name, y = Estimate, color = comp)) +
    geom_hline(yintercept = 0,color = "red") +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = low, ymax = high)) +
    ggtitle("Fixed Mixture Component Effects") +
    ylab(xlab_str) +
    coord_flip()

  pleft = plot_grid(psize, pmodel,
                    nrow = 2, rel_heights = c(0.3,0.7), align = "v")
  plot_grid(pleft, peffects,
            rel_widths = c(0.4,0.6))

}
