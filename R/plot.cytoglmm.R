#' Plot fixded coefficients of random effects model
#'
#' @import ggplot2
#' @import tibble
#' @import magrittr
#' @import dplyr
#' @import cowplot
#' @export
#'
plot.cytoglmm = function(fit,order = FALSE) {

  if(class(fit) != "cytoglmm")
    stop("Input needs to be a cytoglmm object computed by cytoglmm function.")
  if(!is.factor(fit$glmmfit$y))
    stop("Currently only plotting results for logistic regression.")

  summ = summary(fit$glmmfit)

  # random effects
  stdev = sqrt(diag(summ$varcor[[1]])[-1])
  tb_random = tibble(protein_name = names(stdev),
                     stdev = stdev)

  # fixed effects
  alpha = 0.05
  ci = qnorm(1-alpha/2)
  tb_coeff = summ$coefficient
  tb_coeff = tb_coeff[-1,]
  tb_coeff %<>%
    as.data.frame %>%
    rownames_to_column(var = "protein_name") %>%
    as.tibble
  tb_coeff %<>%
    mutate(high = tb_coeff$Estimate+ci*tb_coeff$`Std. Error`,
           low = tb_coeff$Estimate-ci*tb_coeff$`Std. Error`)

  # order proteins according to coefficients
  if(order) {
    ind = sort.int(tb_coeff$Estimate,index.return=TRUE)$ix
    reordered_names = tb_coeff$protein_name[ind]
    tb_coeff$protein_name %<>% factor(levels = reordered_names)
    tb_random$protein_name %<>% factor(levels = reordered_names)
  }

  # plotting
  prandom = ggplot(tb_random, aes(x = stdev, y = protein_name)) +
    geom_point(size = 2) +
    ggtitle("Random Effects") +
    xlab("standard deviation") +
    theme(axis.title.y = element_blank())
  xlab_str = fit$df_samples_subset %>%
    pull(fit$condition) %>%
    levels %>%
    paste(collapse = " <-> ")
  pcoef = ggplot(tb_coeff, aes(x = Estimate, y = protein_name)) +
    geom_vline(xintercept = 0,color = "red") +
    geom_point(size = 2) +
    geom_errorbarh(aes(xmin = low, xmax = high)) +
    ggtitle("Fixed Effects") +
    xlab(xlab_str) +
    theme(axis.title.y = element_blank())

  plot_grid(prandom,pcoef)

}
