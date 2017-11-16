#' Plot fixded coefficients of random effects model
#'
#' @import ggplot2
#' @import tibble
#' @import magrittr
#' @import dplyr
#' @import cowplot
#' @export
#'
plot.cytoglmm = function(fit) {

  if(class(fit) != "cytoglmm")
    stop("Input needs to be a cytoglmm object computed by cytoglmm function.")
  if(!is.factor(fit$glmmfit$y))
    stop("Currently only plotting results for logistic regression.")

  summ = summary(fit$glmmfit)

  # random effects
  stdev = sqrt(diag(summ$varcor$donor)[-1])
  tb_random = tibble(protein_name = names(stdev),
                     stdev = stdev)
  prandom = ggplot(tb_random, aes(x = protein_name, y = stdev)) +
    geom_point(size = 2) +
    ggtitle("Random Effects") +
    ylab("standard deviation") +
    coord_flip()

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
  xlab_str = fit$df_samples_subset %>%
    pull(fit$condition) %>%
    levels %>%
    paste(collapse = " <-> ")
  pcoef = ggplot(tb_coeff, aes(x = protein_name, y = Estimate)) +
    geom_hline(yintercept = 0,color = "red") +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = low, ymax = high)) +
    ggtitle("Fixed Effects") +
    ylab(xlab_str) +
    coord_flip()

  plot_grid(prandom,pcoef)

}
