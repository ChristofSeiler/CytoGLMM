#' Plot fixded coefficients of group-specific fixed effects model
#'
#' @import ggplot2
#' @import tibble
#' @import magrittr
#' @import dplyr
#' @import cowplot
#' @export
#'
plot.cytogroup = function(fit) {

  if(class(fit) != "cytogroup")
    stop("Input needs to be a cytogroup object computed by cytogroup function.")

  summ = summary(fit$groupfit)
  xlab_str = fit$df_samples_subset %>%
    pull(fit$condition) %>%
    levels %>%
    paste(collapse = " <-> ")

  # donor-specific effects
  tb_coeff = tibble(names = names(fit$groupfit$coefficients),
                    value = fit$groupfit$coefficients)
  tb_coeff %<>%
    mutate(protein_name = sapply(str_split(tb_coeff$names,":"),function(x) x[1]),
           donor = sapply(str_split(tb_coeff$names,":"),function(x) x[2]))
  tb_coeff %<>% dplyr::filter(!is.na(donor))
  pdonor = ggplot(tb_coeff, aes(x = protein_name, y = value, group = donor, color = protein_name)) +
    geom_hline(yintercept = 0,color = "red") +
    geom_jitter(size = 1,height = 0.001,alpha = 0.5) +
    ggtitle("Donor-Specific Effects") +
    ylab(xlab_str) +
    coord_flip() +
    theme(legend.position="none")

  # fixed effects
  alpha = 0.05
  ci = qnorm(1-alpha/2)
  tb_coeff = summ$coefficient
  tb_coeff = tb_coeff[fit$protein_names,]
  tb_coeff %<>%
    as.data.frame %>%
    rownames_to_column(var = "protein_name") %>%
    as.tibble
  tb_coeff$Estimate %<>% as.character %>% as.numeric
  tb_coeff$`Std. Error` %<>% as.character %>% as.numeric
  tb_coeff$`z value` %<>% as.character %>% as.numeric
  tb_coeff$`Pr(>|z|)` %<>% as.character %>% as.numeric
  tb_coeff %<>%
    mutate(high = tb_coeff$Estimate+ci*tb_coeff$`Std. Error`,
           low = tb_coeff$Estimate-ci*tb_coeff$`Std. Error`)

  pcoef = ggplot(tb_coeff, aes(x = protein_name, y = Estimate)) +
    geom_hline(yintercept = 0,color = "red") +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = low, ymax = high)) +
    ggtitle("Fixed Effects") +
    ylab(xlab_str) +
    coord_flip()

  plot_grid(pdonor,pcoef)

}
