#' Plot all components of mixture regression
#'
#' @aliases plot.cytoflexmix
#' @method plot cytoflexmix
#'
#' @import ggplot2
#' @import tibble
#' @import magrittr
#' @import dplyr
#' @import flexmix
#' @importFrom cowplot plot_grid
#' @importFrom caret cluster
#' @export
#'
#' @param x A \code{cytoflexmix} class
#' @param k Number of clusters
#' @param separate create two separate \code{\link[ggplot2]{ggplot2}} objects
#' @param ... Other parameters
#' @return \code{\link[ggplot2]{ggplot2}} object
#'
#' @examples
#' set.seed(23)
#' df = generate_data()
#' protein_names = names(df)[3:12]
#' df = dplyr::mutate_at(df, protein_names, function(x) asinh(x/5))
#' mix_fit = CytoGLMM::cytoflexmix(df,
#'                                 protein_names = protein_names,
#'                                 condition = "condition",
#'                                 group = "donor",
#'                                 ks = 2)
#' plot(mix_fit)
plot.cytoflexmix = function(x, k = NULL, separate = FALSE, ...) {

  if(!is(x, "cytoflexmix"))
    stop("Input needs to be a cytoflexmix object computed by cytoflexmix function.")

  # plot selection criteria
  tb_sel = tibble(
    id = seq(x$flexmixfits),
    k = vapply(x$flexmixfits,function(fit) fit@components %>% length, numeric(1)),
    BIC = vapply(x$flexmixfits, BIC, numeric(1))
    )
  # select best model
  best_id = tb_sel$id[which.min(tb_sel$BIC)]
  if(!is.null(k)) best_id = k

  # plot cell to cluster assignments
  tb = tibble(
    group = pull(x$df_samples_subset,x$group),
    cluster = x$flexmixfits[[best_id]]@cluster
  )
  tb$cluster %<>% as.factor
  tb_tally = tb %>%
    group_by(group,cluster) %>%
    tally() %>%
    arrange(cluster)
  tb_tally$group %<>% factor(levels = rev(tb_tally$group))
  pceltocluster = ggplot(tb_tally, aes(x = group, y = n, fill = cluster)) +
    geom_bar(stat="identity", position = "dodge") +
    xlab(x$group) +
    ylab("number of cells") +
    coord_flip() +
    ggtitle("Cluster Assigment")

  # plot component-wise coefficients
  xlab_str = x$df_samples_subset %>%
    pull(x$condition) %>%
    levels %>%
    paste(collapse = " <-> ")

  best_refit = refit(x$flexmixfits[[best_id]],method = "mstep")
  alpha = 0.05
  ci = qnorm(1-alpha/2)
  tb_coeff_all = lapply(seq(tb_sel$k[best_id]),function(comp) {
    summ = summary(best_refit@components[[1]][[comp]])
    tb_coeff = summ$coefficient
    tb_coeff = tb_coeff[x$protein_names,]
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

  if(separate) {
    return(list(pceltocluster=pceltocluster,peffects=peffects))
  } else {
    plot_grid(pceltocluster,peffects)
  }

}
