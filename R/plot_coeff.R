#' Helper function to plot regression coeffcient
#'
#' @import ggplot2
#' @import tibble
#' @import magrittr
#' @import dplyr
#' @import cowplot
#'
plot_coeff = function(tb,title_str, title_str_right, xlab_str, redline = 0,
                      order = FALSE, separate = FALSE) {

  tb$run %<>% as.factor

  # median over seeds
  tb_summary = tb %>%
    group_by(protein_name) %>%
    summarize(
      coeff_median = median(coeff),
      min_median = quantile(coeff,probs = 0.05),
      max_median = quantile(coeff,probs = 0.95)
    )
  if(order) {
    ind = sort.int(tb_summary$coeff_median,index.return=TRUE)$ix
    reordered_names = tb_summary$protein_name[ind]
    tb_summary$protein_name %<>% factor(levels = reordered_names)
    tb$protein_name %<>% factor(levels = reordered_names)
  }

  # plot all bootstrap runs
  pall = ggplot(tb, aes(x = coeff, y = protein_name, color = protein_name)) +
    geom_vline(xintercept = redline,color = "red") +
    geom_jitter(size = 1,width = 0.001,alpha = 0.5) +
    ggtitle(title_str) +
    xlab(xlab_str) +
    theme(legend.position="none",
          axis.title.y = element_blank())

  # plot summary
  psummary = ggplot(tb_summary, aes(x = coeff_median, y = protein_name)) +
    geom_vline(xintercept = redline,color = "red") +
    geom_point(size = 2) +
    geom_errorbarh(aes(xmin = min_median, xmax = max_median)) +
    ggtitle(title_str_right) +
    xlab(xlab_str) +
    theme(axis.title.y = element_blank())

  # combine
  if(separate) {
    return(list(pall=pall,psummary=psummary))
  } else {
    psummary
  }

}
