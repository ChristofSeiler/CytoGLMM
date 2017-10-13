#' Extact coefficient matrix A
#'
#' @import ggplot2
#' @import tibble
#' @import magrittr
#' @import dplyr
#' @import cowplot
#' @export
#'
plot_coeff = function(tb,title_str,xlab_str) {

  tb$run %<>% as.factor

  # median over seeds
  tb_summary = tb %>%
    group_by(protein_name) %>%
    summarize(
      coeff_median = median(coeff),
      min_median = quantile(coeff,probs = 0.05),
      max_median = quantile(coeff,probs = 0.95)
    )

  # plot all bootstrap runs
  pall = ggplot(tb, aes(x = protein_name, y = coeff,color = protein_name)) +
    geom_hline(yintercept = 0,color = "red") +
    geom_jitter(size = 1,height = 0.001,alpha = 0.5) +
    #geom_errorbar(aes(ymin = min, ymax = max)) +
    ggtitle(title_str) +
    ylab(xlab_str) +
    coord_flip() +
    theme(legend.position="none")

  # plot summary
  psummary = ggplot(tb_summary, aes(x = protein_name, y = coeff_median)) +
    geom_hline(yintercept = 0,color = "red") +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = min_median, ymax = max_median)) +
    ggtitle("") +
    ylab(xlab_str) +
    coord_flip()

  # combine
  plot_grid(pall,psummary)
}
