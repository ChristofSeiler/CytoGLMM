#' Plot two association matrices for both stimulation conditions.
#'
#' @import ggplot2
#' @import reshape2
#' @import gridExtra
#' @export
#'
plot_pairwise_mi = function(fit,
                            df_samples,
                            protein_names) {
  # compute posterior median of pi and beta
  pi = rstan::extract(fit,pars = "pi")[[1]]
  pi_median = apply(X = pi,MARGIN = c(2,3,4),FUN = function(x) quantile(x, probs = 0.5))
  beta = rstan::extract(fit,pars = "beta")[[1]]
  beta_median = apply(X = beta,MARGIN = c(2,3),FUN = function(x) quantile(x, probs = 0.5))
  # compute pairwise mutual informaiton for both conditions
  x_list = list(x_0,x_1)
  M_list = lapply(x_list,function(x) {
    eta = beta_median %*% x
    theta = exp(eta)/sum(exp(eta))
    M = pairwise_mi(theta,pi_median)
    rownames(M) = colnames(M) = protein_names
    M
  })
  # plot both conditions side-by-side
  titles = rownames(contrasts(df_samples$condition))
  get_upper_tri = function(cormat) {
    cormat[lower.tri(cormat)] = NA
    diag(cormat) = 0 # doesn't give zero exactly
    cormat
  }
  ggobj_list = lapply(1:length(M_list),function(i) {
    M_long = melt(get_upper_tri(M_list[[i]]), na.rm = TRUE)
    ggplot(data = M_long, aes(Var2, Var1, fill = value)) +
      geom_tile(color = "white") +
      scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                           midpoint = 0, space = "Lab",
                           limit = range(cbind(M_list[[1]],M_list[[2]])),
                           name = "Mutual information") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            legend.position="bottom") +
      coord_fixed() +
      ggtitle(titles[i])
  })
  legend = get_legend(ggobj_list[[1]])
  ggobj_list = lapply(ggobj_list,function(ggobj) {
    ggobj + theme(legend.position = "none")
  })
  grid.arrange(ggobj_list[[1]], ggobj_list[[2]], legend,
               ncol=2, nrow = 2,
               layout_matrix = rbind(c(1,2), c(3,3)),
               widths = c(2.7, 2.7), heights = c(2.5, 0.3))
}
