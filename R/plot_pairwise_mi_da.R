#' Plot association matrix difference between stimulation conditions.
#'
#' @import ggplot2
#' @import reshape2
#' @import gridExtra
#' @import doParallel
#' @export
#'
plot_pairwise_mi_da = function(fit,
                               df_samples_binned,
                               protein_names) {
  registerDoParallel()
  # compute posterior median of pi and beta
  pi = rstan::extract(fit,pars = "pi")[[1]]
  beta = rstan::extract(fit,pars = "beta")[[1]]
  # compare pairwise mutual informaiton between conditions
  x_list = list(x_0,x_1)
  M_da_list = foreach(i = 1:dim(pi)[1]) %dopar% {
    M_list = lapply(x_list,function(x) {
      eta = beta[i,,] %*% x
      theta = exp(eta)/sum(exp(eta))
      M = pairwise_mi(theta,pi[i,,,])
      rownames(M) = colnames(M) = protein_names
      M
    })
    ifelse(M_list[[1]] > M_list[[2]],1,0)
  }
  M_da = Reduce('+', M_da_list)
  M_probs = ifelse(test = M_da == 0,
                   yes = 1/dim(pi)[1],
                   no = M_da/dim(pi)[1])
  # plot probabilities
  title = paste0("P(",
                 paste(rownames(contrasts(df_samples_binned$condition)),
                       collapse = " > "),
                 ")")
  get_upper_tri = function(cormat) {
    cormat[lower.tri(cormat)] = NA
    diag(cormat) = 0 # doesn't give zero exactly
    cormat
  }
  M_long = melt(get_upper_tri(M_probs), na.rm = TRUE)
  ggplot(data = M_long, aes(Var2, Var1, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                         midpoint = 0, space = "Lab",
                         limit = c(0,1),
                         name = "Posterior\nprobability") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1),
          axis.title.x=element_blank(),
          axis.title.y=element_blank()) +
    coord_fixed() +
    ggtitle(title)
}
