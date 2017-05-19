#' Plot association matrix difference between stimulation conditions.
#'
#' @import ggplot2
#' @import reshape2
#' @import gridExtra
#' @import foreach
#' @import doParallel
#' @export
#'
plot_pairwise_mi_da = function(fit,
                               df_samples,
                               protein_names) {

  # for parallelizing over posterior samples
  registerDoParallel()

  # compute posterior median of pi and beta
  pi = rstan::extract(fit,pars = "pi")[[1]]
  beta = rstan::extract(fit,pars = "beta")[[1]]

  # compare pairwise mutual information between conditions
  x_list = list(c(1,0),c(1,1))
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

  # posterior expected FDR
  # see Mitra, Mueller, and Ji (2016), Bayesian Graphical Models for Differential Pathways
  fdr_threshold = function(FDR) {
    P = M_da / dim(pi)[1]
    seq_kappa = seq(0.5,1,1/(length(M_da_list)/10))
    seq_kappa = kappa_seq[-length(seq_kappa)]
    df_kappa = sapply(seq_kappa,function(kappa) {
      ps = P[lower.tri(P)]
      I = ifelse(test = ps>kappa,yes = 1,no = 0)
      FDR_kappa = sum((1-ps)*I)/sum(I)
      c(kappa=kappa,FDR_kappa=FDR_kappa)
    }) %>% t %>% data.frame
    kappa = df_kappa[which(df_kappa$FDR_kappa <= FDR)[1],"kappa"]
    P_thres = ifelse(P > kappa,yes = FDR,no = NA)
    P_thres
  }
  fdrs = c(0.1,0.01)
  otherwise = paste0(">",fdrs[1])
  fdr_list = lapply(fdrs,fdr_threshold)
  M_fdr = do.call(pmin,c(fdr_list,na.rm = TRUE))
  M_fdr[is.na(M_fdr)] = otherwise

  # plot pairwise matrix
  title = paste0(paste(rownames(contrasts(df_samples$condition)),
                       collapse = " > ")," | data")
  get_upper_tri = function(cormat) {
    cormat[lower.tri(cormat)] = NA
    diag(cormat) = otherwise
    cormat
  }
  M_long = melt(get_upper_tri(M_fdr), na.rm = TRUE)
  M_long$value = factor(M_long$value,levels = c(otherwise,fdrs))
  ggplot(data = M_long, aes(Var2, Var1, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_brewer(name = "FDR",palette = "Greens") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1),
          axis.title.x=element_blank(),
          axis.title.y=element_blank()) +
    coord_fixed() +
    ggtitle(title)
}
