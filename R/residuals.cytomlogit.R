#' Extact donor random effect z
#'
#' @import ggplot2
#' @import tibble
#' @import magrittr
#' @import dplyr
#' @import tidyr
#' @export
#'
residuals.cytomlogit = function(fit,dispersion = TRUE) {

  donors = fit$df_samples_subset %>%
    group_by_("donor",fit$condition) %>%
    tally()
  df_sub = fit$df_samples_subset %>%
    group_by(donor) %>%
    sample_n(min(donors$n)) %>%
    ungroup

  model_fit = fit$model_fit_list[[1]]
  Y = df_sub %>%
    select(fit$protein_names) %>%
    as.matrix
  X = model.matrix(as.formula(paste("~",fit$condition)), data = df_sub)
  donor = df_sub$donor %>% as.factor %>% as.numeric
  theta_raw = X %*% t(model_fit$par$A) + model_fit$par$z[donor,]
  n_cells = nrow(Y)
  title_str = "Residuals without Dispersion"
  if(dispersion) {
    title_str = "Residuals with Dispersion"
    dispersion = rnorm(n = n_cells*ncol(Y),
                       sd = rep(model_fit$par$sigma,n_cells)) %>%
      matrix(ncol = n_cells) %>% t
    theta_raw = theta_raw + dispersion
  }
  softmax = function(x) exp(x)/sum(exp(x))
  theta = softmax(theta_raw)
  Ypred = matrix(0,nrow = n_cells,ncol = ncol(Y))
  for(i in 1:n_cells)
    Ypred[i,] = rmultinom(n = 1,
                          size = sum(Y[i,]),
                          prob = theta[i,])
  Yres = Y-Ypred
  non_marker_cols = which(!names(df_sub) %in% fit$protein_names)
  df_residuals = bind_cols(df_sub[,non_marker_cols],as.tibble(Yres))

  df_residuals_long = df_residuals %>%
    gather_("key","value",fit$protein_names)
  df_residuals_long %<>% dplyr::filter(value < 20 & value > -20)

  ggplot(df_residuals_long,aes(value, fill = key)) +
    geom_histogram(bins = 30) +
    facet_wrap(~key) +
    theme(legend.position="none") +
    ggtitle(title_str)

}
