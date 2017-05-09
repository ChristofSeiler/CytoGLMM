#' Plot GLMM fixed effects coefficients.
#'
#' @export
#'
plot_glmm_ml = function(res_glmm_ml,
                        df_samples,
                        title) {
  num_explanatory = nrow(summary(res_glmm_ml)$coefficient)-1
  adj_alpha = 0.05/num_explanatory
  ci = qnorm(1-adj_alpha/2)
  coeff = summary(res_glmm_ml)$coefficient[2:nrow(summary(res_glmm_ml)$coefficient),]
  # add protein name
  coeff = data.frame(coeff,protein_name=rownames(coeff))
  # sort according to z-value
  coeff$protein_name = factor(coeff$protein_name,
                              levels = coeff$protein_name[order(abs(coeff$z.value),
                                                                decreasing = FALSE)])
  coeff = data.frame(coeff,
                     high=coeff$Estimate+ci*coeff$Std..Error,
                     low=coeff$Estimate-ci*coeff$Std..Error)
  response_name = strsplit(as.character(summary(res_glmm_ml)$call[2]),
                           split = "~") %>% unlist %>% .[1]
  #contrasts(df_samples$condition)
  xlab_str = "Log-Odds (H1 <--> H3)"
  ggplot(coeff, aes(x=protein_name, y=Estimate)) +
    geom_hline(yintercept = 0, col = "red") +
    geom_point() +
    geom_segment(mapping=aes(x=protein_name, y=low, xend=protein_name, yend=high)) +
    ylab(xlab_str) +
    xlab("Proteins") +
    labs(title = response_name) +
    theme(plot.title = element_text(hjust = 0.5)) +
    coord_flip()
}
