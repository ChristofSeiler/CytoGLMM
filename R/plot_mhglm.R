#' Plot GLMM fixed effects coefficients.
#'
#' @export
#'
plot.mhglm <- function(res_glmm_ml) {
  num_explanatory = nrow(summary(res_glmm_ml)$coefficient)-1
  adj_alpha = 0.05/num_explanatory
  ci = qnorm(1-adj_alpha/2)
  coeff = summary(res_glmm_ml)$coefficient[2:nrow(summary(res_glmm_ml)$coefficient),]
  # add protein name
  coeff = data.frame(coeff,protein_name=rownames(coeff))
  # sort according to z-value
  coeff$protein_name = factor(coeff$protein_name,
                              levels = coeff$protein_name[order(coeff$z.value,
                                                                decreasing = FALSE)])
  coeff = data.frame(coeff,
                     high=coeff$Estimate+ci*coeff$Std..Error,
                     low=coeff$Estimate-ci*coeff$Std..Error)
  lab = ""
  title = ""
  if(is.factor(res_glmm_ml$y)) {
    condition_names = rownames(contrasts(res_glmm_ml$y))
    lab = paste0("Log-Odds (",condition_names[1]," <--> ",condition_names[2],")")
    title = paste(condition_names,collapse = " vs. ")
  } else {
    lab = "arcsinh transformed counts"
    title = paste0("Prediction of ",
                   strsplit(as.character(res_glmm_ml$call)," ")[[2]][1])
  }
  ggplot(coeff, aes(x=protein_name, y=Estimate)) +
    geom_hline(yintercept = 0, col = "red") +
    geom_point() +
    geom_segment(mapping=aes(x=protein_name, y=low, xend=protein_name, yend=high)) +
    ylab(lab) +
    labs(title = title) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title.y = element_blank()) +
    coord_flip()
}
