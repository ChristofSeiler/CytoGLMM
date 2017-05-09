#' Generalized linear mixed model with maximum likelihood.
#'
#' @import mbest
#' @import doParallel
#' @export
#'
glmm_ml <- function(df_samples,
                    protein_names) {
  registerDoParallel()
  markers_str = paste0(protein_names,collapse = " + ")
  formula_expr = parse(text = paste0("mhglm(",
                                     paste("condition ~",markers_str,"+",paste0("(",markers_str," | donor),")),
                                     "family = binomial(link='logit'),",
                                     "data = df_samples,",
                                     "control = mhglm.control(parallel = TRUE,fit.method = 'firthglm.fit'))"))
  formula_expr
  eval(formula_expr)
}
