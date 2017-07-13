#' Generalized linear mixed model with maximum likelihood.
#'
#' @import mbest
#' @import doParallel
#' @export
#'
glmm_ml <- function(df_samples,
                    protein_names,
                    response,
                    random_var = "donor",
                    cores = detectCores()) {
  registerDoParallel(cores = cores)
  markers_str = paste0(protein_names,collapse = " + ")
  formula_expr = NULL
  if( is.factor(df_samples[,response]) ) {
    formula_expr = parse(text = paste0("mhglm(",
                                       paste(response,"~",markers_str,"+",paste0("(",markers_str," | ",random_var,"),")),
                                       "family = binomial(link='logit'),",
                                       "data = df_samples,",
                                       "control = mhglm.control(parallel = TRUE,fit.method = 'firthglm.fit'))"))
  } else {
    formula_expr = parse(text = paste0("mhglm(",
                                       paste(response,"~",markers_str,"+",paste0("(",markers_str," | ",random_var,"),")),
                                       "family = gaussian(link = 'identity'),",
                                       "data = df_samples,",
                                       "control = mhglm.control(parallel = TRUE,fit.method = 'firthglm.fit'))"))
  }
  eval(formula_expr)
}
