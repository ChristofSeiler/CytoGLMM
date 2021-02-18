#' Generalized linear mixed model with maximum likelihood.
#'
#' @import mbest
#' @import doParallel
#' @import dplyr
#' @export
#'
glmm_moment = function(df_samples,
                       protein_names,
                       response,
                       group = "donor",
                       covariate_names = NULL,
                       num_cores = 1) {
  parallel_str = "FALSE"
  if(num_cores > 1) {
    registerDoParallel(cores = num_cores)
    parallel_str = "TRUE"
  }
  markers_str = paste0(c(protein_names, covariate_names), collapse = " + ")
  formula_expr = NULL
  if( is.factor(pull(df_samples,response)) ) {
    formula_expr = parse(text = paste0("mhglm(",
                                       paste(response,"~",markers_str,"+",paste0("(",markers_str," | ",group,"), ")),
                                       "family = binomial(link='logit'), ",
                                       "data = df_samples, ",
                                       "control = mhglm.control(parallel = ", parallel_str, ", ",
                                       "fit.method = 'firthglm.fit'))"))
  } else {
    formula_expr = parse(text = paste0("mhglm(",
                                       paste(response,"~",markers_str,"+",paste0("(",markers_str," | ",group,"), ")),
                                       "family = gaussian(link = 'identity'), ",
                                       "data = df_samples, ",
                                       "control = mhglm.control(parallel = ", parallel_str, ", ",
                                       "fit.method = 'firthglm.fit'))"))
  }
  eval(formula_expr)
}
