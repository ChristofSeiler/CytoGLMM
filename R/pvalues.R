#' Extract and adjust p-values from GLMM.
#'
#' @export
#'
pvalues = function(fit) {

  if(class(fit) == "mhglm") {

    pvalues_unadj = summary(fit)$coefficients[-1,4]
    pvalues_adj = p.adjust(pvalues_unadj,method = "BH")
    df_pvalues = data.frame(protein_names = names(pvalues_unadj),
                            pvalues_unadj,
                            pvalues_adj)
    df_pvalues$protein_names = as.character(df_pvalues$protein_names)
    df_pvalues = df_pvalues[order(df_pvalues$pvalues_unadj),]
    rownames(df_pvalues) = NULL
    df_pvalues

  } else {

    stop("not implemented yet")

  }

}
