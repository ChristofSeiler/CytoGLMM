#' Heatmap of median marker expression.
#'
#' @import dplyr
#' @import magrittr
#' @import pheatmap
#' @export
#'
plot_heatmap = function(df_samples,
                        sample_info_names,
                        protein_names,
                        arrange_by,
                        rownames) {
  expr_median = df_samples %>%
    group_by(.dots = sample_info_names) %>%
    summarise_at(protein_names,median) %>%
    arrange_(arrange_by)
  df_expr_median = as.data.frame(expr_median[,protein_names])
  rownames(df_expr_median) = pull(expr_median,rownames)
  df_annotation = data.frame(arrange_by = expr_median[,arrange_by])
  rownames(df_annotation) = pull(expr_median,rownames)
  color = colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(100)
  pheatmap(t(df_expr_median),
           color = color,
           clustering_method = "average",
           show_colnames = FALSE,
           cluster_cols = FALSE,
           cluster_rows = TRUE,
           annotation_col = df_annotation)
}
