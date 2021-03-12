#' Heatmap of median marker expression
#'
#' @import dplyr
#' @import magrittr
#' @import RColorBrewer
#' @importFrom stats median
#' @importFrom pheatmap pheatmap
#' @importFrom grDevices colorRampPalette
#' @export
#'
#' @param df_samples Data frame or tibble with proteins counts,
#'   cell condition, and group information
#' @param sample_info_names Column names that contain information about
#'   the cell, e.g. donor, condition, file name, or cell type
#' @param protein_names A vector of column names of protein to use in the
#'   analysis
#' @param arrange_by_1 Column name
#' @param arrange_by_2 Column name
#' @param cluster_cols Apply hierarchical cluster to columns
#' @param fun Summary statistics of marker expression
#' @return \code{\link[pheatmap]{pheatmap}} object
#'
#' @examples
#' set.seed(23)
#' df <- generate_data()
#' protein_names <- names(df)[3:12]
#' df <- dplyr::mutate_at(df, protein_names, function(x) asinh(x/5))
#' CytoGLMM::plot_heatmap(df,
#'                        protein_names = protein_names,
#'                        sample_info_names = c("donor", "condition"),
#'                        arrange_by_1 = "condition")
plot_heatmap <- function(df_samples,
                         sample_info_names,
                         protein_names,
                         arrange_by_1,
                         arrange_by_2 = "",
                         cluster_cols = FALSE,
                         fun = median) {
    expr_median <- df_samples %>%
      group_by(.dots = sample_info_names) %>%
      summarise_at(protein_names,fun) %>%
      arrange_(arrange_by_1) %>%
      as.data.frame
    if(nchar(arrange_by_2) > 0) expr_median %<>% arrange_(arrange_by_2)
    df_expr_median <- as.data.frame(expr_median[,protein_names])
    rownames(df_expr_median) <- seq_len(nrow(expr_median))
    col_names <- arrange_by_1
    if(nchar(arrange_by_2) > 0) col_names <- c(arrange_by_1,arrange_by_2)
    df_annotation <- data.frame(expr_median[,col_names])
    names(df_annotation) <- col_names
    rownames(df_annotation) <- seq_len(nrow(expr_median))
    color <- colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(100)
    pheatmap(t(df_expr_median),
             color = color,
             clustering_method = "average",
             show_colnames = FALSE,
             cluster_cols = cluster_cols,
             cluster_rows = TRUE,
             annotation_col = df_annotation)
}
