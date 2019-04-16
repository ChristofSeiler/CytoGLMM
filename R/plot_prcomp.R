#' Plot PCA of subsampled data using ggplot.
#'
#' @import ggplot2
#' @import dplyr
#' @import magrittr
#' @import cowplot
#' @import factoextra
#' @export
#'
plot_prcomp = function(df_samples,
                       protein_names,
                       color_var = "treatment",
                       subsample_size = 10000,
                       seed = 0xdada,
                       repel = TRUE) {
  set.seed(seed)
  n = min(table(df_samples[,color_var]),subsample_size)
  by_variable = df_samples %>% group_by_(color_var) %>% sample_n(n) %>% ungroup %>% as.data.frame
  res_pca = prcomp(by_variable[,protein_names],scale. = FALSE)
  explained_var = (100*res_pca$sdev^2/sum(res_pca$sdev^2)) %>% round(.,1)
  xlab_str = paste0("PC1 (",explained_var[1],"%)")
  ylab_str = paste0("PC2 (",explained_var[2],"%)")
  p1 = fviz_eig(res_pca,geom="bar") + ylab("Variance in %")
  p2 = fviz_pca_var(res_pca, col.var="contrib",
                    gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                    repel = repel) +
    coord_fixed(ratio = explained_var[2] / explained_var[1]) +
    ggtitle("Markers") +
    xlab(xlab_str) +
    ylab(ylab_str)
  by_variable %<>% add_column(PC1 = res_pca$x[,1])
  by_variable %<>% add_column(PC2 = res_pca$x[,2])
  p3 = ggplot(by_variable,aes(PC1,PC2)) +
    geom_density_2d(aes_string(col = color_var)) +
    coord_fixed(ratio = explained_var[2] / explained_var[1]) +
    xlab(xlab_str) +
    ylab(ylab_str) +
    ggtitle("Cells")
  plot_grid(p3, p1, p2, nrow = 2, rel_widths = c(0.7,0.3))
}
