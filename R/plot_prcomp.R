#' Plot PCA of subsampled data using ggplot.
#'
#' @import ggfortify
#' @export
#'
plot_prcomp <- function(df_samples,
                        protein_names,
                        color_var = "condition",
                        subsample_size = 10000,
                        seed = 1234) {
  set.seed(seed)
  subsample_ids = sample(x = nrow(df_samples),
                         size = subsample_size,
                         replace = FALSE)
  df_samples_subset = df_samples[subsample_ids,]
  table(df_samples_subset$donor)
  res_pca = prcomp(df_samples_subset[,protein_names],scale. = FALSE)
  explained_var = (100*res_pca$sdev^2/sum(res_pca$sdev^2)) %>% round(.,1)
  autoplot(res_pca,
           data = df_samples[subsample_ids,],
           colour = color_var,
           alpha = 0.5,
           loadings = TRUE,
           loadings.label = TRUE,
           loadings.colour = "black",
           loadings.label.colour = "black",
           loadings.label.repel = TRUE,
           xlab = paste0("PC1 (",explained_var[1],"%)"),
           ylab = paste0("PC2 (",explained_var[2],"%)")) +
    coord_fixed(ratio = explained_var[2] / explained_var[1])
}
