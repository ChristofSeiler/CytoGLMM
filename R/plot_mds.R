#' MDS on median marker expression.
#'
#' @import ggplot2
#' @import dplyr
#' @import magrittr
#' @import tibble
#' @export
#'
plot_mds = function(df_samples,
                    protein_names,
                    sample_info_names,
                    color,
                    sample_label = "") {
  expr_median = df_samples %>%
    group_by(.dots = sample_info_names) %>%
    summarise_at(protein_names,median)
  dist_matrix = dist(expr_median[,-seq(sample_info_names)])
  mds_res = cmdscale(dist_matrix,eig = TRUE, k = 2) # k is the number of dim
  explained_var = (100*mds_res$eig[1:2]/sum(mds_res$eig)) %>% round(digits = 1)
  expr_median %<>% bind_cols(tibble(MDS1 = mds_res$points[,1],
                                    MDS2 = mds_res$points[,2]))
  
  # protein_sd = apply(expr_median[,protein_names],2,sd)
  # protein_selection = protein_names[protein_sd != 0]
  # expr_cor = cor(expr_median[,protein_selection],expr_median[,c("MDS1","MDS2")])
  # ranking = apply(expr_cor,1,function(expr) max(abs(expr))) %>%
  #   order(decreasing = TRUE)
  # expr_cor[ranking,]

  gg = ggplot(expr_median, aes_string(x = "MDS1", y = "MDS2",color = color)) +
    geom_point(size = 2) +
    coord_fixed(ratio = explained_var[2] / explained_var[1]) +
    xlab(paste0("MDS1 (",explained_var[1],"%)")) +
    ylab(paste0("MDS2 (",explained_var[2],"%)"))
  if(nchar(sample_label) > 1)
    gg  = gg + geom_label(aes_string(label = sample_label))
  gg
}
