#' MDS on median marker expression
#'
#' @import ggplot2
#' @import dplyr
#' @import magrittr
#' @import tibble
#' @importFrom cowplot plot_grid
#' @export
#'
#' @param df_samples Data frame or tibble with proteins counts,
#'   cell condition, and group information
#' @param protein_names A vector of column names of protein to use in the analysis
#' @param sample_info_names Column names that contain information about
#'   the cell, e.g. donor, condition, file name, or cell type
#' @param color Column name
#' @param sample_label Column name
#' @return \code{\link[cowplot]{cowplot}} object
#'
#' @examples
#' set.seed(23)
#' df = generate_data()
#' protein_names = names(df)[3:12]
#' df = dplyr::mutate_at(df, protein_names, function(x) asinh(x/5))
#' CytoGLMM::plot_mds(df,
#'                    protein_names = protein_names,
#'                    sample_info_names = c("donor", "condition"),
#'                    color = "condition")
plot_mds = function(df_samples,
                    protein_names,
                    sample_info_names,
                    color,
                    sample_label = "") {

  # compute median marker expression
  expr_median = df_samples %>%
    group_by(.dots = sample_info_names) %>%
    summarise_at(protein_names,median) %>%
    as.data.frame
  dist_matrix = dist(expr_median[,-seq(sample_info_names)])
  mds_res = cmdscale(dist_matrix,eig = TRUE, k = 2) # k is the number of dim
  explained_var = (100*mds_res$eig[seq_len(2)]/sum(mds_res$eig)) %>% round(digits = 1)
  expr_median %<>% bind_cols(tibble(MDS1 = mds_res$points[,1],
                                    MDS2 = mds_res$points[,2]))

  # make circle of correlation plot
  protein_sd = apply(expr_median[,protein_names],2,sd)
  # only keep makers that have some variability
  protein_selection = protein_names[protein_sd != 0]
  # correlations between variables and MDS axes
  expr_cor = cor(expr_median[,protein_selection],expr_median[,c("MDS1","MDS2")]) %>% as.tibble
  expr_cor %<>% add_column(protein_selection)
  # add arrows coordinates
  expr_cor %<>% add_column(x0 = rep(0,nrow(expr_cor)))
  expr_cor %<>% add_column(y0 = rep(0,nrow(expr_cor)))
  # prepare circle
  circle = function(center = c(0, 0), npoints = 100) {
    r = 1
    tt = seq(0, 2*pi, length = npoints)
    xx = center[1] + r*cos(tt)
    yy = center[1] + r*sin(tt)
    return(tibble(x = xx, y = yy))
  }
  corcir = circle(c(0, 0), npoints = 100)
  circle_plot = ggplot() + geom_path(data = corcir, aes(x = x, y = y), colour = "gray65") +
    geom_hline(yintercept = 0, colour = "gray65") +
    geom_vline(xintercept = 0, colour = "gray65") +
    xlim(-1.1, 1.1) + ylim(-1.1, 1.1) +
    geom_segment(data = expr_cor, aes(x = x0, y = y0, xend = MDS1, yend = MDS2),
                 colour = "gray65") +
    geom_text(data = expr_cor, aes(x = MDS1, y = MDS2, label = protein_selection)) +
    labs(x = "MDS1") +
    labs(y = "MDS2") +
    coord_fixed()

  # plot MDS
  mds_plot = ggplot(expr_median, aes_string(x = "MDS1", y = "MDS2",color = color)) +
    geom_point(size = 2) +
    coord_fixed(ratio = explained_var[2] / explained_var[1]) +
    xlab(paste0("MDS1 (",explained_var[1],"%)")) +
    ylab(paste0("MDS2 (",explained_var[2],"%)"))
  if(nchar(sample_label) > 1)
    mds_plot  = mds_plot + geom_label(aes_string(label = sample_label))

  plot_grid(mds_plot, circle_plot, nrow = 1, rel_widths = c(0.6,0.4))
}
