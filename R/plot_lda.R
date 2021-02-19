#' LDA on marker expression
#'
#' @import ggplot2
#' @import dplyr
#' @import magrittr
#' @import tibble
#' @importFrom ggrepel geom_text_repel
#' @importFrom MASS lda
#' @export
#'
#' @param df_samples Data frame or tibble with proteins counts,
#'   cell condition, and group information
#' @param protein_names A vector of column names of protein to use in the analysis
#' @param group The column name of the group variable
#' @param cor_scaling_factor Scaling factor of circle of correlations
#' @param arrow_color Color of correlation circle
#' @param marker_color Colors of marker names
#' @param marker_size Size of markerr names
#' @return \code{\link[ggplot2]{ggplot2}} object
#'
#' @examples
#' set.seed(23)
#' df = generate_data()
#' protein_names = names(df)[3:12]
#' df = dplyr::mutate_at(df, protein_names, function(x) asinh(x/5))
#' df$condition = rep(c("A", "B", "C", "D"), each = length(df$condition)/4)
#' CytoGLMM::plot_lda(df,
#'                    protein_names = protein_names,
#'                    group = "condition",
#'                    cor_scaling_factor = 2)
plot_lda = function(df_samples,
                    protein_names,
                    group,
                    cor_scaling_factor = 1,
                    arrow_color = "black",
                    marker_color = "black",
                    marker_size = 5) {

  if(nlevels(factor(dplyr::pull(df_samples, group))) <= 2)
    stop("need more than 2 levels in group variable")

  formula_str = paste(group, "~", paste(protein_names, collapse = " + "))
  expr_lda = lda(as.formula(formula_str), data = df_samples)
  expr_lda_pred = predict(expr_lda)
  # collect all necessary plotting info
  tb = tibble(type = dplyr::pull(df_samples, group),
              LD1 = expr_lda_pred$x[,1],
              LD2 = expr_lda_pred$x[,2],
              pred = expr_lda_pred$class)
  # base plot
  gglda = ggplot(tb, aes(LD1, LD2, color = type)) +
    geom_density_2d() +
    ggtitle("LDA With Marker Correlations") +
    labs(color = group)
  # make circle of correlation plot
  # correlations between variables and MDS axes
  expr_cor = cor(df_samples %>% dplyr::select(protein_names),
                 tb[,c("LD1","LD2")]) %>% as_tibble
  # scaling factor (otherwise too crowded)
  expr_cor = expr_cor * cor_scaling_factor
  expr_cor %<>% add_column(protein_names)
  # add arrows coordinates
  expr_cor %<>% add_column(x0 = rep(0,nrow(expr_cor)))
  expr_cor %<>% add_column(y0 = rep(0,nrow(expr_cor)))
  # add correlation arrows
  gglda = gglda +
    annotate("segment",
             x = expr_cor$x0, xend = expr_cor$LD1,
             y = expr_cor$y0, yend = expr_cor$LD2,
             colour = arrow_color,
             alpha = 1.0,
             arrow = arrow(type = "open", length = unit(0.03, "npc")))
  ## marker names labels
  gglda + geom_text_repel(data = expr_cor,
                          aes(x = LD1, y = LD2,
                              label = protein_names),
                          size = marker_size,
                          color = marker_color)

}
