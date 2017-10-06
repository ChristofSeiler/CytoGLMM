#' Plot 2d histogram.
#'
#' @import ggplot2
#' @import dplyr
#' @import magrittr
#' @import RColorBrewer
#' @import hexbin
#' @export
#'
plot_donors_2d = function(df_samples,
                          condition,
                          protein1,
                          protein2) {

  # define color scale
  colorscale = scale_fill_gradientn(
    colors = rev(brewer.pal(9, "YlGnBu")),
    values = c(0, exp(seq(-5, 0, length.out = 100)))
  )

  # 2d histogram
  ggplot(df_samples,aes_string(x = protein1,y = protein2)) +
    stat_binhex(binwidth = c(0.2, 0.2)) +
    colorscale +
    coord_fixed() +
    facet_wrap(c(condition,"donor"),ncol = 8)
}
