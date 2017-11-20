#' Plot one density plot faceted by donors.
#'
#' @import ggplot2
#' @export
#'
plot_donors <- function(df_samples,
                       protein_name,
                       condition = "condition",
                       group = "donor",
                       density = TRUE) {
  if(density) {
    ggplot(df_samples, aes_string(x = protein_name,y = "..scaled..",color = condition)) +
      geom_density() +
      facet_wrap(group) +
      ggtitle(protein_name) +
      xlab("arcsinh transformed counts")
  } else {
    ggplot(df_samples, aes_string(x = protein_name,fill = condition)) +
      geom_histogram(binwidth=.2, alpha=.5, position="identity") +
      facet_wrap(group) +
      ggtitle(protein_name) +
      xlab("arcsinh transformed counts")
  }
}
