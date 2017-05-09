#' Plot one density plot faceted by donors.
#'
#' @import ggplot2
#' @export
#'
plot_donors <- function(df_samples,
                       protein_name,
                       histogram = TRUE) {
  if(!histogram) {
    ggplot(df_samples, aes_string(x = protein_name, color = "condition")) +
      geom_density() +
      facet_wrap(~ donor) +
      ggtitle(protein_name) +
      xlab("arcsinh transformed counts")
  } else {
    ggplot(df_samples, aes_string(x = protein_name, fill = "condition")) +
      geom_histogram(binwidth=.2, alpha=.5, position="identity") +
      facet_wrap(~ donor) +
      ggtitle(protein_name) +
      xlab("arcsinh transformed counts")
  }
}
