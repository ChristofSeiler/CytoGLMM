#' Plot all markers for one donor.
#'
#' @import ggplot2
#' @import reshape2
#' @export
#'
plot_markers <- function(df_samples,
                        donor_id,
                        protein_names,
                        condition = "condition",
                        density = TRUE) {
  
  df_samples_subset = subset(df_samples,donor == donor_id) %>% 
    select(c(protein_names,condition)) %>% 
    melt(id = condition)
  
  if(density) {
    ggplot(df_samples_subset,aes_string(x = "value",
                                             y = "..scaled..",
                                             color = condition)) +
      geom_density() +
      facet_wrap(~ variable,ncol = 8) +
      ggtitle(donor_id) +
      xlab("arcsinh transformed counts")
  } else {
    ggplot(df_samples_subset,aes_string(x = "value",
                                             fill = condition)) +
      geom_histogram(binwidth=.2, alpha=.5, position="identity") +
      facet_wrap(~ variable,ncol = 8) +
      ggtitle(donor_id) +
      xlab("arcsinh transformed counts")
  }
}
