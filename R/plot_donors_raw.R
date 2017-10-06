#' Plot row counts faceted by donors
#'
#' @import ggplot2
#' @import magrittr
#' @import dplyr
#' @export
#'
plot_donors_raw = function(df_samples_subset,
                           protein_name,
                           min = 0,
                           max = quantile(df_samples_subset %>% pull(protein_name),
                                          probs = 0.99)) {
  ggplot(df_samples_subset %>%
           group_by(response_type,donor) %>%
           count_(protein_name) %>%
           dplyr::filter_(paste(protein_name,">",min,"&",protein_name,"<",max)),
         aes_string(x = protein_name,
                    y = "n",
                    color = "response_type",
                    group = "donor")) +
    geom_line(alpha = 0.5) +
    xlab("protein count") +
    ylab("cell count") +
    ggtitle(protein_name)
}
