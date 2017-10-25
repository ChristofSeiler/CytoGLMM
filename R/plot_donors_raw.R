#' Plot row counts faceted by donors
#'
#' @import ggplot2
#' @import magrittr
#' @import dplyr
#' @export
#'
plot_donors_raw = function(df_samples_subset,
                           condition,
                           protein_name,
                           min = 0,
                           max = quantile(df_samples_subset %>% pull(protein_name),
                                          probs = 0.99)) {

  donors = df_samples_subset %>%
    group_by_("donor",condition) %>%
    tally()
  donors

  tb_raw = df_samples_subset %>%
    dplyr::rename_("condition__" = condition) %>%
    group_by(condition__,donor) %>%
    count_(protein_name) %>%
    dplyr::filter_(paste(protein_name,">",min,"&",protein_name,"<",max))
  tb_raw %<>% left_join(donors,by = "donor")

  ggplot(tb_raw,
         aes_string(x = protein_name,
                    y = "n.x/n.y",
                    color = "condition__",
                    group = "donor")) +
    geom_line(alpha = 0.5) +
    geom_point(alpha = 0.5) +
    xlab("protein count") +
    ylab("cell count / total donor cell count") +
    ggtitle(protein_name) +
    theme(legend.title=element_blank())
}
