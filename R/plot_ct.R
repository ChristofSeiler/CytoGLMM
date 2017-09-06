#' Visualize three-way contigency tables.
#'
#' @import ggplot2
#' @import dplyr
#' @import magrittr
#' @export
#'
plot_ct = function(df_samples,
                   sample_info_names,
                   protein_names,
                   fac) {
  tb = df_samples %>%
    select(c(sample_info_names,protein_names)) %>%
    mutate_at(protein_names,as.numeric)
  tb_long = melt(tb,id.vars = sample_info_names,variable.name = "marker")
  con_tb = tb_long %>%
    group_by(.dots = c("value",sample_info_names)) %>%
    count(marker)
  ggplot(con_tb,aes_string(y = fac,x = "marker")) +
    geom_tile(aes(fill = log(n)),colour = "white") +
    scale_fill_gradient(low = "white", high = "steelblue") +
    facet_wrap(~value,ncol = 1) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
}
