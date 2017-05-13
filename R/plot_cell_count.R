#' Plot cell count across donors.
#'
#' @import ggplot2
#' @export
#'
plot_cell_count <- function(df_samples,condition = "condition") {
  ggplot(data = df_samples, aes_string("donor", fill = condition)) +
    geom_bar(position="dodge") +
    ggtitle("Cell Count")
}
