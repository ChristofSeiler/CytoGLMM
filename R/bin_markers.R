#' Transform marker expression to categorical variable with fixed number of levels (bins).
#'
#' @export
#'
bin_markers = function(marker_expr,num_bins = 8) {
  count_range = range(marker_expr)
  bin_breaks = seq(count_range[1],
                   count_range[2],
                   diff(count_range)/num_bins)
  bin_breaks[1] = -Inf
  bin_breaks[length(bin_breaks)] = Inf
  cut(marker_expr,breaks = bin_breaks,labels = 1:num_bins)
}
