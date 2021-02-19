#' Remove samples based on low cell counts
#'
#' @param df_samples_subset Data frame or tibble with proteins counts,
#'   cell condition, and group information
#' @param condition The column name of the condition variable
#' @param group The column name of the group variable
#' @param unpaired true if unpaired samples were provided as input
#' @param cell_n_min Remove samples that are below this cell counts threshold
#' @return NULL.
#'
remove_samples = function(df_samples_subset,
                          condition,
                          group,
                          unpaired,
                          cell_n_min) {

  cell_count = table(pull(df_samples_subset,group),
                     pull(df_samples_subset,condition))
  include = NULL
  if(cell_n_min < Inf) {
    if(unpaired) {
      exclude = which(rowSums(cell_count) < cell_n_min)
    } else {
      exclude = which(apply(cell_count,1,min) < cell_n_min)
    }
    include = rownames(cell_count)[!rownames(cell_count) %in% names(exclude)]
  } else {
    include = rownames(cell_count)
  }
  df_samples_subset[pull(df_samples_subset,group) %in% include,] %>%
    droplevels()

}
