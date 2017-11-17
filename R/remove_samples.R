#' Remove samples based on low cell counts
#'
remove_samples = function(df_samples_subset,
                          condition,
                          group,
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
    include = rownames(cell_count)[rownames(cell_count) %nin% names(exclude)]
  } else {
    include = rownames(cell_count)
  }
  df_samples_subset[pull(df_samples_subset,group) %in% include,] %>%
    droplevels()

}
