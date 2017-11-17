#' Check if samples match or paired on condition
#'
is_unpaired = function(df_samples_subset,
                       condition,
                       group) {

  cell_count = table(pull(df_samples_subset,group),
                     pull(df_samples_subset,condition))
  unpaired = TRUE
  if(sum(apply(cell_count,1,min) > 0) == nrow(cell_count))
    unpaired = FALSE
  unpaired

}
