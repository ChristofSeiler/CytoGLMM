#' Check if samples match or paired on condition
#'
#' @param df_samples_subset Data frame or tibble with proteins counts,
#'   cell condition, and group information
#' @param condition The column name of the condition variable
#' @param group The column name of the group variable
#'
#' @return A boolean
#'
is_unpaired <- function(df_samples_subset,
                        condition,
                        group) {

  cell_count <- table(pull(df_samples_subset,group),
                      pull(df_samples_subset,condition))
  unpaired <- TRUE
  if(sum(apply(cell_count,1,min) > 0) == nrow(cell_count))
    unpaired <- FALSE
  unpaired

}
