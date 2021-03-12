#' Check if input to cytoxxx function have errors
#'
#' @param cell_n_subsample Subsample samples to have this maximum cell count
#' @param cell_n_min A vector of column names of protein to use in the analysis
#' @param protein_names A vector of column names of protein to use in the
#'   analysis
#' @return NULL.
#'
cyto_check = function(cell_n_subsample,
                      cell_n_min,
                      protein_names) {

  if(cell_n_subsample > cell_n_min)
    stop("cell_n_subsample is larger than cell_n_min")
  if(sum(str_detect(protein_names,"/")) > 0)
    stop("protein names cannot contain '/'")
  starts_with_number <- vapply(
    protein_names,
    function(x) str_locate(x, "[0-9]")[1] == 1, logical(1)
    )
  if(sum(starts_with_number,na.rm = TRUE) > 0)
    stop("protein names cannot start with numbers")
  if(sum(make.names(protein_names) != protein_names) > 0)
    stop("cleanup your protein names (don't use special characters)")

}
