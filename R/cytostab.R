#' Evaluate parameter stability with respect to gating sheme.
#'
#' @import magrittr
#' @import stringr
#' @import strucchange
#' @import tibble
#' @export
#'
cytostab = function(df_samples_subset,
                    protein_names,
                    condition,
                    group = "donor",
                    cell_n_min = Inf,
                    cell_n_subsample = 0) {

  # some error checks
  cyto_check(cell_n_subsample = cell_n_subsample,
             cell_n_min = cell_n_min,
             protein_names = protein_names)

  # are the samples paired?
  unpaired = is_unpaired(df_samples_subset,
                         condition = condition,
                         group = group)

  # remove donors with low cell count
  df_samples_subset = remove_samples(df_samples_subset,
                                     condition = condition,
                                     group = group,
                                     unpaired = unpaired,
                                     cell_n_min = cell_n_min)

  # subsample cells
  if(cell_n_subsample > 0) {
    df_samples_subset %<>%
      group_by_(group,condition) %>%
      sample_n(size = cell_n_subsample) %>%
      ungroup
  }

  # compute generalized empirical M-fluctuation process
  fluc = gefp(paste(condition,"~",
                    paste(protein_names,
                          collapse = " + ")),
              family = binomial,
              data = df_samples_subset)

  # gather common test statistics about process
  tb = tibble(
    max_abs = apply(fluc$process,
                    MARGIN = 2,
                    function(x) max(abs(x))),
    mean = colMeans(fluc$process),
    range = apply(fluc$process,
                  MARGIN = 2,
                  function(x) diff(range(x))),
    protein_name = names(fluc$process)
  )
  tb

}
