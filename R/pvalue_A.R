#' Extact and plot noise term
#'
#' @import ggplot2
#' @import tibble
#' @import magrittr
#' @import dplyr
#' @import cowplot
#' @export
#'
pvalue_A = function(dm_model_list,
                    protein_names = protein_names) {

  # some jobs may fail (because of computing cluster instabilities)
  if(length(dm_model_list) == 0)
    stop("no jobs completed successfully")

  tb_A = extract_A(dm_model_list,protein_names = protein_names)
  obsv = tb_A %>% dplyr::filter(run == 1)
  boot = tb_A %>% dplyr::filter(run > 1)
  tb_A = lapply(unique(boot$run),function(i)
    boot %>%
      dplyr::filter(run == i) %>%
      mutate( coeff_translated = (coeff-obsv$coeff) ) %>%
      add_column( coeff_obsv = obsv$coeff )
    ) %>% bind_rows

  tb_A_summary = tb_A %>%
    group_by(protein_name) %>%
    summarize(
      pvalue = mean(abs(coeff_obsv) <= abs(coeff_translated))
    ) %>%
    mutate(pvalue = ifelse(pvalue > 0,pvalue,1/length(unique(boot$run))) ) %>%
    mutate(pvalue_adj = p.adjust(pvalue,method = "BH")) %>%
    mutate(coeff = obsv$coeff)
  tb_A_summary %>%
    arrange(pvalue)

}
