#' Logistic regression with cases bootstrap
#'
#' @import magrittr
#' @import stringr
#' @import parallel
#' @export
#'
cytoglm = function(df_samples_subset,
                   protein_names,
                   condition,
                   group = "donor",
                   cell_n_min = Inf,
                   cell_n_subsample = 0,
                   num_boot = 100,
                   seed = 0xdada,
                   num_cores = 4) {

  set.seed(seed)

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

  bs = function(seed) {
    set.seed(seed)
    # bootstrap sample
    donor_boot = NULL
    if(unpaired) {
      donor_boot = df_samples_subset %>%
        group_by_(group,condition) %>%
        tally() %>%
        group_by_(condition) %>%
        sample_frac(replace = TRUE) %>%
        ungroup
    } else {
      donor_boot = df_samples_subset %>%
        group_by_(group) %>%
        tally() %>%
        sample_frac(replace = TRUE)
    }
    df_boot = inner_join(donor_boot,
                         df_samples_subset,
                         by = group,
                         suffix = c("",".y")) %>% droplevels

    # logistic regression
    fit_glm = glm(formula = paste(condition,"~",
                                  paste(protein_names,
                                        collapse = " + ")),
                  family = binomial(),
                  data = df_boot)
    tibble(protein_name = protein_names,
           coeff = fit_glm$coefficients[protein_names],
           run = seed)
  }
  tb_coef = mclapply(1:num_boot,bs,mc.cores = num_cores) %>% bind_rows()

  # return cytoglm object
  fit = NULL
  fit$tb_coef = tb_coef
  fit$df_samples_subset = df_samples_subset
  fit$protein_names = protein_names
  fit$condition = condition
  fit$group = group
  fit$cell_n_min = cell_n_min
  fit$cell_n_subsample = cell_n_subsample
  fit$unpaired = unpaired
  fit$num_boot = num_boot
  fit$seed = seed
  fit$num_cores = num_cores
  class(fit) = "cytoglm"
  fit

}
