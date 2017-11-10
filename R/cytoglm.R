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
                   cell_n_min = Inf,
                   cell_n_subsample = 0,
                   unpaired = TRUE,
                   num_boot = 100,
                   seed = 0xdada,
                   num_cores = 4) {

  if(cell_n_subsample > cell_n_min) stop("cell_n_subsample is larger than cell_n_min")
  set.seed(seed)

  # remove donors with low cell count
  include = NULL
  if(cell_n_min < Inf) {
    cell_count = table(df_samples_subset$donor,pull(df_samples_subset,condition))
    exclude = which(rowSums(cell_count) < cell_n_min)
    cell_count[exclude,]
    include = rownames(cell_count)[rownames(cell_count) %nin% names(exclude)]
  } else {
    include = rownames(cell_count)
  }
  df_samples_subset %<>%
    dplyr::filter(donor %in% include) %>%
    droplevels()

  # subsample cells
  if(cell_n_subsample > 0) {
    df_samples_subset %<>%
      group_by_("donor",condition) %>%
      sample_n(size = cell_n_subsample) %>%
      ungroup
  }

  bs = function(seed) {
    set.seed(seed)
    # bootstrap sample
    donor_boot = NULL
    if(unpaired) {
      donor_boot = df_samples_subset %>%
        group_by_("donor",condition) %>%
        tally() %>%
        group_by_(condition) %>%
        sample_frac(replace = TRUE) %>%
        ungroup
    } else {
      donor_boot = df_samples_subset %>%
        group_by(donor) %>%
        tally() %>%
        sample_frac(replace = TRUE)
    }
    df_boot = inner_join(donor_boot,
                         df_samples_subset,
                         by = "donor",
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

  # return cytomlogit object
  fit = NULL
  fit$tb_coef = tb_coef
  fit$df_samples_subset = df_samples_subset
  fit$protein_names = protein_names
  fit$condition = condition
  fit$cell_n_min = cell_n_min
  fit$cell_n_subsample = cell_n_subsample
  fit$unpaired = unpaired
  fit$num_boot = num_boot
  fit$seed = seed
  fit$num_cores = num_cores
  class(fit) = "cytoglm"
  fit

}
