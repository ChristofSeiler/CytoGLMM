#' Logistic mixture regression
#'
#' @import magrittr
#' @import stringr
#' @import flexmix
#' @import BiocParallel
#' @export
#'
cytoflexmix = function(df_samples_subset,
                       protein_names,
                       condition,
                       cell_n_min = Inf,
                       cell_n_subsample = 0,
                       ks = 1:10,
                       seed = 0xdada,
                       num_cores = 4) {

  # some error checks
  if(cell_n_subsample > cell_n_min) stop("cell_n_subsample is larger than cell_n_min")
  if(sum(str_detect(protein_names,"/")) > 0) stop("protein names cannot contain '/'")
  starts_with_number = sapply(protein_names,
         function(x) str_locate(x, "[0-9]")[1] == 1)
  if(sum(starts_with_number,na.rm = TRUE) > 0) stop("protein names cannot start with numbers")
  if(sum(make.names(protein_names) != protein_names) > 0)
    stop("cleanup your protein names (don't use special characters)")
  set.seed(seed)

  # are the samples paired?
  cell_count = table(df_samples_subset$donor,pull(df_samples_subset,condition))
  unpaired = TRUE
  if(sum(apply(cell_count,1,min) > 0) == nrow(cell_count))
    unpaired = FALSE

  # remove donors with low cell count
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

  # create formulas
  df_samples_subset$donor %<>% as.factor
  df_samples_subset %<>% mutate(
    xtreatment = ifelse(df_samples_subset$treatment ==
                          levels(df_samples_subset$treatment)[1],yes = 0,no = 1)
  )
  #varying_formula = paste0("cbind(xtreatment,1-xtreatment) ~ 1 | donor")
  #fixed_formula = paste("~",paste(protein_names, collapse = " + "))
  varying_formula = paste0("cbind(xtreatment,1-xtreatment) ~ (",
                           paste(protein_names, collapse = " + "),
                           ") | donor")

  # find best number of cluster
  param = MulticoreParam(workers = num_cores,
                         tasks = length(ks),
                         progressbar = TRUE,
                         RNGseed = seed)
  flexmixfits = bplapply(ks,
                         function(k) {
                           stepFlexmix(as.formula(varying_formula),
                                       data = df_samples_subset,
                                       model = FLXMRglm(family = "binomial"),
                                       #model = FLXMRglmfix(family = "binomial",
                                       #                    fixed = as.formula(fixed_formula)),
                                       k = k,
                                       nrep = 5)
                           },
                         BPPARAM = param)

  # return cytoflexmix object
  fit = NULL
  fit$flexmixfits = flexmixfits
  fit$df_samples_subset = df_samples_subset
  fit$protein_names = protein_names
  fit$condition = condition
  fit$cell_n_min = cell_n_min
  fit$cell_n_subsample = cell_n_subsample
  fit$seed = seed
  fit$ks = ks
  fit$num_cores = num_cores
  class(fit) = "cytoflexmix"
  fit

}
