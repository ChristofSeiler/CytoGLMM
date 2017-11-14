#' Logistic regression with cases bootstrap
#'
#' @import magrittr
#' @import stringr
#' @import parallel
#' @import flexmix
#' @import cowplot
#' @export
#'
cytoflexmix = function(df_samples_subset,
                       protein_names,
                       condition,
                       cell_n_min = Inf,
                       cell_n_subsample = 0,
                       ks = 1:4,
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

  # mixtures of logistic regressions
  xlab_str = df_samples_subset %>%
    pull(condition) %>%
    levels %>% paste(collapse = " <-> ")
  df_samples_subset$treatment = ifelse(df_samples_subset$treatment == "control",yes = 0,no = 1)
  pnames = paste0("cbind(",condition,",1-",condition,") ~ ",
                 paste(protein_names,
                 collapse = " + "))

  # find best number of cluster
  #fo = as.numeric(as.factor(df_samples_subset$donor))
  fits = mclapply(ks,FUN = function(k) {
    flexmix(as.formula(pnames),
            data = df_samples_subset,
            k = k,
            model = FLXMRglm(family = "binomial"),
            #model = FLXMRglmnet(family = "binomial",foldid = fo),
            control = list(iter.max = 10))
    },mc.cores = num_cores)
  BICs = sapply(fits,BIC)
  best = which.min(BICs)
  pbic = ggplot(tibble(k = ks,BIC = BICs),aes(k,BIC)) +
    geom_point() +
    geom_line() +
    geom_vline(xintercept = best,color = "red") +
    ggtitle("Model Selection")

  selected_fit = fits[[best]]
  tb = lapply(seq(selected_fit@components),
              function(comp_id) {
                selected_fit@components[[comp_id]][[1]]@parameters[1]$coef %>%
                  data.frame %>%
                  t %>%
                  as.tibble
              }
         ) %>% bind_rows
  tb = tb[,protein_names]
  tb %<>% mutate(comp_id = seq(selected_fit@components))
  tb %<>% gather(protein_name,coeff,-comp_id)
  tb$comp_id %<>% as.factor

  pcoef = ggplot(tb,aes(x = protein_name,y = coeff,color = comp_id)) +
    geom_hline(yintercept = 0, color = "black") +
    geom_point(size = 2) +
    ggtitle("Mixture of Regressions Coefficients") +
    ylab(xlab_str) +
    coord_flip()

  plot_grid(pbic,pcoef,rel_widths = c(0.3,0.7))

}
