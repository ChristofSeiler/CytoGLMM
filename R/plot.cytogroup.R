#' Plot fixded coefficients of group-specific fixed effects model
#'
#' @aliases plot.cytogroup
#' @method plot cytogroup
#'
#' @import ggplot2
#' @import tibble
#' @import magrittr
#' @import dplyr
#' @importFrom stats qnorm
#' @importFrom methods is
#' @importFrom cowplot plot_grid
#' @importFrom rlang .data
#' @export
#'
#' @param x A \code{cytoglmm} class
#' @param order Order the markers according to the mangintute of the
#'   coefficients
#' @param separate create two separate \code{\link[ggplot2]{ggplot2}} objects
#' @param ... Other parameters
#' @return \code{\link[ggplot2]{ggplot2}} object
#'
#' @examples
#' set.seed(23)
#' df <- generate_data()
#' protein_names <- names(df)[3:12]
#' df <- dplyr::mutate_at(df, protein_names, function(x) asinh(x/5))
#' group_fit <- CytoGLMM::cytogroup(df,
#'                                  protein_names = protein_names,
#'                                  condition = "condition",
#'                                  group = "donor")
#' plot(group_fit)
plot.cytogroup <- function(x, order = FALSE, separate = FALSE, ...) {

  if(!is(x, "cytogroup"))
    stop("Input needs to be a cytogroup object.")

  summ <- summary(x$groupfit)
  xlab_str <- x$df_samples_subset %>%
    pull(x$condition) %>%
    levels %>%
    paste(collapse = " <-> ")

  # donor-specific effects
  tb_coeff_donor <- tibble(names = names(x$groupfit$coefficients),
                           value = x$groupfit$coefficients)
  tb_coeff_donor %<>%
    mutate(protein_name = vapply(str_split(tb_coeff_donor$names,":"),
                                 function(x) x[1], character(1)),
           donor = vapply(str_split(tb_coeff_donor$names,":"),
                          function(x) x[2], character(1)))
  tb_coeff_donor %<>% dplyr::filter(!is.na(.data$donor))

  # fixed effects
  alpha <- 0.05
  ci <- qnorm(1-alpha/2)
  tb_coeff <- summ$coefficient
  tb_coeff <- tb_coeff[x$protein_names,]
  tb_coeff %<>%
    as.data.frame %>%
    rownames_to_column(var = "protein_name") %>%
    as_tibble
  tb_coeff$Estimate %<>% as.character %>% as.numeric
  tb_coeff$`Std. Error` %<>% as.character %>% as.numeric
  tb_coeff$`z value` %<>% as.character %>% as.numeric
  tb_coeff$`Pr(>|z|)` %<>% as.character %>% as.numeric
  tb_coeff %<>%
    mutate(high = tb_coeff$Estimate+ci*tb_coeff$`Std. Error`,
           low = tb_coeff$Estimate-ci*tb_coeff$`Std. Error`)

  # order proteins according to coefficients
  if(order) {
    ind <- sort.int(tb_coeff$Estimate,index.return=TRUE)$ix
    reordered_names <- tb_coeff$protein_name[ind]
    tb_coeff$protein_name %<>% factor(levels = reordered_names)
    tb_coeff_donor$protein_name %<>% factor(levels = reordered_names)
  }

  # plotting
  pdonor <- ggplot(tb_coeff_donor, aes(x = .data$value,
                                       y = .data$protein_name,
                                       group = .data$donor,
                                       color = .data$protein_name)) +
    geom_vline(xintercept = 0, color = "red") +
    geom_jitter(size = 1, width = 0.001, alpha = 0.5) +
    ggtitle("Donor-Specific Effects") +
    xlab(xlab_str) +
    theme(legend.position="none",
          axis.title.y = element_blank())
  pcoef <- ggplot(tb_coeff, aes(x = .data$Estimate, y = .data$protein_name)) +
    geom_vline(xintercept = 0,color = "red") +
    geom_point(size = 2) +
    geom_errorbarh(aes(xmin = .data$low, xmax = .data$high)) +
    ggtitle("Fixed Effects") +
    xlab(xlab_str) +
    theme(axis.title.y = element_blank())

  if(separate) {
    return(list(pdonor=pdonor,pcoef=pcoef))
  } else {
    plot_grid(pdonor,pcoef)
  }

}
