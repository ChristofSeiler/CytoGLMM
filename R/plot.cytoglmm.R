#' Plot fixded coefficients of random effects model
#'
#' @aliases plot.cytoglmm
#' @method plot cytoglmm
#'
#' @import ggplot2
#' @import tibble
#' @import magrittr
#' @import dplyr
#' @importFrom methods is
#' @importFrom rlang .data
#' @export
#'
#' @param x A \code{cytoglmm} class
#' @param order Order the markers according to the mangintute of the coefficients
#' @param separate create two separate \code{\link[ggplot2]{ggplot2}} objects
#' @param ... Other parameters
#' @return \code{\link[ggplot2]{ggplot2}} object
#'
#' @examples
#' set.seed(23)
#' df <- generate_data()
#' protein_names <- names(df)[3:12]
#' df <- dplyr::mutate_at(df, protein_names, function(x) asinh(x/5))
#' glmm_fit <- CytoGLMM::cytoglmm(df,
#'                                protein_names = protein_names,
#'                                condition = "condition",
#'                                group = "donor")
#' plot(glmm_fit)
plot.cytoglmm <- function(x, order = FALSE, separate = FALSE, ...) {

  if(!is(x, "cytoglmm"))
    stop("Input needs to be a cytoglmm object computed by cytoglmm function.")
  if(!is.factor(x$glmmfit$y))
    stop("Currently only plotting results for logistic regression.")

  summ <- summary(x$glmmfit)

  # random effects
  stdev <- sqrt(diag(summ$varcor[[1]])[-1])
  tb_random <- tibble(protein_name = names(stdev),
                     stdev = stdev)
  tb_random <- tb_random[tb_random$protein_name %in% x$protein_names,]

  # fixed effects
  alpha <- 0.05
  ci <- qnorm(1-alpha/2)
  tb_coeff <- summ$coefficient
  tb_coeff <- tb_coeff[-1,]
  tb_coeff %<>%
    as.data.frame %>%
    rownames_to_column(var = "protein_name") %>%
    as.tibble
  tb_coeff %<>%
    mutate(high = tb_coeff$Estimate+ci*tb_coeff$`Std. Error`,
           low = tb_coeff$Estimate-ci*tb_coeff$`Std. Error`)
  tb_coeff <- tb_coeff[tb_coeff$protein_name %in% x$protein_names,]

  # order proteins according to coefficients
  if(order) {
    ind <- sort.int(tb_coeff$Estimate,index.return=TRUE)$ix
    reordered_names <- tb_coeff$protein_name[ind]
    tb_coeff$protein_name %<>% factor(levels = reordered_names)
    tb_random$protein_name %<>% factor(levels = reordered_names)
  }

  # plotting
  prandom <- ggplot(tb_random, aes(x = .data$stdev, y = .data$protein_name)) +
    geom_point(size = 2) +
    ggtitle("Random Effects") +
    xlab("standard deviation") +
    theme(axis.title.y = element_blank())
  xlab_str <- x$df_samples_subset %>%
    pull(x$condition) %>%
    levels %>%
    paste(collapse = " <-> ")
  pcoef <- ggplot(tb_coeff, aes(x = .data$Estimate, y = .data$protein_name)) +
    geom_vline(xintercept = 0,color = "red") +
    geom_point(size = 2) +
    geom_errorbarh(aes(xmin = .data$low, xmax = .data$high)) +
    ggtitle("cytoglmm") +
    xlab(xlab_str) +
    theme(axis.title.y = element_blank())

  if(separate) {
    return(list(prandom=prandom,pcoef=pcoef))
  } else {
    pcoef
  }

}
