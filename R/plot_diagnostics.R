#' Trace plots for HMC.
#'
#' @import rstan
#' @import reshape2
#' @import ggplot2
#' @import cowplot
#' @export
#'
plot_diagnostics <- function(fit,par_name,num_par = 8) {
  if (!inherits(fit, "stanfit"))
    stop("Not a stanfit object.")
  if (fit@mode != 0)
    stop("Stan model does not contain posterior draws.")

  # keep only num_par paramters to avoid overloaded plots
  param = rstan::extract(fit,
                         pars = par_name,
                         permuted = FALSE,
                         inc_warmup = TRUE)
  par_subset_ids = sample(dim(param)[3],size = num_par) %>% sort
  param_long_during = melt(param[1:fit@stan_args[[1]]$warmup,,par_subset_ids],
                      varnames = c("iteration","parameter"))
  param_long_after = melt(param[(fit@stan_args[[1]]$warmup+1):dim(param)[1],,par_subset_ids],
                    varnames = c("iteration","parameter"))

  p1 = ggplot(param_long_during, aes(x = iteration, y = value, color = parameter)) +
    geom_line() +
    ggtitle("Traceplot During Warmup") +
    theme(legend.position="none",
          axis.title.y=element_blank())
  p2 = ggplot(param_long_after, aes(x = iteration, y = value, color = parameter)) +
    geom_line() +
    ggtitle("Traceplot After Warmup") +
    theme(legend.position="none",
          axis.title.y=element_blank())

  legend_b = get_legend(p1 + theme(legend.position="bottom"))
  prow = plot_grid(p1, p2,align = "h")
  plot_grid(prow, legend_b, ncol = 1, rel_heights = c(1, .1))
}
