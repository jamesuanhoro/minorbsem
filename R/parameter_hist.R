#' Histogram of parameter posterior distribution
#'
#' @description Produce histograms of parameter posterior distribution,
#' option to limit plot by type of parameter.
#' @param object (mbsem) An object of \code{\link{mbsem-class}}
#' returned by \code{\link{minorbsem}}.
#' @param subset (character) Subset of parameters:
#' NULL (Default) showing all estimated parameters;
#' Any other response will be used as regular expressions to subset
#' the parameters. It can be loading names or types of parameters.
#' @param ... additional arguments to relevant bayesplot function
#' @returns bayesplot object
#' @examples
#' \dontrun{
#' fit <- minorbsem("F1 =~ x1 + x2 + x3
#'                   F2 =~ x4 + x5 + x6
#'                   F3 =~ x7 + x8 + x9", HS)
#' parameter_hist(fit)
#' }
#' @export
parameter_hist <- function(object, subset = NULL, ...) {
  stopifnot(inherits(object, "mbsem"))

  param_list <- get_param_plot_list(object@data_list)
  if (!is.null(subset)) {
    param_names <- names(param_list)
    param_list <- param_list[grep(subset, param_names, ignore.case = TRUE)]
  }
  pd_mat <- posterior::as_draws_matrix(object@stan_fit$draws(param_list))
  colnames(pd_mat)[seq_along(param_list)] <- names(param_list)

  result <- bayesplot::mcmc_hist(pd_mat, ...)
  return(result)
}
