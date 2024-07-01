#' Compute residual network
#'
#' @description Interpret the error correlations as a residual network model.
#' @param object (mbsem) An object of \code{\link{mbsem-class}}
#' returned by \code{\link{minorbsem}}.
#' @returns data.frame
#' @examples
#' \dontrun{
#' fit <- minorbsem("F1 =~ x1 + x2 + x3
#'                   F2 =~ x4 + x5 + x6
#'                   F3 =~ x7 + x8 + x9", HS)
#' res_net <- compute_residual_network(fit)
#' p_corr_df <- posterior::summarise_draws(res_net)
#' p_corr_mat <- matrix(p_corr_df$mean, fit@data_list$Ni)
#' p_corr_mat
#' qgraph::qgraph(p_corr_mat, threshold = .05)
#' }
#' @export
compute_residual_network <- function(object) {
  stopifnot(inherits(object, "mbsem"))

  if (object@data_list$method >= 90) {
    stop(paste0(
      "There are no residuals to plot when ",
      "method == \"none\", \"WB\", \"WB-cond\", \"WW\"."
    ))
  }

  resids <- posterior::as_draws_df(object@stan_fit$draws("Resid"))
  len <- ncol(resids) - 3 # exclude special variables for draws
  p <- sqrt(len)
  result <- t(apply(resids, 1, function(x) {
    mat <- matrix(x[1:len], p, p)
    diag(mat) <- 1
    inv_mat <- -1 * stats::cov2cor(MASS::ginv(mat))
    diag(inv_mat) <- 1
    inv_vec <- c(as.numeric(inv_mat), x[(len + 1):(len + 3)])
    return(inv_vec)
  }))

  colnames(result) <- gsub("Resid", "p_corr", colnames(resids))
  result <- posterior::as_draws_df(as.data.frame(result))

  return(result)
}
