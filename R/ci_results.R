#' Conditional indepepence estimates for path analysis
#'
#' @description Obtain conditional indepepence estimates for path analysis.
#' Estimates are standardized.
#' @param object (mbsem) An object of \code{\link{mbsem-class}}
#' returned by \code{\link{minorbsem}}.
#' @param interval Confidence interval to select
#' @param summarize (LOGICAL)
#' If TRUE (default): Return posterior summary as data.frame;
#' If FALSE: Return posterior draws as data.frame.
#' @returns Table of parameters related to conditional independence.
#' @examples
#' \dontrun{
#' fit <- minorbpa("x3 ~ x1 + x2\n x4 ~ x3 + x1", HS)
#' ci_table <- ci_results(fit)
#' }
#' @export
ci_results <- function(object, interval = .9, summarize = TRUE) {
  stopifnot(inherits(object, "mbsem"))

  if (object@data_list$method >= 90) {
    stop(paste0(
      "There are no residuals to plot when ",
      "method == \"none\", \"WB\", \"WB-cond\", \"WW\"."
    ))
  }

  if (object@data_list$pa_indicator != 1) {
    stop(paste0(
      "This method only works for objects produced by ",
      "minorbpa()."
    ))
  }

  if (sum(object@data_list$cond_ind_mat) == 0) {
    message(
      "All possible associations are modelled."
    )
    return(NULL)
  }

  ind_names <- rownames(object@data_list$loading_pattern)
  ci_details <- object@data_list$cond_ind_details
  x_var <- ind_names[ci_details[, 1]]
  y_var <- ind_names[ci_details[, 2]]
  z_var <- vector("character", nrow(ci_details))
  for (i in seq_len(nrow(ci_details))) {
    z_var[i] <- paste0(
      ind_names[as.logical(ci_details[i, 3:ncol(ci_details)])],
      collapse = ", "
    )
  }
  var_s <- paste0(x_var, " _||_ ", y_var, " | ", z_var)

  resid_idxs <- paste0("Resid[", ci_details[, 1], ",", ci_details[, 2], "]")

  if (isFALSE(summarize)) {
    result <- posterior::as_draws_df(object@stan_fit$draws(resid_idxs))
    result <- result[
      ,
      c(resid_idxs, ".chain", ".iteration", ".draw")
    ]
    colnames(result)[seq_along(var_s)] <- var_s
  } else {
    resid_idxs_i <- as.integer(factor(resid_idxs, unique(resid_idxs)))
    post_sum <- mbsem_post_sum(object@stan_fit, resid_idxs, interval = interval)
    result <- cbind(relation = var_s, post_sum[resid_idxs_i, -1])
    rownames(result) <- NULL
  }

  return(result)
}
