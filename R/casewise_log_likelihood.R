#' Get casewise log-likelihood when full data is available
#'
#' @description Get casewise log-likelihood for complete data,
#' useful for WAIC, LOOIC, ...
#' @param object (mbsem) An object of \code{\link{mbsem-class}}
#' returned by \code{\link{minorbsem}}.
#' @param include_residuals (LOGICAL) TRUE: Include minor factor
#' residual covariances in model-implied covariance matrix;
#' FALSE: Exclude them. If TRUE, different
#' models fit to the data will hardly be distinguishable.
#' See details below.
#' @returns matrix (posterior iterations BY sample size)
#' containing log-likelihood
#' @details
#' If comparing two models fit using \code{\link{minorbsem}}, it is
#' reasonable to
#' exclude the residual covariance matrix that captures the influences
#' of minor factors when computing the log-likelihood. If they are not
#' excluded, the log-likelihood will be near identical for models with
#' different structures. This is because the minor factor residual
#' covariances capture the degree of misspecification with the hypothesized
#' structure.
#'
#' The option to set \code{include_residuals = TRUE} is included to allow
#' comparison of models fit with minorbsem and other packages.
#' When the influence of minor factors is non-trivial, one can expect
#' that models fit with minorbsem will have better fit to the data
#' since minorbsem models simultaneously model the degree of
#' model misspecification.
#' @examples
#' # Comparing two models using LOOCV
#' fit_1 <- minorbsem("F1 =~ x1 + x2 + x3
#'                     F2 =~ x4 + x5 + x6
#'                     F3 =~ x7 + x8 + x9", HS)
#' # Compute case wise log-likelihood, exclude minor factor residuals
#' ll_mat_1 <- casewise_log_likelihood(fit_1, include_residuals = FALSE)
#' chain_id <- posterior::as_draws_df(fit_1@stan_fit)$.chain
#' loo_1 <- loo::loo(
#'   ll_mat_1,
#'   r_eff = loo::relative_eff(ll_mat_1, chain_id = chain_id)
#' )
#' print(loo_1)
#' fit_2 <- minorbsem("F1 =~ x1 + x2 + x3 + x9
#'                     F2 =~ x4 + x5 + x6
#'                     F3 =~ x7 + x8 + x9", HS)
#' # Compute case wise log-likelihood, exclude minor factor residuals
#' ll_mat_2 <- casewise_log_likelihood(fit_2, include_residuals = FALSE)
#' loo_2 <- loo::loo(
#'   ll_mat_2,
#'   r_eff = loo::relative_eff(ll_mat_2, chain_id = chain_id)
#' )
#' print(loo_2)
#' # Compare both models
#' print(loo::loo_compare(loo_1, loo_2), simplify = FALSE)
#' @export
casewise_log_likelihood <- function(object, include_residuals = FALSE) {
  # data list must have full data
  data_list <- object@data_list

  if (data_list$has_data != 1) {
    stop("Cannot compute casewise log-likelihood without full data.")
  }

  post_mat <- posterior::as_draws_matrix(object@stan_fit)

  m_vcov_list <- create_model_implied_vcov(
    post_mat, data_list,
    include_residuals = include_residuals
  )

  y_dat <- t(t(data_list$Y) - colMeans(data_list$Y))
  mu <- rep(0, data_list$Ni)
  result <- t(apply(m_vcov_list, 2, function(sigma) {
    m_vcov <- matrix(sigma, nrow = data_list$Ni, ncol = data_list$Ni)
    ll <- mvtnorm::dmvnorm(y_dat, mean = mu, sigma = m_vcov, log = TRUE)
    return(ll)
  }))

  return(result)
}
