#' @importFrom methods show
methods::setMethod(
  "show",
  "mbsem",
  function(object) {
    pretty_print_summary(object)
  }
)

#' @export
methods::setGeneric("plot")
methods::setMethod(
  "plot",
  "mbsem",
  function(x, type = "hist", param_type = c("rm", "co", "lo", "fc", "rsq")) {
    stopifnot(type %in% c("hist", "trace"))

    if (type == "hist") {
      parameter_hist(x, param_type = param_type)
    } else if (type == "trace") {
      parameter_trace(x, param_type = param_type)
    }
  }
)

methods::setMethod(
  "fitted",
  "mbsem",
  function(object, include_residuals = TRUE) {
    post_mat <- posterior::as_draws_matrix(object@stan_fit)

    m_vcov_mat <- create_mi_vcov_ll(
      mat = post_mat,
      data_list = object@data_list,
      include_residuals = include_residuals,
      return_ll = FALSE
    )

    m_vcov_mat <- t(m_vcov_mat)

    n_iter <- nrow(m_vcov_mat)

    m_vcov_mat <- matrix(m_vcov_mat, nrow = n_iter)

    return(m_vcov_mat)
  }
)

methods::setMethod(
  "residuals",
  "mbsem",
  function(object, standardized = TRUE) {
    if (object@data_list$method >= 90) {
      warning(paste0(
        "There are no residuals to plot when ",
        "method == \"none\", \"WB\", \"WB-cond\"."
      ))
    }

    resid_mat <- posterior::subset_draws(
      posterior::as_draws_matrix(object@stan_fit),
      variable = "Resid"
    )

    n_iter <- nrow(resid_mat)
    resid_mat <- matrix(resid_mat, nrow = n_iter)

    if (isFALSE(standardized)) {
      k <- object@data_list$Ni
      m_vcov_mat <- fitted(object, include_residuals = FALSE)
      resid_mat <- t(sapply(1:n_iter, function(i) {
        resids <- matrix(resid_mat[i, ], nrow = k)
        m_vcov <- matrix(m_vcov_mat[i, ], nrow = k)
        sd_mat <- diag(sqrt(diag(m_vcov)))
        resids <- sd_mat %*% resids %*% sd_mat
      }))
    }

    return(resid_mat)
  }
)

methods::setMethod(
  "logLik",
  "mbsem",
  function(object, include_residuals = FALSE) {
    casewise_log_likelihood(object, include_residuals)
  }
)

methods::setValidity("mbsempriors", function(object) {
  if (any(c(object@sl_par, object@rs_par, object@sc_par) <= 0)) {
    "@sl_par, @rs_par and @sc_par must all be greater than 0"
  } else if (any(c(object@lkj_shape, object@rc_par, object@fc_par) < 1)) {
    "@lkj_shape, @rc_par and @fc_par must all be at least 1"
  } else {
    TRUE
  }
})
