#' @importFrom methods show
methods::setMethod(
  "show",
  "mbsem",
  function(object) {
    pretty_print_summary(object)
  }
)

methods::setMethod(
  "fitted",
  "mbsem",
  function(object) {
    omega_mat <- posterior::as_draws_matrix(object@stan_fit$draws("Omega"))

    return(omega_mat)
  }
)

methods::setMethod(
  "residuals",
  "mbsem",
  function(object, standardized = TRUE) {
    data_list <- object@data_list

    if (data_list$method >= 90) {
      warning(paste0(
        "There are no residuals to plot when ",
        "method == \"none\", \"WB\", \"WB-cond\", \"WW\"."
      ))
    }

    resid_mat <- posterior::as_draws_matrix(object@stan_fit$draws("Resid"))

    n_iter <- nrow(resid_mat)
    resid_mat <- matrix(resid_mat, nrow = n_iter)

    if (isFALSE(standardized)) {
      k <- data_list$Ni
      m_vcov_mat <- fitted(object)
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
  function(object) {
    data_list <- object@data_list

    if (data_list$ret_ll == 0) {
      err_msg <- paste0(
        "No log-likelihood returned. ",
        "You must set \"compute_ll = TRUE\" when calling ?minorbsem. ",
        "If you already did that, see the compute_ll parameter ",
        "description under ?minorbsem."
      )
      stop(err_msg)
    }
    posterior::as_draws_matrix(object@stan_fit$draws("log_lik"))
  }
)

methods::setValidity("mbsempriors", function(object) {
  if (any(
    c(
      object@sl_par, object@rs_par, object@sc_par
    ) <= 0
  )) {
    paste0(
      "@sl_par, @rs_par, and @sc_par ",
      "must all be greater than 0"
    )
  } else if (any(c(object@lkj_shape, object@rc_par, object@fc_par) < 1)) {
    "@lkj_shape, @rc_par and @fc_par must all be at least 1"
  } else {
    TRUE
  }
})
