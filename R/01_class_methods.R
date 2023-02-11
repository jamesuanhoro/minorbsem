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
  "residuals",
  "mbsem",
  function(object, type = "matrix") {
    stopifnot(type %in% c("table", "matrix", "range"))

    if (type == "table") {
      object@minor_factor_matrix
    } else {
      plot_residuals(object, type)
    }
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
