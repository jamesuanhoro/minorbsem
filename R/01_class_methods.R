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
