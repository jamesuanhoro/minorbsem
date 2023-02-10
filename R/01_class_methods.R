#' @importFrom methods show
methods::setMethod(
  "show",
  "mbsem_object",
  function(object) {
    pretty_print_summary(object)
  }
)

methods::setMethod(
  "hist",
  "mbsem_object",
  function(x, param_type = c("rm", "co", "lo", "fc", "fv")) {
    parameter_hist(x, param_type = param_type)
  }
)

methods::setMethod(
  "residuals",
  "mbsem_object",
  function(object, method = "matrix") {
    stopifnot(method %in% c("table", "matrix", "range"))

    if (method == "table") {
      object@minor_factor_matrix
    } else {
      plot_residuals(object, method)
    }
  }
)

methods::setMethod(
  "logLik",
  "mbsem_object",
  function(object, include_residuals = FALSE) {
    casewise_log_likelihood(object, include_residuals)
  }
)
