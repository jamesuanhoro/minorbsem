#' Clean up Stan fit helper function
#' @description A function that cleans up model returned by Stan
#' @param stan_fit Stan fit
#' @param data_list Data list object passed to Stan
#' @param interval (real in (0, 1)) Credible interval to return.
#' @param priors An object of \code{\link{mbsempriors-class}}.
#' See \code{\link{new_mbsempriors}} for more information.
#' @returns An object of \code{\link{mbsem-class}}
#' @keywords internal
clean_up_stan_fit <- function(
    stan_fit,
    data_list,
    interval = .9,
    priors) {
  major_parameters <- create_major_params(
    stan_fit = stan_fit, data_list = data_list, interval = interval
  )

  minor_factor_matrix <- mbsem_post_sum(
    stan_fit = stan_fit, variable = "Resid", interval = interval
  )

  mbsem_result <- new_mbsem()
  mbsem_result <- methods::initialize(
    mbsem_result,
    major_parameters = major_parameters,
    minor_factor_matrix = minor_factor_matrix,
    data_list = data_list,
    priors = priors,
    stan_fit = stan_fit,
    version = minorbsem_version()
  )

  return(mbsem_result)
}
