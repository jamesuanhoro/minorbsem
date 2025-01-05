#' Fit Bayesian path analysis models
#'
#' @description Fit Bayesian path anlaysis models with tests of
#' conditional independence.
#' @inheritParams minorbsem
#' @returns An object of \code{\link{mbsem-class}}
#' @export
minorbpa <- function(
    model = NULL,
    data = NULL,
    sample_cov = NULL,
    sample_nobs = NULL,
    data_list = NULL,
    method = "normal",
    orthogonal = FALSE,
    correlation = FALSE,
    centered = TRUE,
    seed = 12345,
    warmup = 1000,
    sampling = 1000,
    refresh = (warmup + sampling) / 10,
    adapt_delta = .9,
    max_treedepth = 10,
    chains = 3,
    ncores = max(parallel::detectCores() - 2, 1),
    priors = new_mbsempriors(),
    show = TRUE,
    show_messages = TRUE,
    compute_ll = FALSE,
    acov_mat = NULL,
    ret_data_list = FALSE) {
  message("Processing user input ...")

  if (is.null(data_list)) {
    data_list <- user_input_process(
      model, data, sample_cov, sample_nobs, method,
      orthogonal, simple_struc = TRUE, correlation, centered, priors,
      compute_ll, acov_mat, pa = TRUE
    )
  }

  if (isTRUE(ret_data_list)) {
    return(data_list)
  }

  message("User input fully processed :)\n Now to modeling.")

  stan_fit <- target_fitter(
    data_list, seed, warmup, sampling, refresh,
    adapt_delta, max_treedepth, chains, ncores, show_messages, pa = TRUE
  )

  mbsem_results <- clean_up_stan_fit(
    stan_fit = stan_fit, data_list = data_list, priors = priors
  )
  if (show) {
    show(mbsem_results)
  }

  return(mbsem_results)
}
