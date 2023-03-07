#' Fit random-effects Bayesian meta-analytic CFAs with minor factors assumed.
#'
#' @description A function to fit random-effects Bayesian meta-analytic CFAs
#' with minor factors assumed \insertCite{uanhoro_hierarchical_2022}{minorbsem}.
#' Builds off work by \insertCite{wu_quantifying_2015;textual}{minorbsem}.
#' @param model A description of the user-specified model, lavaan syntax.
#' @param sample_cov (list of matrices) sample variance-covariance matrices.
#' The rownames and/or colnames must contain the observed variable names.
#' @param sample_nobs (vector of positive integer) Number of observations
#' for each study.
#' @param method (character) One of "normal", "lasso", "logistic",
#' "GDP", or "none". See details below.
#' @inheritParams minorbsem
#' @returns An object of \code{\link{mbsem-class}}
#' @details
#' There are different methods for estimating models in this package:
#'
#' - \code{normal}: under belief that minor factor influences are
#' on average zero with continuous deviations away from zero.
#' - \code{lasso}: under belief that minor factor influences are largely
#' zero with a small number of non-zero residual covariances.
#' - \code{logistic}: for similar belief as normal but more readily
#' accomodates extreme outliers.
#' - \code{GDP}: to mimic a global-local approach, i.e.
#' attempt to shrink near 0 residual covariances to 0
#' with minimal shrinking for larger residual covariances.
#' - \code{none}: if intending to ignore the influence of minor factors.
#' @references \insertAllCited{}
#' @keywords internal
meta_mbcfa <- function(
    model = NULL,
    sample_cov = NULL,
    sample_nobs = NULL,
    method = "normal",
    orthogonal = FALSE,
    simple_struc = TRUE,
    seed = 12345,
    warmup = 1000,
    sampling = 1000,
    refresh = (warmup + sampling) / 10,
    adapt_delta = .9,
    max_treedepth = 10,
    chains = 4,
    ncores = max(parallel::detectCores() - 2, 1),
    priors = new_mbsempriors(),
    show = TRUE,
    show_messages = TRUE) {
  # THIS IS A WIP!!
  message("Processing user input ...")

  # Model cannot be NULL
  user_input_check("model", model)

  # Priors must be class mbsempriors
  user_input_check("priors", priors)

  # method must be valid
  user_input_check("method-meta", method)

  # CmdStan path must be set
  tryCatch(cmdstanr::cmdstan_path(),
    error = function(e) {
      stop(paste0(
        "Error: CmdStan path has not been set yet.", " ",
        "See ?cmdstanr::set_cmdstan_path()."
      ))
    }
  )

  return(NULL)
}
