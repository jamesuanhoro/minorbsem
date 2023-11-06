#' Fit Bayesian SEMs with minor factors assumed
#'
#' @description The major function to fit models assuming the influence of
#' minor factors \insertCite{uanhoro_modeling_2023}{minorbsem}.
#' @param model A description of the user-specified model, lavaan syntax.
#' @param data An optional data frame containing the observed variables used in
#' the model.
#' @param sample_cov (matrix) sample variance-covariance matrix.
#' The rownames and/or colnames must contain the observed variable names.
#' @param sample_nobs (positive integer) Number of observations if the full
#' data frame
#' is missing and only sample covariance matrix is given.
#' @param method (character) One of "normal", "lasso", "logistic",
#' "GDP", "WB", "WB-cond", "WW", or "none". See details below.
#' @param orthogonal (logical) constrain factors orthogonal, must be TRUE to fit
#' bifactor models.
#' @param simple_struc (LOGICAL) Only relevant for CFAs.
#' If TRUE: assume simple structure;
#' If FALSE: estimate all cross-loadings using generalized
#' double Pareto priors.
#' @param seed (positive integer) seed, set to obtain replicable results.
#' @param warmup (positive integer) The number of warmup iterations to run per
#' chain.
#' @param sampling (positive integer) The number of post-warmup iterations to
#' run per chain, retained for inference.
#' @param refresh (positive integer) How often to print the status
#' of the sampler.
#' @param adapt_delta (real in (0, 1)) Increase to resolve divergent
#' transitions.
#' @param max_treedepth (positive integer) Increase to resolve problems with
#' maximum tree depth.
#' @param chains (positive integer) The number of Markov chains to run.
#' @param ncores (positive integer) The number of chains to run in parallel.
#' @param priors An object of \code{\link{mbsempriors-class}}.
#' See \code{\link{new_mbsempriors}} for more information.
#' @param show (Logical) If TRUE, show table of results, if FALSE, do not
#' show table of results. As an example, use FALSE for simulation studies.
#' @param show_messages (Logical) If TRUE, show messages from Stan sampler,
#' if FALSE, hide messages.
#' @param compute_ll (Logical) If TRUE, compute log-likelihood,
#' if FALSE, do not. This may be useful for cross-validation. This argument is
#' ignored when:
#' (i) the full dataset is not provided;
#' (ii) the method is WB, use WB-cond instead.
#' @returns An object of \code{\link{mbsem-class}}
#' @details
#' CFAs assume standardized factors.
#' Latent variable regression models print results with standardized loadings.
#'
#' There are different methods for estimating models in this package:
#'
#' - \code{normal}: under belief that minor factor influences are
#' on average zero with continuous deviations away from zero
#' \insertCite{uanhoro_modeling_2023}{minorbsem}.
#' - \code{lasso}: under belief that minor factor influences are largely
#' zero with a small number of non-zero residual covariances.
#' - \code{logistic}: for similar belief as normal but more readily
#' accomodates extreme outliers.
#' - \code{GDP}: to mimic a global-local approach, i.e.
#' attempt to shrink near 0 residual covariances to 0
#' with minimal shrinking for larger residual covariances
#' \insertCite{armagan_generalized_2013}{minorbsem}.
#' - \code{WB}: to model the covariance matrix hierarchically under
#' assumptions of adventitiuous error
#' \insertCite{wu_quantifying_2015}{minorbsem};
#' does NOT allow for computation of casewise log-likelihoods and LOO-CV.
#' - \code{WB-cond}: same as WB but estimates the
#' "population covariance matrix",
#' allowing for computation of casewise log-likelihoods and LOO-CV.
#' - \code{WW}: A variation on WB-cond, but assumes the population
#' covariance matrix is Wishart as opposed to inverse-Wishart;
#' - \code{none}: if intending to ignore the influence of minor factors.
#' @examples
#' \dontrun{
#' minorbsem("# latent variable definitions
#'            F1 =~ x1 + x2 + x3
#'            F2 =~ x4 + x5 + x6
#'            F3 =~ x7 + x8 + x9", HS)
#' minorbsem("# latent variable definitions
#'            ind60 =~ x1 + x2 + x3
#'            dem60 =~ y1 + y2 + y3 + y4
#'            dem65 =~ y5 + y6 + y7 + y8
#'            # latent regressions
#'            dem60 ~ ind60
#'            dem65 ~ ind60 + dem60", PD)
#' }
#' @references \insertAllCited{}
#' @export
minorbsem <- function(
    model = NULL,
    data = NULL,
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
    chains = 3,
    ncores = max(parallel::detectCores() - 2, 1),
    priors = new_mbsempriors(),
    show = TRUE,
    show_messages = TRUE,
    compute_ll = FALSE) {
  message("Processing user input ...")

  # Model cannot be NULL
  user_input_check("model", model)

  # Priors must be class mbsempriors
  user_input_check("priors", priors)

  # method must be valid
  user_input_check("method", method)

  # Must provide either data or sample_cov and sample_nobs
  user_input_check("data", data, sample_cov, sample_nobs)

  # Run lavaan fit
  if (!is.null(data)) {
    lav_fit <- lavaan::cfa(
      model,
      data = data,
      std.lv = TRUE,
      do.fit = FALSE,
      likelihood = "wishart",
      do.fit = FALSE,
      ceq.simple = FALSE,
      orthogonal = orthogonal
    )
  } else {
    lav_fit <- lavaan::cfa(
      model,
      sample.cov = sample_cov, sample.nobs = sample_nobs,
      std.lv = TRUE,
      do.fit = FALSE,
      likelihood = "wishart",
      do.fit = FALSE,
      ceq.simple = FALSE,
      orthogonal = orthogonal
    )
  }

  # Obtain data list for Stan
  data_list <- create_data_list(
    lavaan_object = lav_fit,
    method = method,
    simple_struc = simple_struc,
    priors = priors,
    compute_ll = compute_ll,
    model = model,
    orthogonal = orthogonal
  )

  message("User input fully processed :)\n Now to modeling.")

  stan_fit <- target_fitter(
    data_list, seed, warmup, sampling, refresh,
    adapt_delta, max_treedepth, chains, ncores, show_messages
  )

  mbsem_results <- clean_up_stan_fit(
    stan_fit = stan_fit, data_list = data_list, priors = priors
  )
  if (show) {
    show(mbsem_results)
  }

  return(mbsem_results)
}
