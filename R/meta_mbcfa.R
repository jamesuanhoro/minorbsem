#' Fit random-effects Bayesian meta-analytic CFAs with minor factors assumed.
#'
#' @description A function to fit fixed- or
#' random-effects Bayesian meta-analytic CFAs
#' with minor factors assumed \insertCite{uanhoro_hierarchical_2022}{minorbsem}.
#' Builds off work by \insertCite{wu_quantifying_2015;textual}{minorbsem}.
#' Does not yet accomodate moderators, and covariance matrices must be complete.
#' This will change in the near future.
#' @param model A description of the user-specified model, lavaan syntax.
#' @param data An optional data frame containing the observed variables used in
#' the model.
#' @param group An optional string identifying the grouping variable in
#' the data object.
#' @param sample_cov (list of matrices) sample variance-covariance matrices.
#' The rownames and/or colnames must contain the observed variable names.
#' For now, assumes there are no missing elements in the covariance matrices.
#' @param sample_nobs (vector of positive integer) Number of observations
#' for each study.
#' @param method (character) One of "normal", "lasso", "logistic",
#' "GDP", or "none". See details below.
#' @param type (character) One of "fe", "re", or "dep" for fixed-effects,
#' random-effects, and dependent-samples MASEM respectively.
#' The "dep" argument is experimental, see details below.
#' @param cluster An optional integer vector identifying the cluster each group
#' belongs to.
#' Asssume there are five groups, the first three belong to cluster 1
#' and the last two belong to cluster 2,
#' then the argument would be: \code{cluster = c(1, 1, 1, 2, 2)}.
#' This feature is experimental, see details below.
#' @inheritParams minorbsem
#' @returns An object of \code{\link{mbsem-class}}
#' @details
#' CFAs assume standardized factors.
#' Latent variable regression models are not yet implemented.
#'
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
#' with minimal shrinking for larger residual covariances
#' \insertCite{armagan_generalized_2013}{minorbsem}.
#' - \code{none}: if intending to ignore the influence of minor factors.
#'
#' When \code{type = "dep"}, the user must supply the cluster IDs, see cluster
#' parameter documentation above. However, this feature is experimental and only
#' available when \code{target = "cmdstan"}.
#' Additionally, the cluster inputs are not validated.
#' @examples
#' \dontrun{
#' model_syntax <- "# latent variable definitions
#' F1 =~ JP1 + JP2 + JP3
#' F2 =~ JN1 + JN2 + JN4 + JN4
#' F3 =~ TD1 + TD2"
#' meta_mbcfa(
#'   model_syntax,
#'   sample_cov = issp89$data, sample_nobs = issp89$n
#' )
#' model_syntax <- paste0(
#'   "distress =~ ", paste0("x", 1:14, collapse = " + "), "\n",
#'   "anxiety =~ ", paste0("x", seq(1, 14, 2), collapse = " + "), "\n",
#'   "depression =~ ", paste0("x", seq(2, 14, 2), collapse = " + ")
#' )
#' meta_mbcfa(
#'   model_syntax,
#'   sample_cov = Norton13$data, sample_nobs = Norton13$n, orthogonal = TRUE
#' )
#' }
#' @references \insertAllCited{}
#' @export
meta_mbcfa <- function(
    model = NULL,
    data = NULL,
    group = NULL,
    sample_cov = NULL,
    sample_nobs = NULL,
    method = "normal",
    type = "re",
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
    show_messages = TRUE,
    cluster = NULL,
    target = "rstan") {
  message("Processing user input ...")

  # Model cannot be NULL
  user_input_check("model", model)

  # Priors must be class mbsempriors
  user_input_check("priors", priors)

  # method must be valid
  user_input_check("method-meta", method)

  # type must be valid
  # Future: when moderators are added, warn user that fixed-effects
  # ignores moderators.
  user_input_check("type-meta", type)

  # Must provide either data and group or sample_cov and sample_nobs
  user_input_check("data-meta", data, group, sample_cov, sample_nobs)

  # target must be valid
  user_input_check("target", target)

  # check for cluster when type = "dep"
  user_input_check("meta-cluster", type, target, cluster)

  # Run lavaan fit
  if (!is.null(data)) {
    lav_fit <- lavaan::cfa(
      model,
      data = data, group = group, std.lv = TRUE,
      do.fit = FALSE, se = "none", test = "none", orthogonal = orthogonal
    )
  } else {
    lav_fit <- lavaan::cfa(
      model,
      sample.cov = sample_cov, sample.nobs = sample_nobs, std.lv = TRUE,
      do.fit = FALSE, se = "none", test = "none", orthogonal = orthogonal
    )
  }

  # Obtain data list for Stan
  data_list <- create_data_list_meta(
    lavaan_object = lav_fit,
    method = method,
    type = type,
    simple_struc = simple_struc,
    priors = priors,
    cluster = cluster
  )

  message("User input fully processed :)\n Now to modeling.")

  stan_fit <- target_fitter(
    target, data_list, seed, warmup, sampling, refresh,
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
