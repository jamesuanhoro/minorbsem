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
#' @param type (character) One of "fe" or "re" for fixed-effects MASEM
#' or random-effects MASEM respectively.
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
    show_messages = TRUE) {
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
    priors = priors
  )

  message("User input fully processed :)\n Now to modeling.")

  if (data_list$sem_indicator == 0) {
    # mod_resid <- stanmodels$meta_cfa_resid_rs
    mod_resid <- list()
  }

  message("Fitting Stan model ...")

  stan_fit <- rstan::sampling(
    mod_resid,
    data = data_list,
    chains = chains,
    cores = ncores,
    seed = seed,
    warmup = warmup,
    iter = warmup + sampling,
    refresh = refresh,
    init = function() {
      list(
        resids = rep(1e-3, (data_list$Ni^2 - data_list$Ni) / 2)
      )
    },
    control = list(
      adapt_delta = adapt_delta,
      max_treedepth = max_treedepth
    ),
    show_messages = show_messages
  )

  mbsem_results <- clean_up_stan_fit(
    stan_fit = stan_fit, data_list = data_list, priors = priors
  )
  if (show) {
    show(mbsem_results)
  }

  return(mbsem_results)
}
