#' Fit random-effects Bayesian meta-analytic CFAs with minor factors assumed.
#'
#' @description A function to fit random-effects Bayesian meta-analytic CFAs
#' with minor factors assumed \insertCite{uanhoro_hierarchical_2022}{minorbsem}.
#' Builds off work by \insertCite{wu_quantifying_2015;textual}{minorbsem}.
#' @param data An optional data frame containing the observed variables used in
#' the model.
#' @param group An optional string identifying the grouping variable in
#' the data object.
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
#' @examples
#' \dontrun{
#' meta_mbcfa("# latent variable definitions
#'             F1 =~ JP1 + JP2 + JP3
#'             F2 =~ JN1 + JN2 + JN4 + JN4
#'             F3 =~ TD1 + TD2",
#'   sample_cov = issp89$data, sample_nobs = issp89$n
#' )
#' }
#' @references \insertAllCited{}
#' @keywords internal
meta_mbcfa <- function(
    model = NULL,
    data = NULL,
    group = NULL,
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
  message("Processing user input ...")

  # Model cannot be NULL
  user_input_check("model", model)

  # Priors must be class mbsempriors
  user_input_check("priors", priors)

  # method must be valid
  user_input_check("method-meta", method)

  # Must provide either data and group or sample_cov and sample_nobs
  user_input_check("data-meta", data, group, sample_cov, sample_nobs)

  # CmdStan path must be set
  tryCatch(cmdstanr::cmdstan_path(),
    error = function(e) {
      stop(paste0(
        "Error: CmdStan path has not been set yet.", " ",
        "See ?cmdstanr::set_cmdstan_path()."
      ))
    }
  )

  # Run lavaan fit
  if (!is.null(data)) {
    suppressWarnings(lav_fit <- lavaan::cfa(
      model,
      data = data, group = group, std.lv = TRUE,
      se = "none", test = "none", orthogonal = orthogonal
    ))
  } else {
    suppressWarnings(lav_fit <- lavaan::cfa(
      model,
      sample.cov = sample_cov, sample.nobs = sample_nobs, std.lv = TRUE,
      se = "none", test = "none", orthogonal = orthogonal
    ))
  }

  # Obtain data list for Stan
  data_list <- create_data_list_meta(
    lavaan_object = lav_fit,
    method = method,
    simple_struc = simple_struc,
    priors = priors
  )

  message("User input fully processed :)\n Now to modeling.")

  message(paste0(
    "Compiling Stan code ...\n",
    "This takes a while the first time you run a CFA ",
    "and the first time you run an SEM"
  ))

  if (data_list$sem_indicator == 0) {
    mod_resid <- cmdstanr::cmdstan_model(
      system.file("Stan/meta_cfa_resid.stan", package = "minorbsem"),
      stanc_options = list("O1")
    )
  } else if (data_list$sem_indicator == 1) {
    mod_resid <- cmdstanr::cmdstan_model(
      system.file("Stan/meta_sem_resid.stan", package = "minorbsem"),
      stanc_options = list("O1")
    )
  }

  message("Fitting Stan model ...")

  stan_fit <- mod_resid$sample(
    data = data_list,
    seed = seed,
    iter_warmup = warmup,
    iter_sampling = sampling,
    refresh = refresh,
    init = function() {
      list(
        resids = rep(1e-3, data_list$Ni^2 - data_list$Ni)
      )
    },
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    chains = chains,
    parallel_chains = ncores,
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
