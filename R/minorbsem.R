#' Fit Bayesian SEMs with minor factors assumed
#'
#' @description The major function to fit models
#' @param model A description of the user-specified model, lavaan syntax.
#' @param data An optional data frame containing the observed variables used in
#' the model.
#' @param sample_cov (matrix) sample variance-covariance matrix.
#' The rownames and/or colnames must contain the observed variable names.
#' @param sample_nobs (positive integer) Number of observations if the full
#' data frame
#' is missing and only sample covariance matrix is given.
#' @param method (character) One of "normal", "lasso", "logistic",
#' "GDP" or "none".
#' Select "normal" under belief that minor factor influences are
#' on average zero with continuous deviations away from zero.
#' Select "lasso" under belief that minor factor influences are
#' indeed zero with a small number of non-zero residual covariances.
#' Select "logistic" for similar belief as normal but more readily
#' accomodates extreme outliers.
#' Select "GDP" to mimic a global-local approach, i.e.
#' attempt to shrink near 0 residual covariances to 0
#' with minimal shrinking for larger residual covariances.
#' Select "none" if intending to ignore the influence of minor factors.
#' @param orthogonal (logical) constrain factors orthogonal, must be TRUE to fit
#' bifactor models.
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
#' @returns An object of \code{\link{mbsem-class}}
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
#' @export
minorbsem <- function(
    model = NULL,
    data = NULL,
    sample_cov = NULL,
    sample_nobs = NULL,
    method = "normal",
    orthogonal = FALSE,
    seed = 12345,
    warmup = 500,
    sampling = 500,
    refresh = (warmup + sampling) / 10,
    adapt_delta = .9,
    max_treedepth = 10,
    chains = 3,
    ncores = max(parallel::detectCores() - 2, 1),
    priors = new_mbsempriors(),
    show = TRUE,
    show_messages = TRUE) {
  message("Processing user input ...")

  # Model cannot be NULL
  user_input_check("model", model)

  # Priors must be class mbsempriors
  user_input_check("model", priors)

  # method must be valid
  user_input_check("method", method)

  # Must provide either data or sample_cov and sample_nobs
  user_input_check("data", data, sample_cov, sample_nobs)

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
    lav_fit <- lavaan::cfa(
      model,
      data = data,
      std.lv = TRUE,
      orthogonal = orthogonal
    )
  } else {
    lav_fit <- lavaan::cfa(
      model,
      sample.cov = sample_cov, sample.nobs = sample_nobs,
      std.lv = TRUE,
      orthogonal = orthogonal
    )
  }

  # Obtain data list for Stan
  data_list <- create_data_list(lav_fit, method, priors)

  message("User input fully processed :)\n Now to modeling.")

  message(paste0(
    "Compiling Stan code ...\n",
    "This takes a while the first time you run a CFA ",
    "and the first time you run an SEM"
  ))

  if (data_list$sem_indicator == 0) {
    mod_resid <- cmdstanr::cmdstan_model(
      system.file("Stan/cfa_resid.stan", package = "minorbsem"),
      stanc_options = list("O1")
    )
  } else if (data_list$sem_indicator == 1) {
    mod_resid <- cmdstanr::cmdstan_model(
      system.file("Stan/sem_resid.stan", package = "minorbsem"),
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

  mbsem_results <- clean_up_stan_fit(stan_fit, data_list, priors)
  if (show) {
    show(mbsem_results)
  }

  return(mbsem_results)
}
