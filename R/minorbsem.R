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
#' @param orthogonal (logical) constrain factors orthogonal, must be TRUE to fit
#' bifactor models.
#' @param seed (positive integer) seed, set to obtain replicable results.
#' @param warmup (positive integer) The number of warmup iterations to run per
#' chain.
#' @param sampling (positive integer) The number of post-warmup iterations to
#' run per chain, retained for inference.
#' @param adapt_delta (real in (0, 1)) Increase to resolve divergent
#' transitions.
#' @param max_treedepth (positive integer) Increase to resolve problems with
#' maximum tree depth.
#' @param chains (positive integer) The number of Markov chains to run.
#' @param ncores (positive integer) The number of chains to run in parallel.
#' @param lkj_shape (positive real) The shape parameter of the LKJ-prior on the
#' interfactor correlation matrix.
#' @param sl_par (positive real) The scale parameter of the
#' Student-t(df = 3, loc = 0) prior on the hyper-parameter of the standard
#' deviation of loadings.
#' @param rs_par (positive real) The scale parameter of the
#' Student-t(df = 3, loc = 0) prior on the residual standard deviations.
#' @param rc_par (positive real) The shape parameter of the Beta(rc_par, rc_par)
#' prior on the residual error correlations.
#' @param sc_par (positive real) The scale parameter of the
#' Student-t(df = 3, loc = 0) prior on the hyper-parameter of the standard
#' deviations of coefficients; SD(coefs) vary by outcome.
#' @returns A list
#' @examples
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
#' @export
minorbsem <- function(
    model = NULL,
    data = NULL,
    sample_cov = NULL,
    sample_nobs = NULL,
    orthogonal = FALSE,
    seed = 12345,
    warmup = 500,
    sampling = 500,
    adapt_delta = .9,
    max_treedepth = 10,
    chains = 3,
    ncores = max(parallel::detectCores() - 2, 1),
    lkj_shape = 2.0,
    sl_par = 1.0,
    rs_par = 2.5,
    rc_par = 2.0,
    sc_par = 1.0) {
  message("Processing user input ...")

  # Model cannot be NULL
  if (is.null(model)) {
    stop("Model cannot be null")
  }

  # Must provide either data or sample_cov and sample_nobs
  if (is.null(data) & (is.null(sample_cov) | is.null(sample_nobs))) {
    stop("User must provide either:\n\t(i) a dataset or\n\t(ii) sample covariance and sample size")
  }

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
  data_list <- create_data_list(lav_fit, lkj_shape, sl_par, rs_par, rc_par, sc_par)

  message("User input fully processed :)\n Now to modeling.")

  message("Compiling Stan code ...")

  # TODO: This should be a package-level global config setting up by the user
  cmdstanr::set_cmdstan_path("~/cmdstan/")

  if (data_list$sem_indicator == 0) {
    mod_resid <- cmdstanr::cmdstan_model(
      "src/inst/cfa_resid_nrm.stan",
      stanc_options = list("O1")
    )
  } else if (data_list$sem_indicator == 1) {
    mod_resid <- cmdstanr::cmdstan_model(
      "src/inst/sem_resid_nrm.stan",
      stanc_options = list("O1")
    )
  }

  message("Fitting Stan model ...")

  stan_fit <- mod_resid$sample(
    data = data_list,
    seed = seed,
    iter_warmup = warmup,
    iter_sampling = sampling,
    refresh = (warmup + sampling) / 10,
    init = function() {
      list(
        resids = rep(1e-3, data_list$Ni^2 - data_list$Ni)
      )
    },
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    chains = chains,
    parallel_chains = ncores
  )

  clean_results <- clean_up_stan_fit(stan_fit, data_list)
  pretty_print_summary(clean_results)

  return(clean_results)
}
