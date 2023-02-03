minorbsem <- function(
    model = NULL,
    data = NULL,
    sample_cov = NULL,
    sample_nobs = NULL,
    orthogonal = FALSE,
    seed = 12345,
    warmup = 1000,
    sampling = 1000,
    adapt_delta = .9,
    max_treedepth = 10,
    chains = 3,
    ncores = parallel::detectCores() - 2,
    lkj_shape = 2,
    sl_par = 1,
    rs_par = 2.5,
    rc_par = 2.0) {
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
  data_list <- create_data_list(lav_fit, lkj_shape, sl_par, rs_par, rc_par)

  message("User input fully processed :)\n Now to modeling.")

  message("Compiling Stan code ...")

  # TODO: This should be a global object set up by the user
  cmdstanr::set_cmdstan_path("~/cmdstan/")
  mod_resid <- cmdstanr::cmdstan_model(
    "src/cfa_resid_nrm.stan",
    stanc_options = list("O1")
  )

  message("Fitting Stan code ...")

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
  return(stan_fit)
}
