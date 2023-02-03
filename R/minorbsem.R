minorbsem <- function (
    model = NULL,
    data = NULL,
    sample_cov = NULL,
    sample_nobs = NULL,
    orthogonal = FALSE
) {

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
      model, data = data,
      std.lv = TRUE,
      orthogonal = orthogonal)
  } else {
    lav_fit <- lavaan::cfa(
      model, sample.cov = sample_cov, sample.nobs = sample_nobs,
      std.lv = TRUE,
      orthogonal = orthogonal)
  }

  # Obtain data list for Stan
  data_list <- create_data_list(lav_fit)

  message("User input fully processed :)\n Now to modeling.")

  message("Compiling Stan code ...")

  return(data_list)
}
