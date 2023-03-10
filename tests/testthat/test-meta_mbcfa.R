test_that("Random method works for meta-CFA on issp89", {
  method <- random_method_selection(meta = TRUE)
  model_syntax <- "# latent variable definitions
    F1 =~ JP1 + JP2 + JP3
    F2 =~ JN1 + JN2 + JN4 + JN4
    F3 =~ TD1 + TD2"

  expect_error(fit <- meta_mbcfa(
    model_syntax,
    sample_cov = issp89$data[1:3], sample_nobs = issp89$n[1:3],
    orthogonal = sample(c(TRUE, FALSE), 1),
    simple_struc = sample(c(TRUE, FALSE), 1),
    warmup = 500, sampling = 500, chains = 3,
    method = method, refresh = 0, show_messages = FALSE
  ), NA)
  expect_true(all(slotNames(fit) %in% c(
    "major_parameters", "minor_factor_matrix", "data_list",
    "priors", "stan_fit", "version"
  )))
  mbsem_test_hist(fit)
  mbsem_test_trace(fit)
  mbsem_test_plot_residuals(fit, method)
  expect_error(
    kbl <- pretty_print_summary(fit),
    NA
  )
  mbsem_test_kbls_shared(kbl, method, meta = TRUE)
  expect_true(regexpr("Residual variances", kbl, ignore.case = TRUE) > 0)
  testthat::expect_error(
    casewise_log_likelihood(fit),
    "Not yet implemented for meta-analysis models"
  )
})

test_that("Random method works for meta-CFA on Norton13", {
  method <- random_method_selection(meta = TRUE)
  model_syntax <- paste0(
    "distress =~ ", paste0("x", 1:14, collapse = " + "), "\n",
    "anxiety =~ ", paste0("x", seq(1, 14, 2), collapse = " + "), "\n",
    "depression =~ ", paste0("x", seq(2, 14, 2), collapse = " + ")
  )

  expect_error(fit <- meta_mbcfa(
    model_syntax,
    sample_cov = Norton13$data[1:3], sample_nobs = Norton13$n[1:3],
    orthogonal = TRUE, simple_struc = TRUE,
    warmup = 500, sampling = 500, chains = 3,
    method = method, refresh = 0, show_messages = FALSE
  ), NA)
  expect_true(all(slotNames(fit) %in% c(
    "major_parameters", "minor_factor_matrix", "data_list",
    "priors", "stan_fit", "version"
  )))
  mbsem_test_hist(fit)
  mbsem_test_trace(fit)
  mbsem_test_plot_residuals(fit, method)
  expect_error(
    kbl <- pretty_print_summary(fit),
    NA
  )
  mbsem_test_kbls_shared(kbl, method, meta = TRUE)
  expect_true(regexpr("Residual variances", kbl, ignore.case = TRUE) > 0)
  testthat::expect_error(
    casewise_log_likelihood(fit),
    "Not yet implemented for meta-analysis models"
  )
})
