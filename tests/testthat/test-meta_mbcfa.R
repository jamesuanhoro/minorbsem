mbsem_test_meta_ll <- function(fit, method) {
  residuals <- sample(c(TRUE, FALSE), 1)

  if (method_hash(method) >= 90 && isTRUE(residuals)) {
    warn_msg <- paste0(
      "include_residuals = TRUE is ignored when ",
      "minorbsem method == \"none\", \"WB\", \"WB-cond\"."
    )
    testthat::expect_warning(
      ll_mat <- casewise_log_likelihood(
        fit,
        include_residuals = residuals,
        use_armadillo = !residuals
      ),
      warn_msg
    )
  } else {
    testthat::expect_error(
      ll_mat <- casewise_log_likelihood(
        fit,
        include_residuals = residuals,
        use_armadillo = !residuals
      ),
      NA
    )
  }
  testthat::expect_equal(ncol(ll_mat), fit@data_list$Ng)
  testthat::expect_equal(nrow(ll_mat), 500 * 1)
  testthat::expect_true(!is.na(sum(ll_mat)))
  testthat::expect_true(sum(abs(ll_mat)) > 0)
}

test_that("Random method works for meta-CFA on issp89", {
  print(type <- sample(c("fe", "re"), 1))
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
    warmup = 500, sampling = 500, chains = 1,
    method = method, type = type,
    refresh = 0, show_messages = FALSE
  ), NA)
  expect_true(all(slotNames(fit) %in% c(
    "major_parameters", "minor_factor_matrix", "data_list",
    "priors", "stan_fit", "version"
  )))
  mbsem_test_hist(fit)
  mbsem_test_trace(fit)
  mbsem_test_plot_residuals(fit, method)
  expect_error(
    print_out <- paste0(
      capture.output(pretty_print_summary(fit)),
      collapse = "\n"
    ),
    NA
  )
  mbsem_test_pp_shared(print_out, method, meta = TRUE)
  expect_true(grepl("Residual variances", print_out, ignore.case = TRUE))
  mbsem_test_meta_ll(fit, method)
})

test_that("Random method works for meta-CFA on Norton13", {
  print(type <- sample(c("fe", "re"), 1))
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
    warmup = 500, sampling = 500, chains = 1,
    method = method, type = type,
    refresh = 0, show_messages = FALSE
  ), NA)
  expect_true(all(slotNames(fit) %in% c(
    "major_parameters", "minor_factor_matrix", "data_list",
    "priors", "stan_fit", "version"
  )))
  mbsem_test_hist(fit)
  mbsem_test_trace(fit)
  mbsem_test_plot_residuals(fit, method)
  expect_error(
    print_out <- paste0(
      capture.output(pretty_print_summary(fit)),
      collapse = "\n"
    ),
    NA
  )
  mbsem_test_pp_shared(print_out, method, meta = TRUE)
  expect_true(grepl("Residual variances", print_out, ignore.case = TRUE))
  mbsem_test_meta_ll(fit, method)
})
