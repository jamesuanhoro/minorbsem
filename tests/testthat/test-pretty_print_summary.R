test_that("Random method prints correctly for CFA", {
  skip_if_not_installed("cmdstanr")
  expect_error(
    print_out <- capture_output(
      pretty_print_summary(fit_cfa), width = 300
    ),
    NA
  )
  mbsem_test_pp_shared(print_out, method_cfa)
  expect_true(grepl("Residual variances", print_out, ignore.case = TRUE))
  if (!orthogonal_cfa) {
    expect_true(grepl(
      "Inter-factor correlations",
      print_out,
      ignore.case = TRUE
    ))
  }
  if (syntax_cfa_idx == 3) {
    expect_true(grepl(
      "Error correlations",
      print_out,
      ignore.case = TRUE
    ))
  }
  expect_error(
    print_out <- capture_output(
      pretty_print_summary(fit_cfa, simple = FALSE), width = 300
    ),
    NA
  )
  mbsem_test_pp_shared(print_out, method_cfa, simple = FALSE)
})

test_that("Random method prints correctly for SEM", {
  skip_if_not_installed("cmdstanr")
  expect_error(
    print_out <- capture_output(
      pretty_print_summary(fit_sem), width = 300
    ),
    NA
  )
  mbsem_test_pp_shared(print_out, method_sem)
  expect_true(grepl(
    "Latent regression coefficients",
    print_out,
    ignore.case = TRUE
  ))
  expect_true(grepl("R square", print_out, ignore.case = TRUE))
  expect_error(
    print_out <- capture_output(
      pretty_print_summary(fit_sem, simple = FALSE), width = 300
    ),
    NA
  )
  mbsem_test_pp_shared(print_out, method_sem, simple = FALSE)
})
