test_that("Random method prints correctly for CFA", {
  expect_error(
    print_out <- paste0(
      capture_output(pretty_print_summary(fit_cfa), width = 300),
      collapse = "\n"
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
})

test_that("Random method prints correctly for SEM", {
  expect_error(
    print_out <- paste0(
      capture_output(pretty_print_summary(fit_sem), width = 300),
      collapse = "\n"
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
})
