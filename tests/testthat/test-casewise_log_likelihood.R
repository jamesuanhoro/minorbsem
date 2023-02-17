test_that("Random method (any case) works for CFA", {
  method <- random_method_selection()
  model_syntaxes <- c(
    "F1 =~ x1 + x2 + x3\n F2 =~ x4 + x5 + x6\n F3 =~ x7 + x8 + x9",
    "F1 =~ x1 + x2 + x3 + x9\n F2 =~ x4 + x5 + x6\n F3 =~ x7 + x8 + x9",
    "F1 =~ x1 + x2 + x3 + x9\n F2 =~ x4 + x5 + x6\n F3 =~ x7 + x8 + x9
      x2 ~~ x7\n  x3 ~~ x5"
  )
  model_syntax <- sample(model_syntaxes, 1)
  fit <- minorbsem(
    model_syntax, HS,
    warmup = 500, sampling = 500, chains = 3,
    method = method, refresh = 0, show_messages = FALSE
  )
  residuals <- sample(c(TRUE, FALSE), 1)
  if (tolower(method) == "none" && isTRUE(residuals)) {
    warn_msg <- paste0(
      "include_residuals = TRUE is ignored when ",
      "minorbsem method == \"none\". "
    )
    expect_warning(
      ll_mat <- casewise_log_likelihood(
        fit,
        include_residuals = residuals
      ),
      warn_msg
    )
  } else {
    expect_error(
      ll_mat <- casewise_log_likelihood(
        fit,
        include_residuals = residuals
      ),
      NA
    )
  }
  expect_equal(ncol(ll_mat), nrow(HS))
  expect_equal(nrow(ll_mat), 500 * 3)
  expect_true(!is.na(sum(ll_mat)))
})

test_that("Random method (any case) works for SEM", {
  method <- random_method_selection()
  model_syntax <- "
  ind60 =~ x1 + x2 + x3\n dem60 =~ y1 + y2 + y3 + y4\n
  dem65 =~ y5 + y6 + y7 + y8\n dem60 ~ ind60\n dem65 ~ ind60 + dem60"
  fit <- minorbsem(
    model_syntax, PD,
    warmup = 500, sampling = 500, chains = 3,
    method = method, refresh = 0, show_messages = FALSE
  )
  residuals <- sample(c(TRUE, FALSE), 1)
  if (tolower(method) == "none" && isTRUE(residuals)) {
    warn_msg <- paste0(
      "include_residuals = TRUE is ignored when ",
      "minorbsem method == \"none\". "
    )
    expect_warning(
      ll_mat <- casewise_log_likelihood(
        fit,
        include_residuals = residuals
      ),
      warn_msg
    )
  } else {
    expect_error(
      ll_mat <- casewise_log_likelihood(
        fit,
        include_residuals = residuals
      ),
      NA
    )
  }
  expect_equal(ncol(ll_mat), nrow(PD))
  expect_equal(nrow(ll_mat), 500 * 3)
  expect_true(!is.na(sum(ll_mat)))
})
