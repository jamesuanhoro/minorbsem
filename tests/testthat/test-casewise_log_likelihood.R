mbsem_test_ll_1 <- function(fit, method) {
  residuals <- sample(c(TRUE, FALSE), 1)
  if (tolower(method) == "none" && isTRUE(residuals)) {
    warn_msg <- paste0(
      "include_residuals = TRUE is ignored when ",
      "minorbsem method == \"none\". "
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
  testthat::expect_equal(ncol(ll_mat), nrow(fit@data_list$Y))
  testthat::expect_equal(nrow(ll_mat), 500 * 3)
  testthat::expect_true(!is.na(sum(ll_mat)))
  testthat::expect_true(sum(abs(ll_mat)) > 0)
}

mbsem_test_ll_2 <- function(fit, method) {
  residuals <- sample(c(TRUE, FALSE), 1)
  if (tolower(method) == "none" && isTRUE(residuals)) {
    warn_msg <- paste0(
      "include_residuals = TRUE is ignored when ",
      "minorbsem method == \"none\". "
    )
    testthat::expect_warning(
      ll_mat_arma <- casewise_log_likelihood(
        fit,
        include_residuals = residuals,
        use_armadillo = TRUE
      ),
      warn_msg
    )
    testthat::expect_warning(
      ll_mat_base <- casewise_log_likelihood(
        fit,
        include_residuals = residuals,
        use_armadillo = FALSE
      ),
      warn_msg
    )
    testthat::expect_true(sum(abs(ll_mat_arma - ll_mat_base)) < 1e-5)
  } else {
    ll_mat_arma <- casewise_log_likelihood(
      fit,
      include_residuals = residuals,
      use_armadillo = TRUE
    )
    ll_mat_base <- casewise_log_likelihood(
      fit,
      include_residuals = residuals,
      use_armadillo = FALSE
    )
    testthat::expect_true(sum(abs(ll_mat_arma - ll_mat_base)) < 1e-5)
  }
}

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
    orthogonal = sample(c(TRUE, FALSE), 1),
    simple_struc = sample(c(TRUE, FALSE), 1),
    warmup = 500, sampling = 500, chains = 3,
    method = method, refresh = 0, show_messages = FALSE
  )
  mbsem_test_ll_1(fit, method)
})

test_that("Random method (any case) works for SEM", {
  method <- random_method_selection()
  model_syntax <- "
  ind60 =~ x1 + x2 + x3\n dem60 =~ y1 + y2 + y3 + y4\n
  dem65 =~ y5 + y6 + y7 + y8\n dem60 ~ ind60\n dem65 ~ ind60 + dem60"
  fit <- minorbsem(
    model_syntax, PD,
    orthogonal = sample(c(TRUE, FALSE), 1),
    simple_struc = sample(c(TRUE, FALSE), 1),
    warmup = 500, sampling = 500, chains = 3,
    method = method, refresh = 0, show_messages = FALSE
  )
  mbsem_test_ll_1(fit, method)
})

test_that("CFA: Different LL methods are equal", {
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
    orthogonal = sample(c(TRUE, FALSE), 1),
    simple_struc = sample(c(TRUE, FALSE), 1),
    warmup = 500, sampling = 500, chains = 3,
    method = method, refresh = 0, show_messages = FALSE
  )
  mbsem_test_ll_2(fit, method)
})

test_that("SEM: Different LL methods are equal", {
  method <- random_method_selection()
  model_syntax <- "
  ind60 =~ x1 + x2 + x3\n dem60 =~ y1 + y2 + y3 + y4\n
  dem65 =~ y5 + y6 + y7 + y8\n dem60 ~ ind60\n dem65 ~ ind60 + dem60"
  fit <- minorbsem(
    model_syntax, PD,
    orthogonal = sample(c(TRUE, FALSE), 1),
    simple_struc = sample(c(TRUE, FALSE), 1),
    warmup = 500, sampling = 500, chains = 3,
    method = method, refresh = 0, show_messages = FALSE
  )
  mbsem_test_ll_2(fit, method)
})
