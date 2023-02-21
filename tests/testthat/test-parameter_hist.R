mbsem_test_hist <- function(fit) {
  testthat::expect_error(
    gg <- parameter_hist(
      fit,
      param_type = c(
        "rm", "lo", "ev", "co", "rc", "fc", "rsq",
        "re"
      )
    ),
    NA
  )
  testthat::expect_true(inherits(gg, "ggplot"))
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
    method = method, refresh = 0, show_messages = FALSE
  )
  mbsem_test_hist(fit)
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
    method = method, refresh = 0, show_messages = FALSE
  )
  mbsem_test_hist(fit)
})
