update_cmdstan_loc("~/cmdstan/")

test_that("Random method (any case) works for CFA", {
  method <- random_method_selection()
  model_syntaxes <- c(
    "F1 =~ x1 + x2 + x3\n F2 =~ x4 + x5 + x6\n F3 =~ x7 + x8 + x9",
    "F1 =~ x1 + x2 + x3 + x9\n F2 =~ x4 + x5 + x6\n F3 =~ x7 + x8 + x9",
    "F1 =~ x1 + x2 + x3 + x9\n F2 =~ x4 + x5 + x6\n F3 =~ x7 + x8 + x9
      x2 ~~ x7\n  x3 ~~ x5"
  )
  syntax_idx <- sample(length(model_syntaxes), 1)
  model_syntax <- model_syntaxes[syntax_idx]
  fit <- minorbsem(
    model_syntax, HS,
    method = method, refresh = 0, show_messages = FALSE
  )
  expect_error(
    kbl <- pretty_print_summary(fit),
    NA
  )
  expect_true(inherits(kbl, "kableExtra"))
  expect_true(regexpr(method, kbl, ignore.case = TRUE) > 0)
  expect_true(regexpr("Factor loadings", kbl, ignore.case = TRUE) > 0)
  expect_true(regexpr("Residual variances", kbl, ignore.case = TRUE) > 0)
  expect_true(regexpr(
    "Inter-factor correlations",
    kbl,
    ignore.case = TRUE
  ) > 0)
  if (syntax_idx == 3) {
    expect_true(regexpr(
      "Error correlations",
      kbl,
      ignore.case = TRUE
    ) > 0)
  }
})

test_that("Random method (any case) works for SEM", {
  method <- random_method_selection()
  model_syntax <- "
  ind60 =~ x1 + x2 + x3\n dem60 =~ y1 + y2 + y3 + y4\n
  dem65 =~ y5 + y6 + y7 + y8\n dem60 ~ ind60\n dem65 ~ ind60 + dem60"
  fit <- minorbsem(
    model_syntax, PD,
    method = method, refresh = 0, show_messages = FALSE
  )
  expect_error(
    kbl <- pretty_print_summary(fit),
    NA
  )
  expect_true(inherits(kbl, "kableExtra"))
  expect_true(regexpr(method, kbl, ignore.case = TRUE) > 0)
  expect_true(regexpr(
    "Latent regression coefficients",
    kbl,
    ignore.case = TRUE
  ) > 0)
  expect_true(regexpr("Factor loadings", kbl, ignore.case = TRUE) > 0)
  expect_true(regexpr("R square", kbl, ignore.case = TRUE) > 0)
})

update_cmdstan_loc("")