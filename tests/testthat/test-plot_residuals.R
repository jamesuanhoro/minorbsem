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
    method = method, refresh = 0, show_messages = FALSE
  )
  type <- sample(c("matrix", "range"), 1)
  if (tolower(method) == "none") {
    expect_error(
      gg <- plot_residuals(fit, type = type),
      "There are no residuals to plot when method = \"none\""
    )
  } else {
    expect_error(
      gg <- plot_residuals(fit, type = type),
      NA
    )
    expect_true(inherits(gg, "ggplot"))
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
  type <- sample(c("matrix", "range"), 1)
  if (tolower(method) == "none") {
    expect_error(
      gg <- plot_residuals(fit, type = type),
      "There are no residuals to plot when method = \"none\""
    )
  } else {
    expect_error(
      gg <- plot_residuals(fit, type = type),
      NA
    )
    expect_true(inherits(gg, "ggplot"))
  }
})
