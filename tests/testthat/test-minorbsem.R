test_that("Random method (any case) works for CFA", {
  method <- random_method_selection()
  model_syntaxes <- c(
    "F1 =~ x1 + x2 + x3\n F2 =~ x4 + x5 + x6\n F3 =~ x7 + x8 + x9",
    "F1 =~ x1 + x2 + x3 + x9\n F2 =~ x4 + x5 + x6\n F3 =~ x7 + x8 + x9",
    "F1 =~ x1 + x2 + x3 + x9\n F2 =~ x4 + x5 + x6\n F3 =~ x7 + x8 + x9
      x2 ~~ x7\n  x3 ~~ x5"
  )
  model_syntax <- sample(model_syntaxes, 1)
  expect_error(fit <- minorbsem(
    model_syntax, HS,
    method = method, refresh = 0, show_messages = FALSE
  ), NA)
  expect_true(all(slotNames(fit) %in% c(
    "major_parameters", "minor_factor_matrix", "data_list",
    "priors", "stan_fit", "version"
  )))
})

test_that("Random method (any case) works for SEM", {
  method <- random_method_selection()
  model_syntax <- "
  ind60 =~ x1 + x2 + x3\n dem60 =~ y1 + y2 + y3 + y4\n
  dem65 =~ y5 + y6 + y7 + y8\n dem60 ~ ind60\n dem65 ~ ind60 + dem60"
  expect_error(fit <- minorbsem(
    model_syntax, PD,
    method = method, refresh = 0, show_messages = FALSE
  ), NA)
  expect_true(all(slotNames(fit) %in% c(
    "major_parameters", "minor_factor_matrix", "data_list",
    "priors", "stan_fit", "version"
  )))
})

test_that("Noncreated method fails", {
  method <- paste0(method_hash(), collapse = "")
  method <- strsplit(paste0(method_hash(), collapse = ""), "")[[1]]
  len_method <- length(method)
  method <- tolower(paste0(
    sample(method)[1:sample(len_method, 1)],
    collapse = ""
  ))
  model_syntaxes <- c(
    "F1 =~ x1 + x2 + x3\n F2 =~ x4 + x5 + x6\n F3 =~ x7 + x8 + x9",
    "F1 =~ x1 + x2 + x3 + x9\n F2 =~ x4 + x5 + x6\n F3 =~ x7 + x8 + x9",
    "F1 =~ x1 + x2 + x3 + x9\n F2 =~ x4 + x5 + x6\n F3 =~ x7 + x8 + x9
      x2 ~~ x7\n  x3 ~~ x5"
  )
  model_syntax <- sample(model_syntaxes)[1]
  err_msg <- paste0(
    "method must be one of the following: ",
    paste0("\"", method_hash(), "\"", collapse = ", ")
  )
  expect_error(fit <- minorbsem(
    model_syntax, HS,
    method = method, refresh = 0, show_messages = FALSE
  ), err_msg)
})
