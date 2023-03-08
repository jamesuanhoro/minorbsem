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

mbsem_test_trace <- function(fit) {
  testthat::expect_error(
    gg <- parameter_trace(
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

mbsem_test_plot_residuals <- function(fit, method) {
  type <- sample(c("matrix", "range"), 1)
  if (method_hash(method) >= 90) {
    testthat::expect_error(
      gg <- plot_residuals(fit, type = type),
      paste0(
        "There are no residuals to plot when ",
        "method == \"none\", \"WB\", \"WB-cond\"."
      )
    )
  } else {
    testthat::expect_error(
      gg <- plot_residuals(fit, type = type),
      NA
    )
    testthat::expect_true(inherits(gg, "ggplot"))
  }
}

mbsem_test_kbls_shared <- function(kbl, method, meta = FALSE) {
  testthat::expect_true(inherits(kbl, "kableExtra"))
  testthat::expect_true(regexpr(method, kbl, ignore.case = TRUE) > 0)
  testthat::expect_true(regexpr("Goodness of fit", kbl, ignore.case = TRUE) > 0)
  testthat::expect_true(regexpr("RMSE", kbl, ignore.case = TRUE) > 0)
  testthat::expect_true(regexpr("Factor loadings", kbl, ignore.case = TRUE) > 0)
  if (isTRUE(meta)) {
    testthat::expect_true(regexpr("RMSEA", kbl, ignore.case = TRUE) > 0)
  } else {
    testthat::expect_true(regexpr("PPP", kbl, ignore.case = TRUE) > 0)
  }
}
