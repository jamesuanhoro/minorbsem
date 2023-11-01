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
        "method == \"none\", \"WB\", \"WB-cond\", \"WW\"."
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

mbsem_test_pp_shared <- function(print_out, method) {
  testthat::expect_true(grepl(method, print_out, ignore.case = TRUE))
  testthat::expect_true(grepl("Goodness of fit", print_out, ignore.case = TRUE))
  testthat::expect_true(grepl("RMSE", print_out, ignore.case = TRUE))
  testthat::expect_true(grepl("Factor loadings", print_out, ignore.case = TRUE))
  testthat::expect_true(grepl("PPP", print_out, ignore.case = TRUE))
}

dat_cov <- function(dat = "HS") {
  if (dat == "HS") {
    dat <- HS[, paste0("x", 1:9)] # nolint
  } else if (dat == "PD") {
    dat <- PD # nolint
  }

  input <- list(data = NULL, cov = NULL, nobs = NULL)
  if (sample(2, 1) == 1) {
    input$data <- dat
  } else {
    input$cov <- cov(dat)
    input$nobs <- nrow(dat)
  }

  return(input)
}
