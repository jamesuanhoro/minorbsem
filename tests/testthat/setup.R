mbsem_test_hist <- function(fit) {
  # do default
  testthat::expect_error(
    gg <- parameter_hist(fit),
    NA
  )
  testthat::expect_true(inherits(gg, "ggplot"))
  # do all
  testthat::expect_error(
    gg <- parameter_hist(fit, param_type = "all"),
    NA
  )
  testthat::expect_true(inherits(gg, "ggplot"))
}

mbsem_test_trace <- function(fit) {
  # do default
  testthat::expect_error(
    gg <- parameter_trace(fit),
    NA
  )
  # do all
  testthat::expect_error(
    gg <- parameter_trace(fit, param_type = "all"),
    NA
  )
  testthat::expect_true(inherits(gg, "ggplot"))
}

mbsem_test_plot_residuals <- function(fit, method) {
  if (method_hash(method) >= 90) {
    testthat::expect_error(
      gg <- plot_residuals(fit, type = sample(c("matrix", "range"), 1)),
      paste0(
        "There are no residuals to plot when ",
        "method == \"none\", \"WB\", \"WB-cond\", \"WW\"."
      )
    )
  } else {
    testthat::expect_error(
      # do default
      gg <- plot_residuals(fit),
      NA
    )
    testthat::expect_true(inherits(gg, "ggplot"))
    testthat::expect_error(
      # do other
      gg <- plot_residuals(fit, type = "range"),
      NA
    )
    testthat::expect_true(inherits(gg, "ggplot"))
  }
}

mbsem_test_pp_shared <- function(print_out, method, simple = TRUE) {
  testthat::expect_true(grepl(method, print_out, ignore.case = TRUE))
  testthat::expect_true(grepl("Goodness of fit", print_out, ignore.case = TRUE))
  testthat::expect_true(grepl("RMSE", print_out, ignore.case = TRUE))
  testthat::expect_true(grepl("Factor loadings", print_out, ignore.case = TRUE))
  testthat::expect_true(grepl("PPP", print_out, ignore.case = TRUE))
  if (isFALSE(simple)) {
    testthat::expect_true(grepl("Location", print_out, ignore.case = TRUE))
    testthat::expect_true(grepl("Dispersion", print_out, ignore.case = TRUE))
    testthat::expect_true(grepl("convergence", print_out, ignore.case = TRUE))
    testthat::expect_true(grepl("median", print_out, ignore.case = TRUE))
    testthat::expect_true(grepl("mad", print_out, ignore.case = TRUE))
  }
}

dat_cov <- function(dat = "HS", data_must = FALSE) {
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

  if (isTRUE(data_must)) {
    input$data <- dat
  }

  return(input)
}

# Doing CFA here then sharing it across tests, instead
# of repeatedly fitting the model
method_cfa <- random_method_selection()
model_cfa_syntaxes <- c(
  "F1 =~ x1 + x2 + x3\n F2 =~ x4 + x5 + x6\n F3 =~ x7 + x8 + x9",
  "F1 =~ x1 + x2 + x3 + x9\n F2 =~ x4 + x5 + x6\n F3 =~ x7 + x8 + x9",
  "F1 =~ x1 + x2 + x3 + x9\n F2 =~ x4 + x5 + x6\n F3 =~ x7 + x8 + x9
    x2 ~~ x7\n  x3 ~~ x5"
)
syntax_cfa_idx <- sample(length(model_cfa_syntaxes), 1)
model_cfa_syntax <- model_cfa_syntaxes[syntax_cfa_idx]
cfa_dat <- dat_cov("HS")
orthogonal_cfa <- sample(c(TRUE, FALSE), 1)
expect_error(fit_cfa <- minorbsem(
  model_cfa_syntax,
  data = cfa_dat$dat, sample_cov = cfa_dat$cov, sample_nobs = cfa_dat$nobs,
  orthogonal = orthogonal_cfa,
  simple_struc = sample(c(TRUE, FALSE), 1),
  warmup = 500, sampling = 500, chains = 1,
  method = method_cfa, refresh = 0, show_messages = FALSE
), NA)

# Doing SEM here then sharing across tests
method_sem <- random_method_selection()
while (method_hash(method_cfa) >= 90 && method_hash(method_sem) >= 90) {
  method_sem <- random_method_selection()
}
model_sem_syntaxes <- c(
  "ind60 =~ x1 + x2 + x3\n dem60 =~ y1 + y2 + y3 + y4
  dem65 =~ y5 + y6 + y7 + y8\n dem60 ~ ind60\n dem65 ~ ind60 + dem60",
  "ind60 =~ x1 + x2 + x3\n dem60 =~ y1 + y2 + y3 + y4
  dem65 =~ y5 + y6 + y7 + y8\n dem65 ~ ind60 + dem60",
  "ind60 =~ x1 + x2 + x3\n dem60 =~ y1 + y2 + y3 + y4
  dem65 =~ y5 + y6 + y7 + y8\n dem65 ~ ind60 + dem60
  y1 ~~ y5\n  y4 ~~ y8"
)
syntax_sem_idx <- sample(length(model_sem_syntaxes), 1)
model_sem_syntax <- model_sem_syntaxes[syntax_sem_idx]
sem_dat <- dat_cov("PD")
while (!is.null(sem_dat$dat) && !is.null(cfa_dat$dat)) {
  sem_dat <- dat_cov("PD")
}
orthogonal_sem <- !orthogonal_cfa
expect_error(fit_sem <- minorbsem(
  model_sem_syntax,
  data = sem_dat$dat, sample_cov = sem_dat$cov, sample_nobs = sem_dat$nobs,
  orthogonal = orthogonal_sem,
  warmup = 500, sampling = 500, chains = 1,
  method = method_sem, refresh = 0, show_messages = FALSE
), NA)
