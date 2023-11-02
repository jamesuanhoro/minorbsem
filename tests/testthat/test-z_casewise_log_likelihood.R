mbsem_test_ll_1 <- function(fit, method) {
  residuals <- sample(c(TRUE, FALSE), 1)

  if (method_hash(method) == 90) {
    err_msg <- paste0(
      "Cannot compute casewise log-likelihood when ",
      "method = \"WB\"; set method to \"WB-cond\"."
    )
    testthat::expect_error(
      casewise_log_likelihood(
        fit,
        include_residuals = residuals,
        use_armadillo = !residuals
      ),
      err_msg
    )
  } else {
    if (method_hash(method) >= 90 && isTRUE(residuals)) {
      warn_msg <- paste0(
        "include_residuals = TRUE is ignored when ",
        "minorbsem method == \"none\", \"WB\", \"WB-cond\", \"WW\"."
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
    testthat::expect_equal(nrow(ll_mat), 500 * 1)
    testthat::expect_true(!is.na(sum(ll_mat)))
    testthat::expect_true(sum(abs(ll_mat)) > 0)
  }
}

mbsem_test_ll_2 <- function(fit, method) {
  residuals <- sample(c(TRUE, FALSE), 1)

  if (method_hash(method) == 90) {
    err_msg <- paste0(
      "Cannot compute casewise log-likelihood when ",
      "method = \"WB\"; set method to \"WB-cond\"."
    )
    testthat::expect_error(
      casewise_log_likelihood(
        fit,
        include_residuals = residuals,
        use_armadillo = !residuals
      ),
      err_msg
    )
  } else {
    if (method_hash(method) >= 90 && isTRUE(residuals)) {
      warn_msg <- paste0(
        "include_residuals = TRUE is ignored when ",
        "minorbsem method == \"none\", \"WB\", \"WB-cond\", \"WW\"."
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
}

if (fit_cfa@data_list$has_data == 0) {
  cfa_dat <- dat_cov("HS", data_must = TRUE)
  fit_cfa <- minorbsem(
    model_cfa_syntax,
    data = cfa_dat$dat, sample_cov = cfa_dat$cov, sample_nobs = cfa_dat$nobs,
    orthogonal = orthogonal_cfa,
    simple_struc = sample(c(TRUE, FALSE), 1),
    warmup = 500, sampling = 500, chains = 1,
    method = method_cfa, refresh = 0, show_messages = FALSE
  )
}

if (fit_sem@data_list$has_data == 0) {
  sem_dat <- dat_cov("PD", data_must = TRUE)
  fit_sem <- minorbsem(
    model_sem_syntax,
    data = sem_dat$dat, sample_cov = sem_dat$cov, sample_nobs = sem_dat$nobs,
    orthogonal = orthogonal_sem,
    warmup = 500, sampling = 500, chains = 1,
    method = method_sem, refresh = 0, show_messages = FALSE
  )
}

test_that("Random method returns log-likelihood for CFA", {
  mbsem_test_ll_1(fit_cfa, method_cfa)
})

test_that("Random method returns log-likelihood for SEM", {
  mbsem_test_ll_1(fit_sem, method_sem)
})

test_that("CFA: Different LL methods are equal", {
  mbsem_test_ll_2(fit_cfa, method_cfa)
})

test_that("SEM: Different LL methods are equal", {
  mbsem_test_ll_2(fit_sem, method_sem)
})
