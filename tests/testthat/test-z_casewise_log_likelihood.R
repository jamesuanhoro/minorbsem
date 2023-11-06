mbsem_test_ll_1 <- function(fit) {
  if (fit@data_list$ret_ll == 0) {
    err_msg <- paste0(
      "No log-likelihood returned. ",
      "You must set \"compute_ll = TRUE\" when calling ?minorbsem. ",
      "If you already did that, see the compute_ll parameter ",
      "description under ?minorbsem."
    )
    testthat::expect_error(
      logLik(fit),
      err_msg,
      fixed = TRUE
    )
  } else {
    testthat::expect_error(
      ll_mat <- logLik(fit),
      NA
    )
    testthat::expect_equal(ncol(ll_mat), nrow(fit@data_list$Y))
    testthat::expect_equal(nrow(ll_mat), 500 * 1)
    testthat::expect_true(!is.na(sum(ll_mat)))
    testthat::expect_true(sum(abs(ll_mat)) > 0)
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
    method = method_cfa, refresh = 0, show_messages = FALSE, compute_ll = TRUE
  )
}

if (fit_sem@data_list$has_data == 0) {
  sem_dat <- dat_cov("PD", data_must = TRUE)
  fit_sem <- minorbsem(
    model_sem_syntax,
    data = sem_dat$dat, sample_cov = sem_dat$cov, sample_nobs = sem_dat$nobs,
    orthogonal = orthogonal_sem,
    warmup = 500, sampling = 500, chains = 1,
    method = method_sem, refresh = 0, show_messages = FALSE, compute_ll = TRUE
  )
}

test_that("Random method returns log-likelihood for CFA", {
  mbsem_test_ll_1(fit_cfa)
})

test_that("Random method returns log-likelihood for SEM", {
  mbsem_test_ll_1(fit_sem)
})
