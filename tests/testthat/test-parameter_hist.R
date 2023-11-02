test_that("Random method (any case) returns histogram for CFA", {
  skip_if_not_installed("cmdstanr")
  mbsem_test_hist(fit_cfa)
})

test_that("Random method (any case) returns histogram for SEM", {
  skip_if_not_installed("cmdstanr")
  mbsem_test_hist(fit_sem)
})
