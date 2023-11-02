test_that("Random method (any case) returns traceplot for CFA", {
  skip_if_not_installed("cmdstanr")
  mbsem_test_trace(fit_cfa)
})

test_that("Random method (any case) returns traceplot for SEM", {
  skip_if_not_installed("cmdstanr")
  mbsem_test_trace(fit_sem)
})
