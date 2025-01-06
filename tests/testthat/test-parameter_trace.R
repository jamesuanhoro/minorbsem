test_that("Random method (any case) returns traceplot for CFA", {
  mbsem_test_trace(fit_cfa)
})

test_that("Random method (any case) returns traceplot for SEM", {
  mbsem_test_trace(fit_sem)
})

test_that("Random method (any case) returns traceplot for PA", {
  mbsem_test_trace(fit_pa)
})
