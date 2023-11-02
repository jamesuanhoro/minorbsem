test_that("Random method (any case) returns histogram for CFA", {
  mbsem_test_hist(fit_cfa)
})

test_that("Random method (any case) returns histogram for SEM", {
  mbsem_test_hist(fit_sem)
})
