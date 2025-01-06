test_that("PA produces CI", {
  mbsem_test_pa_ci(fit_pa, TRUE)
  mbsem_test_pa_ci(fit_pa, FALSE)
})
