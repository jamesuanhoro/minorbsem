test_that("Random method (any case) works for CFA residual plot", {
  mbsem_test_plot_residuals(fit_cfa, method_cfa)
})

test_that("Random method (any case) works for SEM residual plot", {
  mbsem_test_plot_residuals(fit_sem, method_sem)
})
