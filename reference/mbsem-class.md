# Class `"mbsem"`: For representing a model fitted with minorbsem

The `mbsem` class represents a fitted model. It contains summary
statistics for structural parameters, summary statistics for
standardized residual covariances, the data used to fit the model, and
the fitted Stan model.

## Objects from the Class

Objects can be created by calls to the
[`minorbsem`](https://jamesuanhoro.github.io/minorbsem/reference/minorbsem.md)
function.

## Slots

- `major_parameters`::

  Summary statistics for structural parameters

- `minor_factor_matrix`::

  Summary statistics for standardized residual covariances

- `data_list`::

  Data used to fit the model

- `priors`::

  Priors used to fit the model

- `stan_fit`::

  Fitted CmdStan model

- `version`::

  Package version used to fit model

## Methods

- logLik:

  `signature(object = "mbsem")`: Returns the casewise log-likelihood as
  long as the full data is available.

- fitted:

  `signature(object = "mbsem")`: Returns the posterior distribution of
  the model-implied covariance matrix as a \#(iterations) by \#(items
  ^ 2) matrix.

- residuals:

  `signature(object = "mbsem", standardized = TRUE)`: Returns the
  posterior distribution of residual covariances reflecting the
  influences of minor factors as a \#(iterations) by \#(items ^ 2)
  matrix. If `standardized = TRUE`, returns standardized residual
  covariances. If `standardized = FALSE`, returns UNstandardized
  residual covariances.

- show:

  `signature(object = "mbsem")`: Pretty printing of model results. See
  [`pretty_print_summary`](https://jamesuanhoro.github.io/minorbsem/reference/pretty_print_summary.md)
  for more printing options.

## See also

[`parameter_hist`](https://jamesuanhoro.github.io/minorbsem/reference/parameter_hist.md),
[`parameter_trace`](https://jamesuanhoro.github.io/minorbsem/reference/parameter_trace.md),
[`plot_residuals`](https://jamesuanhoro.github.io/minorbsem/reference/plot_residuals.md),
[`pretty_print_summary`](https://jamesuanhoro.github.io/minorbsem/reference/pretty_print_summary.md)

## Examples

``` r
  if (FALSE) { # \dontrun{
    fit_1 <- minorbsem("F1 =~ x1 + x2 + x3
                        F2 =~ x4 + x5 + x6
                        F3 =~ x7 + x8 + x9", HS)
    # Model-implied covariances
    mod_imp_cov <- fitted(fit_1)
    # Get average of model-implied covariance matrix
    matrix(colMeans(mod_imp_cov), nrow = fit_1@data_list$Ni)
  } # }
```
