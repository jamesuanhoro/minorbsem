# A class for setting up priors.

A class for setting up priors.

## Slots

- `lkj_shape`:

  (positive real) The shape parameter of the LKJ-prior on the
  interfactor correlation matrix in confirmatory factor models.

- `ml_par`:

  (real) The location parameter of the normal prior on loadings.

- `sl_par`:

  (positive real) The scale parameter of the normal prior on loadings.

- `rs_par`:

  (positive real) The scale parameter of the Student-t(df = 3, loc = 0)
  prior on the residual standard deviations.

- `rc_par`:

  (positive real) The shape parameter of the Beta(rc_par, rc_par) prior
  on the residual error correlations.

- `sc_par`:

  (positive real) The scale parameter of the normal prior on
  coefficients.

- `rm_par`:

  (positive real) The scale parameter of the normal prior on the tau /
  CRMR parameter.
