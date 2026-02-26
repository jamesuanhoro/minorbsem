# Set priors in package

Modify default priors in package.

## Usage

``` r
new_mbsempriors(
  lkj_shape = 2,
  ml_par = 0,
  sl_par = 1,
  rs_par = 1,
  rc_par = 2,
  sc_par = 0.5,
  rm_par = 0.15
)
```

## Arguments

- lkj_shape:

  (positive real) The shape parameter of the LKJ-prior on the
  interfactor correlation matrix in confirmatory factor models.

- ml_par:

  (real) The location parameter of the normal prior on loadings.

- sl_par:

  (positive real) The scale parameter of the normal prior on loadings.

- rs_par:

  (positive real) The scale parameter of the Student-t(df = 3, loc = 0)
  prior on the residual standard deviations.

- rc_par:

  (positive real) The shape parameter of the Beta(rc_par, rc_par) prior
  on the residual error correlations.

- sc_par:

  (positive real) The scale parameter of the normal prior on
  coefficients.

- rm_par:

  (positive real) The scale parameter of the normal prior on the tau /
  CRMR parameter.

## Value

An object of
[`mbsempriors-class`](https://jamesuanhoro.github.io/minorbsem/reference/mbsempriors-class.md)

## Examples

``` r
# Change LKJ shape parameter only
custom_priors <- new_mbsempriors(lkj_shape = 1.0)
if (FALSE) { # \dontrun{
minorbsem("F1 =~ x1 + x2 + x3
           F2 =~ x4 + x5 + x6
           F3 =~ x7 + x8 + x9", HS,
  priors = custom_priors
)
} # }
```
