# Visualize standardized residual covariances

Visualize distribution of standardized residual covariances

## Usage

``` r
plot_residuals(object, type = "matrix")
```

## Arguments

- object:

  (mbsem) An object of
  [`mbsem-class`](https://jamesuanhoro.github.io/minorbsem/reference/mbsem-class.md)
  returned by
  [`minorbsem`](https://jamesuanhoro.github.io/minorbsem/reference/minorbsem.md).

- type:

  (string) Either: "range" for lineranges (by default) or "matrix" for a
  matrix with point estimates in the lower half.

## Value

ggplot object

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- minorbsem("F1 =~ x1 + x2 + x3
                  F2 =~ x4 + x5 + x6
                  F3 =~ x7 + x8 + x9", HS)
plot_residuals(fit)
} # }
```
