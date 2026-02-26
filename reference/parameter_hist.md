# Histogram of parameter posterior distribution

Produce histograms of parameter posterior distribution, option to limit
plot by type of parameter.

## Usage

``` r
parameter_hist(object, subset = NULL, ...)
```

## Arguments

- object:

  (mbsem) An object of
  [`mbsem-class`](https://jamesuanhoro.github.io/minorbsem/reference/mbsem-class.md)
  returned by
  [`minorbsem`](https://jamesuanhoro.github.io/minorbsem/reference/minorbsem.md).

- subset:

  (character) Subset of parameters: NULL (Default) showing all estimated
  parameters; Any other response will be used as regular expressions to
  subset the parameters. It can be loading names or types of parameters.

- ...:

  additional arguments to relevant bayesplot function

## Value

bayesplot object

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- minorbsem("F1 =~ x1 + x2 + x3
                  F2 =~ x4 + x5 + x6
                  F3 =~ x7 + x8 + x9", HS)
parameter_hist(fit)
} # }
```
