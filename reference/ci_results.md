# Conditional indepepence estimates for path analysis

Obtain conditional indepepence estimates for path analysis. Estimates
are standardized.

## Usage

``` r
ci_results(object, interval = 0.9, summarize = TRUE)
```

## Arguments

- object:

  (mbsem) An object of
  [`mbsem-class`](https://jamesuanhoro.github.io/minorbsem/reference/mbsem-class.md)
  returned by
  [`minorbsem`](https://jamesuanhoro.github.io/minorbsem/reference/minorbsem.md).

- interval:

  Confidence interval to select

- summarize:

  (LOGICAL) If TRUE (default): Return posterior summary as data.frame;
  If FALSE: Return posterior draws as data.frame.

## Value

Table of parameters related to conditional independence.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- minorbpa("x3 ~ x1 + x2\n x4 ~ x3 + x1", HS)
ci_table <- ci_results(fit)
} # }
```
