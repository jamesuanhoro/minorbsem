# Pretty print model results

Nice printing of model results, optionally produces HTML document

## Usage

``` r
pretty_print_summary(
  object,
  interval = 0.9,
  digits = 3,
  simple = TRUE,
  save_html = NULL
)
```

## Arguments

- object:

  (mbsem) An object of
  [`mbsem-class`](https://jamesuanhoro.github.io/minorbsem/reference/mbsem-class.md)
  returned by
  [`minorbsem`](https://jamesuanhoro.github.io/minorbsem/reference/minorbsem.md).

- interval:

  (real in (0, 1)) Credible interval to return.

- digits:

  (positive integer) Number of decimal places to print in table

- simple:

  (Logical) TRUE to produce table with less information about
  parameters; FALSE: produces table with more information

- save_html:

  (string) Optional file name to save table as HTML

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- minorbsem("F1 =~ x1 + x2 + x3
                  F2 =~ x4 + x5 + x6
                  F3 =~ x7 + x8 + x9", HS)
pretty_print_summary(fit)
pretty_print_summary(fit, simple = FALSE)
} # }
```
