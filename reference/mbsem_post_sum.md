# Posterior summary helper function

A function that slightly modifies the default summary function in
posterior package

## Usage

``` r
mbsem_post_sum(stan_fit, variable, interval = 0.9, major = FALSE)
```

## Arguments

- stan_fit:

  Fitted Stan object

- variable:

  Variable(s) to search for in Stan

- interval:

  Confidence interval to select

- major:

  If TRUE, add some preamble for printing the major parameters table.

## Value

Summary of posterior draws
