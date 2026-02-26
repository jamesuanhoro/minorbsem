# Create major parameters helper function

A function that creates the table of major parameters

## Usage

``` r
create_major_params(stan_fit, data_list, interval = 0.9)
```

## Arguments

- stan_fit:

  Fitted Stan object

- data_list:

  Data list object passed to Stan

- interval:

  Confidence interval to select

## Value

Summary of posterior draws
