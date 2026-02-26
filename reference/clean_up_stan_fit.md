# Clean up Stan fit helper function

A function that cleans up model returned by Stan

## Usage

``` r
clean_up_stan_fit(stan_fit, data_list, interval = 0.9, priors)
```

## Arguments

- stan_fit:

  Stan fit

- data_list:

  Data list object passed to Stan

- interval:

  (real in (0, 1)) Credible interval to return.

- priors:

  An object of
  [`mbsempriors-class`](https://jamesuanhoro.github.io/minorbsem/reference/mbsempriors-class.md).
  See
  [`new_mbsempriors`](https://jamesuanhoro.github.io/minorbsem/reference/new_mbsempriors.md)
  for more information.

## Value

An object of
[`mbsem-class`](https://jamesuanhoro.github.io/minorbsem/reference/mbsem-class.md)
