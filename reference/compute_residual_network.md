# Compute residual network

Interpret the error correlations as a residual network model.

## Usage

``` r
compute_residual_network(object)
```

## Arguments

- object:

  (mbsem) An object of
  [`mbsem-class`](https://jamesuanhoro.github.io/minorbsem/reference/mbsem-class.md)
  returned by
  [`minorbsem`](https://jamesuanhoro.github.io/minorbsem/reference/minorbsem.md).

## Value

A data.frame containing posterior samples of the partial correlation
matrix.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- minorbsem("F1 =~ x1 + x2 + x3
                  F2 =~ x4 + x5 + x6
                  F3 =~ x7 + x8 + x9", HS)
res_net <- compute_residual_network(fit)
p_corr_df <- posterior::summarise_draws(res_net)
n_items <- sqrt(nrow(p_corr_df))
p_corr_mat <- matrix(p_corr_df$mean, n_items)
p_corr_mat
qgraph::qgraph(p_corr_mat, layout = "spring")

# Complete Gaussian graphical model via a unidimensional model
# with all loadings set to zero.
fit <- minorbsem(paste0("F =~ ", paste0("0 * x", 1:9, collapse = " + ")), HS)
res_net <- compute_residual_network(fit)
p_corr_df <- posterior::summarise_draws(res_net)
n_items <- sqrt(nrow(p_corr_df))
p_corr_mat <- matrix(p_corr_df$mean, n_items)
p_corr_mat
qgraph::qgraph(p_corr_mat, layout = "spring")
} # }
```
