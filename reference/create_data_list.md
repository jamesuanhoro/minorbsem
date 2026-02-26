# Stan data helper function

A function that creates data list object passed to Stan

## Usage

``` r
create_data_list(
  lavaan_object = NULL,
  method = "normal",
  simple_struc = TRUE,
  correlation = correlation,
  priors = NULL,
  compute_ll = FALSE,
  partab = NULL,
  centered = TRUE,
  acov_mat = NULL,
  old_names = NULL
)
```

## Arguments

- lavaan_object:

  lavaan fit object of corresponding model

- method:

  (character) One of "normal", "lasso", "logistic", "GDP", "WB",
  "WB-cond", "WW", or "none". See details below.

- simple_struc:

  (LOGICAL) Only relevant for CFAs. If TRUE: assume simple structure; If
  FALSE: estimate all cross-loadings using generalized double Pareto
  priors.

- correlation:

  (LOGICAL) If TRUE: perform correlation structure analysis based on
  logarithm of a matrix transformation (Archakov and Hansen 2021) ; If
  FALSE (default): perform covariance structure analysis.

- priors:

  An object of
  [`mbsempriors-class`](https://jamesuanhoro.github.io/minorbsem/reference/mbsempriors-class.md).
  See
  [`new_mbsempriors`](https://jamesuanhoro.github.io/minorbsem/reference/new_mbsempriors.md)
  for more information.

- compute_ll:

  (Logical) If TRUE, compute log-likelihood, if FALSE, do not. This may
  be useful for cross-validation. This argument is ignored when: (i) the
  full dataset is not provided; (ii) the method is WB, use WB-cond
  instead.

- partab:

  lavaanify result of corresponding model

- centered:

  (LOGICAL) Only relevant for WB-cond and WW methods when
  `correlation = TRUE`. If TRUE (default): Use a centered
  parameterization; If FALSE: Use a non-centered parameterization.

- acov_mat:

  (Optional) Asymptotic variance matrix of lower triangular half
  (column-order) of the correlation matrix to be used for correlation
  structure analysis. This parameter is useful if importing polychoric
  or meta-analytic SEM pooled correlation matrix.

- old_names:

  (Optional) Variable name order of original correlation matrix, used to
  reorder acov_mat.

## Value

Data list object used in fitting Stan model
