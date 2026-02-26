# User input processing function

A function that processes.

## Usage

``` r
user_input_process(
  model = NULL,
  data = NULL,
  sample_cov = NULL,
  sample_nobs = NULL,
  method = "normal",
  orthogonal = FALSE,
  simple_struc = TRUE,
  correlation = FALSE,
  centered = TRUE,
  priors = new_mbsempriors(),
  compute_ll = FALSE,
  acov_mat = NULL,
  pa = FALSE
)
```

## Arguments

- model:

  A description of the user-specified model, lavaan syntax.

- data:

  An optional data frame containing the observed variables used in the
  model.

- sample_cov:

  (matrix) sample variance-covariance matrix. The rownames and/or
  colnames must contain the observed variable names.

- sample_nobs:

  (positive integer) Number of observations if the full data frame is
  missing and only sample covariance matrix is given.

- method:

  (character) One of "normal", "lasso", "logistic", "GDP", "WB",
  "WB-cond", "WW", or "none". See details below.

- orthogonal:

  (logical) constrain factors orthogonal, must be TRUE to fit bifactor
  models.

- simple_struc:

  (LOGICAL) Only relevant for CFAs. If TRUE: assume simple structure; If
  FALSE: estimate all cross-loadings using generalized double Pareto
  priors.

- correlation:

  (LOGICAL) If TRUE: perform correlation structure analysis based on
  logarithm of a matrix transformation (Archakov and Hansen 2021) ; If
  FALSE (default): perform covariance structure analysis.

- centered:

  (LOGICAL) Only relevant for WB-cond and WW methods when
  `correlation = TRUE`. If TRUE (default): Use a centered
  parameterization; If FALSE: Use a non-centered parameterization.

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

- acov_mat:

  (Optional) Asymptotic variance matrix of lower triangular half
  (column-order) of the correlation matrix to be used for correlation
  structure analysis. This parameter is useful if importing polychoric
  or meta-analytic SEM pooled correlation matrix.

- pa:

  (LOGICAL) If TRUE: Path-analytic model; If FALSE (default): Generic
  model.
