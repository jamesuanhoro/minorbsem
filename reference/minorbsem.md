# Fit Bayesian SEMs with minor factors assumed

The major function to fit models assuming the influence of minor factors
(Uanhoro 2023) .

## Usage

``` r
minorbsem(
  model = NULL,
  data = NULL,
  sample_cov = NULL,
  sample_nobs = NULL,
  data_list = NULL,
  method = "normal",
  orthogonal = FALSE,
  simple_struc = TRUE,
  correlation = FALSE,
  centered = TRUE,
  seed = 12345,
  warmup = 1000,
  sampling = 1000,
  refresh = (warmup + sampling)/10,
  adapt_delta = 0.9,
  max_treedepth = 10,
  chains = 3,
  ncores = max(parallel::detectCores() - 2, 1),
  priors = new_mbsempriors(),
  show = TRUE,
  show_messages = TRUE,
  compute_ll = FALSE,
  acov_mat = NULL,
  ret_data_list = FALSE
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

- data_list:

  (list) A modified version of the data_list returned by minorbsem. Can
  be used to modify specific priors, see example below.

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

- seed:

  (positive integer) seed, set to obtain replicable results.

- warmup:

  (positive integer) The number of warmup iterations to run per chain.

- sampling:

  (positive integer) The number of post-warmup iterations to run per
  chain, retained for inference.

- refresh:

  (positive integer) How often to print the status of the sampler.

- adapt_delta:

  (real in (0, 1)) Increase to resolve divergent transitions.

- max_treedepth:

  (positive integer) Increase to resolve problems with maximum tree
  depth.

- chains:

  (positive integer) The number of Markov chains to run.

- ncores:

  (positive integer) The number of chains to run in parallel.

- priors:

  An object of
  [`mbsempriors-class`](https://jamesuanhoro.github.io/minorbsem/reference/mbsempriors-class.md).
  See
  [`new_mbsempriors`](https://jamesuanhoro.github.io/minorbsem/reference/new_mbsempriors.md)
  for more information.

- show:

  (Logical) If TRUE, show table of results, if FALSE, do not show table
  of results. As an example, use FALSE for simulation studies.

- show_messages:

  (Logical) If TRUE, show messages from Stan sampler, if FALSE, hide
  messages.

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

- ret_data_list:

  (LOGICAL) If TRUE, returns the `data_list` and `prior` objects, see
  example. If FALSE (default), fits the model given user inputs.

## Value

An object of
[`mbsem-class`](https://jamesuanhoro.github.io/minorbsem/reference/mbsem-class.md)

## Details

CFAs assume standardized factors. Latent variable regression models
print results with standardized loadings.

There are different methods for estimating models in this package:

- `normal`: under belief that minor factor influences are on average
  zero with continuous deviations away from zero (Uanhoro 2023) .

- `lasso`: under belief that minor factor influences are largely zero
  with a small number of non-zero residual covariances.

- `logistic`: for similar belief as normal but more readily accomodates
  extreme outliers.

- `GDP`: to mimic a global-local approach, i.e. attempt to shrink near 0
  residual covariances to 0 with minimal shrinking for larger residual
  covariances (Armagan et al. 2013) .

- `WB`: to model the covariance matrix hierarchically under assumptions
  of adventitiuous error (Wu and Browne 2015) ; does NOT allow for
  computation of casewise log-likelihoods and LOO-CV.

- `WB-cond`: same as WB but estimates the "population covariance
  matrix", allowing for computation of casewise log-likelihoods and
  LOO-CV.

- `WW`: A variation on WB-cond, but assumes the population covariance
  matrix is Wishart as opposed to inverse-Wishart;

- `none`: if intending to ignore the influence of minor factors.

WB-cond and WW are equivalent for correlation structure analysis.

## References

Archakov I, Hansen PR (2021). “A New Parametrization of Correlation
Matrices.” *Econometrica*, **89**(4), 1699–1715. ISSN 1468-0262,
[doi:10.3982/ECTA16910](https://doi.org/10.3982/ECTA16910) ,
<https://onlinelibrary.wiley.com/doi/abs/10.3982/ECTA16910>.  
  
Armagan A, Dunson DB, Lee J (2013). “Generalized double Pareto
shrinkage.” *Statistica Sinica*, **23**(1), 119–143. ISSN 1017-0405,
<https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3903426/>.  
  
Uanhoro JO (2023). “Modeling Misspecification as a Parameter in Bayesian
Structural Equation Models.” *Educational and Psychological
Measurement*, **0**(0), 00131644231165306. ISSN 0013-1644,
[doi:10.1177/00131644231165306](https://doi.org/10.1177/00131644231165306)
, <https://journals.sagepub.com/doi/abs/10.1177/00131644231165306>.  
  
Wu H, Browne MW (2015). “Quantifying Adventitious Error in a Covariance
Structure as a Random Effect.” *Psychometrika*, **80**(3), 571–600. ISSN
1860-0980,
[doi:10.1007/s11336-015-9451-3](https://doi.org/10.1007/s11336-015-9451-3)
.

## Examples

``` r
if (FALSE) { # \dontrun{
mod_cfa <- minorbsem(
  "# latent variable definitions
  F1 =~ x1 + x2 + x3
  F2 =~ x4 + x5 + x6
  F3 =~ x7 + x8 + x9", HS
)
new_pd <- PD
apply(PD, 2, sd) # first 8 variables have relatively large SDs
new_pd[, 1:8] <- new_pd[, 1:8] / 3 # move SDs closer to 1
mod_sem <- minorbsem(
  "# latent variable definitions
  ind60 =~ x1 + x2 + x3
  dem60 =~ y1 + y2 + y3 + y4
  dem65 =~ y5 + y6 + y7 + y8
  # latent regressions
  dem60 ~ ind60
  dem65 ~ ind60 + dem60", new_pd
)
mod_dl <- minorbsem(
  "# latent variable definitions
  F1 =~ x1 + x2 + x3
  F2 =~ x4 + x5 + x6
  F3 =~ x7 + x8 + x9", HS,
  ret_data_list = TRUE
)
mod_dl$load_est # prior mean for loadings
mod_dl$load_se # prior sd for loadings
mod_dl$load_se[9, 3] <- .75 # set prior SD for F3 =~ x9 to .75
# fit model with updated prior
mod <- minorbsem(data_list = mod_dl)
} # }
```
