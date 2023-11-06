
# minorbsem

<!-- badges: start -->

[![Project Status: Active The project has reached a stable, usable state
and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
![GitHub](https://img.shields.io/github/license/jamesuanhoro/minorbsem)
[![Codecov test
coverage](https://codecov.io/gh/jamesuanhoro/minorbsem/branch/master/graph/badge.svg)](https://app.codecov.io/gh/jamesuanhoro/minorbsem?branch=master)
[![minimal R
version](https://img.shields.io/badge/R%3E%3D-3.4.0-6666ff.svg)](https://cran.r-project.org/)
![GitHub R package
version](https://img.shields.io/github/r-package/v/jamesuanhoro/minorbsem)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/minorbsem)](https://cran.r-project.org/package=minorbsem)
![GitHub last
commit](https://img.shields.io/github/last-commit/jamesuanhoro/minorbsem)
[![R-CMD-check](https://github.com/jamesuanhoro/minorbsem/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jamesuanhoro/minorbsem/actions/workflows/R-CMD-check.yaml)
[![minorbsem status
badge](https://jamesuanhoro.r-universe.dev/badges/minorbsem)](https://jamesuanhoro.r-universe.dev)
[![JOSS-paper](https://joss.theoj.org/papers/c0cd5b1a2d66bbf21fb00d237f646180/status.svg)](https://joss.theoj.org/papers/c0cd5b1a2d66bbf21fb00d237f646180)
<!-- badges: end -->

#### Table of Contents

- [Package overview](#package-overview)
  - [Goals](#goals)
  - [Permitted models and supported data
    types](#permitted-models-and-supported-data-types)
- [Installation](#installation)
- [A reasonably complete
  demonstration](#a-reasonably-complete-demonstration)
  - [Model comparisons](#model-comparisons)
- [Additional examples](#additional-examples)
  - [Different methods to capture the influence of minor
    factors](#different-methods-to-capture-the-influence-of-minor-factors)
- [Contributions are encouraged](#contributions-are-encouraged)
- [Citations](#citations)

## Package overview

Structural equation models (SEMs) rarely reject the null hypothesis that
there is no model misspecification. One explanation for this problem is
that covariance structures are influenced by major factors which we can
hypothesize about and minor factors which we cannot predict a-priori,
e.g. MacCallum and Tucker (1991).

### Goals

The goal of `minorbsem` is to facilitate fitting Bayesian SEMs that
estimate the influence of minor factors on the covariance matrix,
following the approach in Uanhoro (2023). Briefly, the method estimates
all residual covariances with priors that shrink them towards zero, and
the model returns the magnitude of the influence of minor factors.

### Permitted models and supported data types

The package only fits a limited number of model configurations:

- CFA, allowing cross-loadings (which may be automatically estimated),
  correlated errors terms, and fully oblique or orthogonal factors
  (useful for fitting bifactor models)
- SEMs allowing latent regressions (only), cross-loadings, and
  correlated error terms.

However, the package does not currently support fitting:

- MIMIC,
- multi-group models,
- multilevel models, or
- models with specially constrained parameters (e.g., setting two
  parameters equal).

All data are assumed multivariate normal, i.e. no binary, ordinal
models.

## Installation

`minorbsem` is hosted on GitHub, so we need the `remotes` package to
install it. We also need to install the `cmdstanr` package and CmdStan
in order to use Stan.

Instructions:

``` r
install.packages("remotes")  # install remotes

# next install cmdstanr and CmdStan:
install.packages(
  "cmdstanr",
  repos = c("https://mc-stan.org/r-packages/", getOption("repos"))
)
cmdstanr::check_cmdstan_toolchain(fix = TRUE)
cmdstanr::install_cmdstan()

# Then finally minorbsem:
remotes::install_github("jamesuanhoro/minorbsem")
```

## A reasonably complete demonstration

``` r
library(minorbsem)
# Basic Holzinger-Swineford model
syntax_1 <- "
F1 =~ x1 + x2 + x3
F2 =~ x4 + x5 + x6
F3 =~ x7 + x8 + x9"
# Expect a summary table output
fit_1 <- minorbsem(syntax_1, HS)

# Save output table to html file, see: ?pretty_print_summary for more options
pretty_print_summary(fit_1, save_html = "baseline_model.html")

# Histogram of parameters, see: ?parameter_hist for arguments
parameter_hist(fit_1)

# Traceplot of parameters, see: ?parameter_trace for arguments
parameter_trace(fit_1)

# Examine all standardized residual covariances
plot_residuals(fit_1)
plot_residuals(fit_1, type = "range")
```

### Model comparisons

Fit a second model that estimates all cross-loadings and shrinks them to
0 using a global-local prior (`simple_struc = FALSE`).

Then compare both models using the Leave-One-Out (LOO) method, using the
`loo` package: `install.packages("loo")`.

``` r
fit_2 <- minorbsem(syntax_1, HS, simple_struc = FALSE)

# Compute case wise log-likelihood
# exclude minor factor influences or else both models will fit data equally
ll_mat_1 <- casewise_log_likelihood(fit_1, include_residuals = FALSE)
ll_mat_2 <- casewise_log_likelihood(fit_2, include_residuals = FALSE)
chain_id <- posterior::as_draws_df(fit_1@stan_fit)$.chain

# loo for model 1
loo_1 <- loo::loo(
  ll_mat_1,
  r_eff = loo::relative_eff(ll_mat_1, chain_id = chain_id)
)
print(loo_1)

# loo for model 2
loo_2 <- loo::loo(
  ll_mat_2,
  r_eff = loo::relative_eff(ll_mat_2, chain_id = chain_id)
)
print(loo_2)

# Compare both models
print(loo::loo_compare(loo_1, loo_2), simplify = FALSE)
```

## Additional examples

### Different methods to capture the influence of minor factors

Default method above is `method = "normal"` assuming standardized
residual covariances are on average 0 and vary from 0 in continuous
fashion.

``` r
## Fit same model as above but use global-local prior to estimate
# minor factor influences
fit_gdp <- minorbsem(syntax_1, HS, method = "GDP")
plot_residuals(fit_gdp)
parameter_hist(fit_gdp)
parameter_trace(fit_gdp)

# Ignoring minor factor influences
fit_reg <- minorbsem(syntax_1, HS, method = "none")
parameter_hist(fit_reg)
parameter_trace(fit_reg)

# Error!!!: Plotting residuals will give an error message
# since minor factor influences are assumed null
plot_residuals(fit_reg)
```

There are other methods, see details section in `?minorbsem`.

Model comparison via LOO works as above in the [first
example](#model-comparisons).

## Contributions are encouraged

All users of R (or SEM) are invited to submit functions or ideas for
functions.

Feel free to:

- [open an issue](https://github.com/jamesuanhoro/minorbsem/issues/) to
  report a bug or to discuss recommendations;
- submit pull requests to recommend modifications or suggest
  improvements.

You can also email the package maintainer, James Uanhoro (James dot
Uanhoro at unt dot edu). Thank you for helping improve minorbsem :).

## Citations

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-maccallum_representing_1991" class="csl-entry">

MacCallum, Robert C., and Ledyard R. Tucker. 1991. “Representing Sources
of Error in the Common-Factor Model: Implications for Theory and
Practice.” *Psychological Bulletin* 109 (3): 502–11.
<https://doi.org/10.1037/0033-2909.109.3.502>.

</div>

<div id="ref-uanhoro_modeling_2023" class="csl-entry">

Uanhoro, James Ohisei. 2023. “Modeling Misspecification as a Parameter
in Bayesian Structural Equation Models.” *Educational and Psychological
Measurement* 0 (0): 00131644231165306.
<https://doi.org/10.1177/00131644231165306>.

</div>

</div>
