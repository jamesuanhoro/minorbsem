
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
- [Additional examples](#additional-examples)
  - [Different methods to capture the influence of minor
    factors](#different-methods-to-capture-the-influence-of-minor-factors)
  - [Relax simple structure](#relax-simple-structure)
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

The package also allows you set priors on all substantive model
parameters directly. Importantly, prior distributions assume latent
variables have a total variance of 1, even in latent regression models.

### Permitted models and supported data types

The package is able to fit a variety of model configurations:

- CFA, allowing automatically estimated penalized cross-loadings
- Path models with latent and observed variables
  - Any observed variables in a structural model must be represented
    with a single-indicator latent variable with the error variance of
    the observed variable constrained to 0

The package is also able to **analyze correlation structures** using
methods in Archakov and Hansen (2021). This includes polychoric
correlation matrices as long as an asymptotic variance matrix is
provided. The relevant paper is under review at Structural Equation
Modeling.

The package does not support fitting multi-group or multilevel models.
See the [bayesianmasem](https://github.com/jamesuanhoro/bayesianmasem)
package for multi-group factor analysis via meta-analysis methods.

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

## Additional examples

### Different methods to capture the influence of minor factors

Default method above is `method = "normal"` assuming standardized
residual covariances are on average 0 and vary from 0 in continuous
fashion.

``` r
# Fit same model as above but use global-local prior to estimate
# minor factor influences
fit_gdp <- minorbsem(syntax_1, HS, method = "GDP")
plot_residuals(fit_gdp)
parameter_hist(fit_gdp)
parameter_trace(fit_gdp)

# Ignoring minor factor influences
fit_none <- minorbsem(syntax_1, HS, method = "none")
parameter_hist(fit_none)
parameter_trace(fit_none)

# Error!!!: Plotting residuals will give an error message
# since minor factor influences are assumed null
plot_residuals(fit_none)
```

### Relax simple structure

``` r
fit_complex <- minorbsem(syntax_1, HS, simple_struc = FALSE)
```

There are other methods, see details section in `?minorbsem`.

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

<div id="ref-archakov_new_2021" class="csl-entry">

Archakov, Ilya, and Peter Reinhard Hansen. 2021. “A New Parametrization
of Correlation Matrices.” *Econometrica* 89 (4): 1699–1715.
<https://doi.org/10.3982/ECTA16910>.

</div>

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
